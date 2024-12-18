library("RhpcBLASctl")
blas_set_num_threads(1)
omp_set_num_threads(1)

set.seed(1234)
library(bayesianVARs)
library(lubridate)
library(xts)
library(coda)

running_variable <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(is.na(running_variable)) {
  running_variable <- 1L
  cluster <- FALSE
}else{
  cluster <- TRUE
  }

# get data
if(!cluster) {
  data <- readRDS(paste0("data/data_growth.RData"))

}else{
  data <- readRDS(paste0("data/data_growth.RData"))
}


# variables of interest
VoI <- c("GDPC1", "CPIAUCSL", "FEDFUNDS")

# grid
lagorder <- 1:5
seeds <- sample(1:100000,10)
est_end <- est_end <- seq(as.Date("1980-03-01"), as.Date("2020-03-01"), by = "quarter")#ymd(index(data$small["1980-03-01/2020-03-01"])) #2021-03-01 end of estimation windows
lags <- lagorder
sizes <- c(large="large")
models <- c("MP_LIT", "HM", "HS", "DL", "NG", "R2D2", "DL_a", "NG_a", "R2D2_a",
            "NG_star", "R2D2_star", "HS_star", "DL_a_star", "NG_a_star",
            "R2D2_a_star", "SSVS2_f", "SSVS2", "SSVS2_star")
grid_main <- expand.grid(seed=seeds, model=models, priorU="HS", est_end=est_end, p=lags)

lagorder_robustness <- 2:4
models_robustness <- c("DL","NG", "R2D2", "HS", "HM", "MP_LIT", "SSVS2_f",
                       "DL_a_star", "NG_a_star", "R2D2_a_star", "HS_star")
priorU <- c("HM", "FLAT")
grid0 <- expand.grid(seed=seeds, model=models_robustness, priorU=priorU, est_end=est_end,
                     p=lagorder_robustness, stringsAsFactors = FALSE)
gridDL <- expand.grid(seed=seeds, model=c("DL", "DL_a_star"), priorU="DL",
                      est_end=est_end, p=lagorder_robustness, stringsAsFactors = FALSE)
gridNG <- expand.grid(seed=seeds, model=c("NG","NG_a_star"), priorU="NG",
                      est_end=est_end, p=lagorder_robustness, stringsAsFactors = FALSE)
gridR2D2 <- expand.grid(seed=seeds, model=c("R2D2", "R2D2_a_star"), priorU="R2D2",
                        est_end=est_end, p=lagorder_robustness, stringsAsFactors = FALSE)
gridSSVS <- expand.grid(seed=seeds, model="SSVS2_f", priorU="SSVS2_f",
                        est_end=est_end, p=lagorder_robustness, stringsAsFactors = FALSE)
grid_robustness <- rbind(grid0, gridDL, gridNG, gridR2D2, gridSSVS)

grid <- rbind(grid_main, grid_robustness)

# Sampler settings
each <- 10 # number of predictive draws per posterior draw
if(cluster){
  burnin <- 5000
  draws <- 15000
}else{
  burnin <- 1000
  draws <- 2000
}

tol <- 1e-20 # controls that the prior variances do not get too small in order to avoid numerical problems

mysplit <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
rearrange <- numeric(each*draws)
for(r in seq.int(draws)){
  rearrange[((r-1)*each + 1):(r*each)] <- (0:(each-1)*draws) + r
}

for (run in running_variable) {
  # Get model specifications
  model <- as.character(grid$model[run])
  p <- as.integer(grid$p[run])
  seed <- real_seed <- as.integer(grid$seed[run])
  priorU <- as.character(grid$priorU[run])

  # Get data for estimation and predictive evaluation
  est_end <- ymd(grid$est_end[run]) # end of estimation period
  est_period <- paste0("/", est_end) # whole estimation period
  Y_est_raw <- data[["large"]][est_period] # data for estimation
  Y_est <- as.matrix(Y_est_raw)

  h <- 4 # forecast horizon
  eval_start <- est_end + months(3)
  eval_end <- est_end + months(12)
  #last available observation is "2021-06-01"; check maximum possible forecast length
  if(eval_end > ymd("2021-06-01")){
    eval_end <- ymd("2021-06-01")
    #check whether eval_start and eval_end are in the same year
    if(year(eval_end) > year(eval_start)){
      h1 <- 12- month(eval_start) # remaining month in year of eval_start
      h <- (h1 + month(eval_end))/3 + 1
    }else h <- (month(eval_end) - month(eval_start))/3 + 1 # transform difference in month to quarter
  }

  eval_dates <- paste0(eval_start,"/", eval_end)
  Y_obs_raw <- data[["large"]][eval_dates] # observed data for evaluation of h-step ahead predictions
  Y_obs <- as.matrix(Y_obs_raw)

  # create directory where everything will be stored
  folder <- paste0("empirical_application/", model, "/priorU_", priorU ,"/p_", p, "/", est_end, "/",
                   seed, "/")
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  filestostore <- paste0(folder, c("LPL", "LPL_VoI", "LPL_univariate", "MSFE", "MAFE"), ".rds")

  if(all(file.exists(filestostore))){
    cat("Nothing to do!\n")
    next()
  }

  # Prior covariance-variance -----------------------------------------------

  # create objects for each possible argument in specify_prior_phi() and specify_prior_sigma()
  # Prior on elements of triangular matrix U
  cholesky_U_prior <- if(priorU == "HS" | priorU == "SSVS" | priorU == "DL" | priorU == "NG" | priorU == "R2D2"){
    priorU
  }else if(priorU == "HM"){
    "HMP"
  }else if(priorU == "SSVS2_f"){
    "SSVS"
  }else if(priorU == "FLAT"){
    "normal"
  }

  # Hyperparameter settings
  nU <- (ncol(Y_est)^2-ncol(Y_est))/2 # number of free elements in U
  cholesky_DL_a <- "1/n"
  cholesky_DL_tol <- tol
  cholesky_R2D2_a <- 1/(2*nU)
  cholesky_R2D2_b <- 0.5
  cholesky_R2D2_tol <- tol^2
  cholesky_NG_a <- 1/(2*nU)
  cholesky_NG_b <- 0.5
  cholesky_NG_c <- 1/(4*nU)
  cholesky_NG_tol <- tol^2
  cholesky_SSVS_c0 <- if(priorU == "SSVS") 0.1 else if(priorU == "SSVS2_f") 0.001
  cholesky_SSVS_c1 <- if(priorU == "SSVS") 6 else if(priorU == "SSVS2_f") 1
  cholesky_SSVS_p <- .5
  cholesky_HMP_lambda_3 <- c(0.01,0.01)
  cholesky_normal_sds <- sqrt(10)

  # SV prior
  cholesky_heteroscedastic <- TRUE
  cholesky_priorhomoscedastic <-  as.numeric(NA)
  cholesky_priormu <- c(0,100)
  cholesky_priorphi <- c(20, 1.5)
  cholesky_priorsigma2 <- c(.5, .5)
  cholesky_priorh0 = "stationary"
  expert_sv_offset <- if((p == 5 & (model == "SSVS" | model == "SSVS_f" | model == "SSVS_star")) |
                         (p>1 & model == "FLAT")){
    1e-30
  }else{
    0
  }

  prior_sigma <- specify_prior_sigma(data = Y_est, type = "cholesky",
                                     cholesky_U_prior = cholesky_U_prior,
                                     cholesky_U_tol = 0,
                                     cholesky_heteroscedastic = cholesky_heteroscedastic,
                                     cholesky_priormu = cholesky_priormu,
                                     cholesky_priorphi = cholesky_priorphi,
                                     cholesky_priorsigma2 = cholesky_priorsigma2,
                                     cholesky_priorh0 = cholesky_priorh0,
                                     cholesky_priorhomoscedastic = cholesky_priorhomoscedastic,
                                     cholesky_DL_a = cholesky_DL_a,
                                     cholesky_DL_tol = cholesky_DL_tol,
                                     cholesky_R2D2_a = cholesky_R2D2_a,
                                     cholesky_R2D2_b = cholesky_R2D2_b,
                                     cholesky_R2D2_tol = cholesky_R2D2_tol,
                                     cholesky_NG_a = cholesky_NG_a,
                                     cholesky_NG_b = cholesky_NG_b,
                                     cholesky_NG_c = cholesky_NG_c,
                                     cholesky_NG_tol = cholesky_NG_tol,
                                     cholesky_SSVS_c0 = cholesky_SSVS_c0,
                                     cholesky_SSVS_c1 = cholesky_SSVS_c1,
                                     cholesky_SSVS_p = cholesky_SSVS_p,
                                     cholesky_HMP_lambda3 = cholesky_HMP_lambda_3,
                                     cholesky_normal_sds = cholesky_normal_sds,
                                     expert_sv_offset = expert_sv_offset)

  # Prior VAR coefficients --------------------------------------------------

  prior_intercept <- 1000

  priorPHI <- if(model == "HS" | model == "R2D2" | model == "NG" |
                 model == "DL" | model == "SSVS"){
    model
  }else if(model == "FLAT" | model == "MP_LIT"){
    "normal"
  }else if(model == "HM"){
    "HMP"
  }else if(model == "SSVS2" | model == "SSVS2_f" | model == "SSVS_star" |
           model == "SSVS2_star" | model == "SSVS_f"){
    "SSVS"
  }else if(model == "HS_star" | model == "HS_plus"){
    "HS"
  }else if(model == "NG_a" | model == "NG_star" | model == "NG_a_star"){
    "NG"
  }else if(model == "DL_a" | model == "DL_a_star" | model == "DL_plus" |
           model == "DL_plus_a_star" | model == "DL_plus_star"){
    "DL"
  }else if(model == "R2D2_a" | model == "R2D2_star" | model == "R2D2_a_star"){
    "R2D2"
  }

  priormean <- 0

  # Hyperparameter settings
  cpe <- p*ncol(Y_est) # coeffiecients per equation without intercept

  # discrete distribution with masses proportional to dexp(1)
  xx <- seq(1/(ncol(Y_est)*cpe), 1, len = 1000)
  dxx <- dexp(xx,1)
  dxxn <- dxx/sum(dxx)
  a_mat <- cbind(xx, dxxn)

  DL_a <- if(model == "DL" | model == "DL_plus" | model == "DL_plus_star"){
    "1/K"
  }else if(model == "DL_a" | model == "DL_a_star" | model == "DL_plus_a_star"){
    a_mat
  }
  DL_tol <- tol
  R2D2_a <- if(model == "R2D2" | model == "R2D2_star"){
    1/(2*cpe)
  }else if(model == "R2D2_a" | model =="R2D2_a_star"){
    a_mat
  }
  R2D2_b <- 0.5
  R2D2_tol <- 2*tol^2
  NG_a <- if(model == "NG" | model == "NG_star"){
    1/(2*cpe)
  }else if(model == "NG_a" | model == "NG_a_star"){
    a_mat
  }
  NG_b <- 0.5
  NG_c <- if(model == "NG" | model == "NG_star"){
    1/(4*cpe)
  }else if(model == "NG_a" | model == "NG_a_star"){
    "0.5*a"
  }
  NG_tol <- tol^2
  SSVS_c0 <- if(model == "SSVS" | model == "SSVS_star" | model == "SSVS_f"){
    0.1
  }else if(model == "SSVS2" | model == "SSVS2_star" | model == "SSVS2_f"){
    0.01
  }
  SSVS_c1 <- if(model == "SSVS" | model == "SSVS_star" | model == "SSVS_f"){
    10
  }else if(model == "SSVS2" | model == "SSVS2_star" | model == "SSVS2_f"){
    100
  }
  SSVS_semiautomatic <- TRUE
  SSVS_p <- if(model == "SSVS" | model == "SSVS_star" | model == "SSVS2" |
               model == "SSVS2_star"){
    c(1,1)
  }else if(model == "SSVS_f" | model == "SSVS2_f"){
    0.5
  }
  HMP_lambda1 <- c(0.01, 0.01)
  HMP_lambda2 <- c(0.01, 0.01)
  normal_sds <- if(model == "FLAT"){
    sqrt(10)
  }else if(model == "MP_LIT"){
    # original Minnesota prior as in Litterman (1986)
    # OLS variances of univariate AR(6) models for each variable
    sigma_sq <- bayesianVARs:::MP_sigma_sq(Y_est, 6)
    # prior variances (intercept = FALSE, because for comparability the prior
    # variance for the intercept is 1000 for all models)
    LIT_V_i <- bayesianVARs:::get_MP_V_prior(sigma_sq = sigma_sq, p=p,
                                             intercept=FALSE)
    sqrt(LIT_V_i)
  }
  global_grouping <- if(model == "SSVS_star" | model == "SSVS2_star" |
                        model == "HS_star" | model == "NG_star" |
                        model == "NG_a_star" | model == "DL_a_star" |
                        model == "DL_plus_star" | model == "DL_plus_a_star" |
                        model == "R2D2_star" | model == "R2D2_a_star"){
    "olcl-lagwise"
  }else{
    "global"
  }

  prior_phi <- specify_prior_phi(data = Y_est, lags = p, prior = priorPHI,
                                 priormean = priormean, PHI_tol = 0,
                                 DL_a = DL_a, DL_tol = DL_tol, R2D2_a = R2D2_a,
                                 R2D2_b = R2D2_b, R2D2_tol = R2D2_tol,
                                 NG_a = NG_a, NG_b = NG_b, NG_c = NG_c,
                                 NG_tol = NG_tol, SSVS_c0 = SSVS_c0,
                                 SSVS_c1 = SSVS_c1,
                                 SSVS_semiautomatic = SSVS_semiautomatic,
                                 SSVS_p = SSVS_p, HMP_lambda1 = HMP_lambda1,
                                 HMP_lambda2 = HMP_lambda2,
                                 normal_sds = normal_sds,
                                 global_grouping = global_grouping)

    # Estimate model ----------------------------------------------------------

  set.seed(seed = real_seed)
  success <- FALSE
  nr_tries <- 0
  while(!success & nr_tries<10){
    nr_tries <- nr_tries + 1
    mod <- try(
      bvar(data = Y_est, lags = p, draws = draws, burnin = burnin,
           prior_intercept = prior_intercept, prior_phi = prior_phi,
           prior_sigma = prior_sigma, sv_keep = "all"),
      silent = FALSE
    )
    cat("\n", gc(), "\n")
    if(inherits(mod, "try-error")){
      cat("\n", mod, "\n")
      real_seed <- sample(1:1000000,1)
      set.seed(real_seed)
    }else{

      cat("Finished sampling...running convergence checks!\n")
      ## Posterior con/divergence based on posterior of L2 norms
      # variance analysis

      PHI_mat <- matrix(mod$PHI, ncol = draws)
      PHI_norms <- sqrt(colSums(PHI_mat^2))
      # Split single chain into ten parts
      splitind <- mysplit(1:draws,10)
      varofsplittedchains <- unlist(lapply(splitind, FUN = function(x) var(PHI_norms[x])))
      PHI_diagnostic <- var(varofsplittedchains)/var(PHI_norms)

      sv_latent <- mod$logvar[dim(mod$logvar)[1],,]
      sv_latent_norm <- sqrt(colSums(sv_latent^2))
      # Split single chain into ten parts
      varofsplittedchainsSV <- unlist(lapply(splitind, FUN = function(x) var(sv_latent_norm[x])))
      SV_diagnostic <- var(varofsplittedchainsSV)/var(sv_latent_norm)

      # #(Geweke's convergence diagnostic)
      # cat("Convergence plots and info file!\n")
      # # geweke ('equality of means')
      # geweke_PHI_full <- coda::geweke.diag(PHI_norms, frac1 = 1/3, frac2 = 1/3)$z
      # info_df <- data.frame(real_seed=real_seed, nr_tries = nr_tries,
      #                       PHI_diagnostic = PHI_diagnostic,
      #                       SV_diagnostic = SV_diagnostic,
      #                       gk_PHI_full = geweke_PHI_full)

      if(PHI_diagnostic>1 | SV_diagnostic>1){
        cat("no convergence...rerunning the sampler!\n")
        real_seed <- sample(1:1000000,1)
        set.seed(real_seed)
      }else{
        cat("Chain succesfully converged!\n")

        cat("Start predicting...\n")
        pred <- tryCatch(bayesianVARs:::predict.bayesianVARs_bvar(mod, ahead = 1:h, each = each,
                                                                  stable = FALSE,
                                                                  simulate_predictive = TRUE,
                                                                  LPL = TRUE, LPL_VoI = VoI,
                                                                  Y_obs = Y_obs),
                         error = function(e) e)
        cat("\n", gc(), "\n")
        pred$predictions <- pred$predictions[,,rearrange]

        if(!inherits(pred, "error")){
          cat("Checking for outliers...\n")
          ## check for extreme outliers
          outliermat <- matrix(as.logical(NA),h, draws*each)
          PRED_diagnostic <- rep(as.numeric(NA), h)
          for(kk in seq.int(h)){
            dists <- sqrt(colSums(pred$predictions[kk,,]^2))
            outlier <- outliermat[kk,] <- dists > quantile(dists, .75) + 5e03*IQR(dists)
            if(sum(outlier>0)){
              # in case of extreme outliers, remove those predictions and recalculate LPLs
              pred$predictions[kk,,outlier] <- as.numeric(NA)

              numericalnormalizerfull <- max(pred$LPL_draws[kk,!outlier]) - 700
              pred$LPL[kk] <- log(mean(exp(pred$LPL_draws[kk,!outlier]-numericalnormalizerfull))) + numericalnormalizerfull

              numericalnormalizerVoI <- max(pred$LPL_sub_draws[kk,!outlier]) - 700
              pred$LPL_VoI[kk] <- log(mean(exp(pred$LPL_sub_draws[kk,!outlier]-numericalnormalizerVoI))) + numericalnormalizerVoI

              pred$LPL_univariate[kk,] <- log(apply(pred$PL_univariate_draws[kk,,!outlier],1,mean))
            }
            dists <- dists[!outlier]
            splitindPRED <- mysplit(1:length(dists),10)
            varofsplittedchainsPRED <- unlist(lapply(splitindPRED, FUN = function(x) var(dists[x])))
            PRED_diagnostic[kk] <- var(varofsplittedchainsPRED)/var(dists)
          }
          if(any(PRED_diagnostic>1)){
            real_seed <- sample(1:1000000,1)
            set.seed(real_seed)
          }else{
            cat("Succesfully converged!\n")
            success <- TRUE
          }
        }else{
          real_seed <- sample(1:1000000,1)
          set.seed(real_seed)
        }
      }
    }
  }

  if(!inherits(mod, "try-error")) {
    # save all draws for full samples
    if(eval_start == ymd("2020-03-01")){

      saveRDS(mod, file = paste0(folder, "mod.rds"))

    }

    if(!inherits(pred, "error")){
      # save log predictive likelihoods
      saveRDS(pred$LPL, paste0(folder,"LPL.rds"))
      saveRDS(pred$LPL_VoI, paste0(folder,"LPL_VoI.rds"))
      saveRDS(pred$LPL_univariate, paste0(folder,"LPL_univariate.rds"))
      cat("Succesfully saved predictive scores!\n")

    }
  }
}
