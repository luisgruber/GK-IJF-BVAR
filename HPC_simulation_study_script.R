library("RhpcBLASctl")
blas_set_num_threads(1)
omp_set_num_threads(1)

running_variable <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(running_variable)) running_variable <- 5L

mysplit <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

set.seed(16428)
library(stochvol)
library(bayesianVARs)

# SIM ---------------------------------------------------------------------

data_sim <- function(T, M, p, PHI_setting, L_setting, SV=TRUE ,SV_setting=NULL, PHI_in =NULL, L_in = NULL) {

  companion <- function(PHI, p, intercept = FALSE){

    if(intercept) X <- PHI[-nrow(PHI),] else X <- PHI
    M <- ncol(X)

    X2 <- diag((p-1)*M)
    X3 <- matrix(0, M, (p-1)*M)

    companion <- cbind(X, rbind(X2,X3))

    return(companion)

  }

  if(is.null(PHI_setting$intercept)){
    K <- p*M
  }else  K <- p*M+1

  # Coefficients
  if(is.null(PHI_in)){
    stable <- FALSE
    while(stable == FALSE){
      PHI <- matrix(0, K,M)
      if(!is.null(PHI_setting$intercept)){

        intercept <- rbinom(M,1, PHI_setting$intercept$p)
        PHI[K,which(intercept==1)] <- rnorm(length(which(intercept==1)),PHI_setting$intercept$mu, PHI_setting$intercept$sigma)

      }

      for (lag in seq.int(p)) {
        diagonal <- rbinom(M, 1, PHI_setting$diagonal$p/lag^PHI_setting$lag_penalizer)
        diag(PHI[(M*(lag-1)+1):(M*lag), ])[which(diagonal==1)] <- rnorm(length(which(diagonal==1)), PHI_setting$diagonal$mu, PHI_setting$diagonal$sigma)

        upper <- rbinom((M^2-M)/2, 1, PHI_setting$offdiagonal$p/lag^PHI_setting$lag_penalizer)
        PHI[(M*(lag-1)+1):(M*lag),][upper.tri(PHI[(M*(lag-1)+1):(M*lag),])][which(upper==1)] <- rnorm(length(which(upper==1)), PHI_setting$offdiagonal$mu, PHI_setting$offdiagonal$sigma)

        lower <- rbinom((M^2-M)/2, 1, PHI_setting$offdiagonal$p/lag^PHI_setting$lag_penalizer)
        PHI[(M*(lag-1)+1):(M*lag),][lower.tri(PHI[(M*(lag-1)+1):(M*lag),])][which(lower==1)] <- rnorm(length(which(lower==1)),PHI_setting$offdiagonal$mu, PHI_setting$offdiagonal$sigma)
      }
      ## check companion
      PHI_comp <- companion(PHI, p, !is.null(PHI_setting$intercept))
      if(max(Mod(eigen(PHI_comp, only.values = TRUE)$values)) < 0.99){
        stable <- TRUE
      }

    }
  }else{
    if(nrow(PHI_in)==K | nrow(PHI_in)==(K+1)){
      if(nrow(PHI_in)==(K+1)){
        PHI_setting$intercept <- TRUE
      }
    }else{
      stop("dim(PHI_in) does not coincide with p and M!")
    }
    PHI <- PHI_in
  }


  # L
  if(is.null(L_in)){
    L <- diag(M)
    for(i in 2:(M)){
      #L[i,1:(i-1)] <- rnorm(i-1)
      ind <- rbinom(i-1,1,L_setting$p)
      L[which(ind==1), i] <- rnorm(1,L_setting$mu, L_setting$sigma)
    }
  }else{
    if(!identical(as.numeric((dim(L_in))), c(M,M))){
      stop("dim(L_in) does not coincide with M!")
    }
    if(any(L_in[lower.tri(L_in)]!=0)){
      stop("L_in must be an upper triangular matrix with ones on the diagonal!")
    }
    if(!all(diag(L_in)==1)){
      stop("L_in must be an upper triangular matrix with ones on the diagonal!")
    }
    L <- L_in
  }

  L_inv <- backsolve(L, diag(M))

  if(SV==TRUE){
    # simulate volas
    vol <- matrix(as.numeric(NA), T,M)
    for (m in seq.int(M)) {
      mu <- SV_setting$mu
      phi <- runif(1, SV_setting$phi_range[1], SV_setting$phi_range[2])
      sig <- runif(1, SV_setting$sigma_range[1], SV_setting$sigma_range[2])
      sv_sim <- svsim(T, mu = mu, phi = phi, sigma = sig)
      vol[,m] <- sv_sim$vol
    }
  }else if(SV==FALSE){

    D <- diag(abs(rnorm(M)))
    Sigma <- t(L_inv) %*% D %*% L_inv
  }

  # Y
  Y <- matrix(as.numeric(NA), T+p, M)
  if(SV==TRUE){
    Y[1:p,] <- rnorm(p*M,0,exp(SV_setting$mu))
  }else{
    Y[1:p,] <- rnorm(p*M,0,1)
  }

  for (t in (1+p):(T+p)) {

    if(SV==TRUE){
      Sigma_t <- t(L_inv) %*% diag(vol[t-p,]^2) %*% L_inv
    }else if(SV==FALSE){
      Sigma_t <- Sigma
    }

    if(is.null(PHI_setting$intercept)){

      meanvec <- as.vector(Y[(t-p):(t-1),])%*%PHI

    }else meanvec <- c(as.vector(Y[(t-p):(t-1),]), 1)%*%PHI

    Y[t,] <- mvtnorm::rmvnorm(1, mean = meanvec, sigma = Sigma_t)

  }

  return(list(PHI=PHI, L=L, Y=Y[-(1:p),]))

}

# Cluster -----------------------------------------------------------------

seeds <- sample(1:10000, 20)
observations <- c(50,100,200)
size <- c(20)
models <- c("HM","DL", "R2D2", "NG", "HS",  "MP_LIT", "SSVS2_f", "SSVS2" ,"DL_a",
            "DL_a_star", "R2D2_star", "NG_star", "HS_star", "SSVS2_star",
            "R2D2_a","R2D2_a_star", "NG_a", "NG_a_star")
scenarios <- c("sparse","medium", "dense")
lags <- 1
grid <- expand.grid(model=models,observations=observations, scenario=scenarios,
                    seed=seeds,  size=size, p=lags, intercept = FALSE)


p <- 1

burnin <- 1000
draws <- 10000

PHI_setting = list(sparse = list(intercept = NULL,
                                 diagonal = list(p=0.8, mu=0.3, sigma=0.3),
                                 offdiagonal = list(p=0.01, mu=0.3, sigma=0.3),
                                 lag_penalizer = 1),
                   medium = list(intercept = NULL,
                                 diagonal = list(p=0.8, mu=0.15, sigma=0.15),
                                 offdiagonal = list(p=0.1, mu=0.1, sigma=0.1),
                                 lag_penalizer = 1),
                   dense = list(intercept = NULL,
                                diagonal = list(p=0.8, mu=0.15, sigma=0.15),
                                offdiagonal = list(p=0.8, mu=0.01, sigma=0.01),
                                lag_penalizer = 1))
L_setting = list(sparse = list(p=0.01, mu=0.001, sigma = 0.001),
                 medium = list(p=0.1, mu=0.001, sigma = 0.001),
                 dense = list(p=0.8, mu=0.001, sigma = 0.001))
SV_setting = list(mu = -10, phi_range = c(0.85, 0.98), sigma_range= c(0.1,0.3))

tol <- 1e-20

for (run in running_variable) {

  M <- as.integer(grid[run, "size"])
  T <- as.integer(grid[run, "observations"])
  scenario <- as.character(grid[run, "scenario"])
  seed <- as.integer(grid[run, "seed"])
  model <- as.character(grid[run, "model"])
  p <- as.integer(grid[run, "p"])
  intercept <- as.logical(grid[run, "intercept"])

  if(intercept){

    PHI_setting[[scenario]]$intercept <- list(p=0.1, mu=0.1, sigma=0.1)

  }

  folder <- paste0("sim_study/T_", T, "/", scenario, "/", seed, "/", model, "/p_", p,
                   "/")
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  file <- paste0(folder, "all_scores.rds")

  if (file.exists(file)) {
    cat("Nothing to do!\n")
    next()
  }

  #Seed
  set.seed(seed = seed)

  #data sim
  data <- data_sim(T, M, p=1, PHI_setting[[scenario]], L_setting[["medium"]],
                   SV = TRUE, SV_setting = SV_setting)
  PHI <- data$PHI
  saveRDS(PHI,file = paste0(folder, "PHItrue.rds"))
  L <- data$L
  Y <- data$Y

  # Prior covariance-variance -----------------------------------------------

  prior_sigma <- specify_prior_sigma(data = Y, type = "cholesky",
                                     cholesky_U_prior = "HS")

  # Prior VAR coefficients --------------------------------------------------

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
  cpe <- p*M # coeffiecients per equation without intercept

  # discrete distribution with masses proportional to dexp(1)
  xx <- seq(1/(ncol(Y)*cpe), 1, len = 1000)
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
    sigma_sq <- bayesianVARs:::MP_sigma_sq(Y, 6)
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

  prior_phi <- specify_prior_phi(data = Y, lags = p, prior = priorPHI,
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
  # Estimate ---------------------------------------------------------------------

  mod <- NULL
  attributes(mod) <- list(class="try-error")
  try_count <- 0
  while(inherits(mod, "try-error") & try_count<=20) {
    try_count <- try_count + 1

    mod <- try(
      bvar(data = Y, lags = p, draws = draws, burnin = burnin,
           prior_intercept = intercept, prior_phi = prior_phi,
           prior_sigma = prior_sigma, sv_keep = "all"),
      silent = FALSE
    )
    if(!inherits(mod, "try-error")) {
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
      if( PHI_diagnostic>1 | SV_diagnostic>1 ){
        attributes(mod) <- list(class="try-error")
        set.seed(sample(1:10000,1))
      }
    }

  }
  if(!inherits(mod, "try-error")) {
    PHImean <- apply(mod$PHI, 1:2, mean)
    RMSE <- sqrt(mean((PHImean-PHI)^2))
    saveRDS(RMSE, file = paste0(folder, "RMSE.rds"))
  }
}

