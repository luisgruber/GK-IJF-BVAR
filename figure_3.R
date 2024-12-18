library(stochvol)
library(bayesianVARs)
library(mvtnorm)
library(colorspace)

dire <- "figures/"
if(!dir.exists(dire)) dir.create(dire)

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
plotPHI <- function(object, summary = "median", colorbar = TRUE, ylabels = NULL, xlabels = NULL,
                    add_numbers = FALSE, zlim = NULL,main="",...){

  optionals <- list(...)

  PHI <- apply(object, 1:2, function(x) do.call(what = summary,
                                                args = list(x)))
  PHI_star <- t(apply(PHI, 1, rev)) # image orders differently

  if(add_numbers){
    alpha <- .5
  }else{
    alpha <- 1
  }

  if(summary %in% c("median", "mean")){
    colspace <- colorspace::diverge_hcl(1001, alpha = alpha, palette = "Blue-Red")
    if(is.null(zlim)){
      zlim <- c(-max(abs(PHI)),max(abs(PHI)))
    }
    colbreaks <- seq(zlim[1]*1.001, zlim[2]*1.001, len=1002)#[-1]
  }else if(summary %in% c("sd", "var", "IQR")){
    colspace <- colorspace::sequential_hcl(1001, alpha = alpha, rev = TRUE,
                                           palette = "Reds 2")#colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)#colorspace::sequential_hcl(5, h = 0, c = c(100, 0), l = 65, rev = TRUE, power = 1, alpha = alpha) #colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)
    if(is.null(zlim)){
      zlim <- c(0,max(abs(PHI)))
    }
  }

  M <- ncol(PHI)
  Kp <- nrow(PHI)
  p <- floor(Kp/M)
  if(colorbar){
    oldpar <- par(no.readonly = TRUE)

    if(!is.null(optionals$layoutmat)){
      if(optionals$layoutmat){
        on.exit(par(oldpar), add = TRUE)
        mat <- matrix(c(rep(1,9),2),nrow = 10)
        layout(mat)
      }
    }else{
      #oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)
      mat <- matrix(c(rep(1,9),2),nrow = 10)
      layout(mat)
    }
    mymar <- oldpar$mar
    mymar[1] <- 0.1
    par(mar=mymar)
  }
  image((PHI_star), zlim = zlim ,xaxt = "n", yaxt="n", col = colspace,
        bty="n", main = main)
  if((Kp-1)>M){
    abline(v = seq(M,Kp-1, M)/(Kp-1) - 0.5/(Kp-1), xpd = FALSE)
  }
  if(is.null(xlabels)){
    axis(3, labels = rownames(PHI), at = seq(0,1, length.out=Kp),
         tick = FALSE, las=2, cex.axis=.75)
  }else if(xlabels=="lags"){
    text(seq(M,Kp, M)/(Kp-1) - 0.5/(Kp-1) - M/(2*(Kp-1)), 1 + 0.5/(M-1) ,
         pos = 3, xpd = TRUE, labels = paste0("Lag ", 1:p))
  }
  if(is.null(ylabels)){
    axis(2, at = seq(1,0, length.out=M), labels = colnames(PHI),
         tick = FALSE, las = 2, cex.axis=.75)
  }

  if(add_numbers){
    text(rep(seq(0,1, length.out=Kp),M),
         sort(rep(seq(1,0, length.out=M),Kp), decreasing = TRUE),
         round(as.vector(PHI),3), cex=.75)
  }

  if(colorbar){
    mymar2 <- oldpar$mar
    mymar2[3] <- .1
    mymar2[1] <- 2.1
    par(mar=mymar2)
    plot.new()
    if(summary %in% c("median", "mean")){

      colbarlim <- length(which(colbreaks<min(PHI_star)))
      #colbarlim <- ceiling((min(PHI_star) - zlim[1]*1.001)/(diff(zlim)+2*zlim[2]/1000)*1001)
      adjustedcolspace <- colspace[-c(1:(colbarlim-1))]

      if(sign(min(PHI_star))!=sign(max(PHI_star))){
        llabels <- round(c(colbreaks[colbarlim],0,tail(colbreaks,1)),2)
        #llabels <- round(sort(c(range(PHI_star),0)),2)
        #
        at <- sort(c(0,1,(0-colbreaks[colbarlim])/(tail(colbreaks,1)-colbreaks[colbarlim])))
      }else{
        llabels <- round(range(PHI_star),2)
        at <- c(0,1)
      }
    }else if(summary %in% c("sd", "var", "IQR")){

      adjustedcolspace <- colspace
      llabels <- round(zlim,2)
      at <- c(0,1)

    }

    len <- length(adjustedcolspace)
    xxx <- seq(0,1,length=len+1)
    rect(xxx[1:len], rep(0, len), xxx[-1],
         rep(1, len), col = adjustedcolspace, border = adjustedcolspace)

    axis(1, labels = llabels, at = at,
         tick = FALSE, las=1, cex.axis=.9, line=-1.2)
  }

}

T <- 100
M <- 20
p <- 1

set.seed(1675)
PHI_in <- diag(M)*runif(M, -.8,.8)
PHI_in[1,20] <- runif(1, -.8,.8)
data <- data_sim(T, M, p,
                 PHI_setting = NULL,
                 L_setting = list(p=0.1, mu=0.001, sigma = 0.001),
                 SV = TRUE,
                 SV_setting = list(mu = -10, phi_range = c(0.85, 0.98), sigma_range= c(0.1,0.3)),
                 PHI_in = PHI_in)
PHI_true <- data$PHI
PHI_true <- array(PHI_true,c(M,M,1))
priorPHI_global <- specify_prior_phi(data = data$Y,lags = 1L,prior = "HS", global_grouping = "global")
priorPHI_semiglobal <- specify_prior_phi(data = data$Y,lags = 1L,prior = "HS", global_grouping = "olcl-lagwise")

priorL <- specify_prior_sigma(data = data$Y, type = "cholesky", cholesky_U_prior =  "HS")

burnin <- 1000L
draws <- 10000L
mod_global <- bvar(data$Y, lags = 1L, prior_intercept = FALSE,
                   draws = draws, burnin = burnin,
                   prior_phi = priorPHI_global,
                   prior_sigma = priorL)
mod_semiglobal <- bvar(data$Y, lags = 1L, prior_intercept = FALSE,
                       draws = draws, burnin = burnin,
                       prior_phi = priorPHI_semiglobal,
                       prior_sigma = priorL)

fileprefix <- paste0(dire, "figure_3")
pdf(paste0(fileprefix, "a.pdf"), height = 5, width = 5)
par(mar = c(1,1,1,1))
plotPHI(PHI_true, add_numbers = FALSE, colorbar = TRUE, #layoutmat = FALSE,
        main = "", xlabels = "", ylabels = "")
mtext("True", side = 1)
dev.off()

pdf(paste0(fileprefix, "b.pdf"), height = 5, width = 5)
par(mar = c(1,1,1,1))
plotPHI(mod_global$PHI, add_numbers = FALSE, colorbar = TRUE, #layoutmat = FALSE,
        main = "", xlabels = "", ylabels = "")
mtext("Posterior median", side = 1)
dev.off()

pdf(paste0(fileprefix, "c.pdf"), height = 5, width = 5)
par(mar = c(1,1,1,1))
plotPHI(mod_global$PHI,"IQR", add_numbers = FALSE, colorbar = TRUE, #layoutmat = FALSE,
        xlabels = "", ylabels = "")
mtext("Posterior interquartile range", side = 1)
dev.off()
pdf(paste0(fileprefix, "d.pdf"), height = 5, width = 5)
par(mar = c(1,1,1,1))
plotPHI(mod_semiglobal$PHI, add_numbers = FALSE, colorbar = TRUE, #layoutmat = FALSE,
        main = "",xlabels = "", ylabels = "")
mtext("Posterior median", side = 1)
dev.off()
pdf(paste0(fileprefix, "e.pdf"), height = 5, width = 5)
par(mar = c(1,1,1,1))
plotPHI(mod_semiglobal$PHI,"IQR", add_numbers = FALSE, colorbar = TRUE, #layoutmat = FALSE,
        xlabels = "", ylabels = "")
mtext("Posterior interquartile range", side = 1)
dev.off()
