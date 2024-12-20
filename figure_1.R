# In order to fully reproduce Figure 1 the Python library 'mpmath' is required.
# One can install and call an isolated python virtual environment via R package reticulate:

# install.packages("reticulate")
# reticulate::py_install("mpmath") # install 'mpmath'

# If PYTHON == FALSE, the density of R2D2 cannot be evaluated and the figure will be incomplete.
PYTHON <- FALSE # either TRUE or FALSE

library(gsl)
if(PYTHON){
  library(reticulate)
  mpmath <- reticulate::import('mpmath')
}

dire <- "figures/"
if(!dir.exists(dire)) dir.create(dire)

ddl <- function(x,a){

  d <- 1/(2^(0.5*(1+a)) * gamma(a)) * abs(x)^(0.5*(a-1)) * besselK(sqrt(2*abs(x)), 1-a)
  d

}
hs_cond <- function(x, zeta){
  -exp(x^2/(2*zeta)) * gsl::expint_Ei(-x^2/(2*zeta))/(sqrt(2*zeta)*pi^(3/2))
}
ng_cond <- function(x,a,zeta){
  2^(1/2-a) * (a/zeta)^(a/2+1/4)/(sqrt(pi)*gamma(a)) * besselK(abs(x)*sqrt(a/zeta), a-1/2)*(abs(x)^(a-1/2))
}

if(PYTHON){
  r2d2_cond <- function(x,a,z){

    d <- rep(0, length(x))
    for(i in seq_along(x)){
      d[i] <- as.numeric(as.character(
        mpmath$meijerg(r_to_py(list(list(),list())),
                       r_to_py(list(c(0.5,0,a-.5),list())),
                       (x[i]^2*a)/(4*z))$real
        #as.numeric(as.character(
        #mpmath$meijerg(r_to_py(list(c(1.5-a,1,0.5),list())),
        #              r_to_py(list(list(),list())),
        #             4*z/(x[i]^2*a))$real
      ))
    }
    d*sqrt(a/(2*z))/(sqrt(2*pi)*gamma(a))
  }
}

#hyperparameters
a_gamma <- .5
a_delta <- a_pi <- a_gamma/2
b <- 0.5
c <- a_delta/2

#x <- seq(-1,1,.001)
#x <- x[x!=0]
x <- c(seq(-1,-.001,.001),-.0001,-.00001,-.00000001)
x <- c(x,abs(rev(x)))
if(PYTHON){
  myr2d2 <- r2d2_cond(x, a_pi, 1)
}


filename <- paste0(dire, if(PYTHON) "figure_1" else "figure_1_incomplete")
pdf(paste0(filename, ".pdf"), width = 16/1.7, height = 7/1.7)
# png(paste0(filename, ".png"), width = 16/1.7, height = 7/1.7, res = 300, units = "in")
par(mfrow=c(1,3), mar = c(4.1,4.1,1.4,1.1), mgp = c(2.3,1,0))
plot(x , ng_cond(x, a_delta, 1),type = "l",col=3, log="y", ylim = c(0.05,(5)), xlab = bquote(phi),
     ylab = bquote("p"~(phi)), lwd = 2,
     xaxs="i", yaxs="i", bty="n", cex.lab=1.2)
lines(x, hs_cond(x, 1),type="l",col=1, lwd = 2)
if(PYTHON){
  lines(x,myr2d2  , col=2, lwd = 2)
}
lines(x, ddl(x, a_gamma),col=4, lwd = 2)
legend("topright", col = c(1,2,3,4), lty=1, legend = c("HS","R2D2","NG","DL"), bty = "n", lwd=2)

x_tail <- 5:20
plot(x_tail, ng_cond(x_tail, a_delta, 1), type = "l", log="y", col=3, xlab = bquote(phi),
     ylab = bquote("p"~(phi)), lwd = 2,ylim = c(1e-06,5e-03),
     xaxs="i", yaxs="i", bty="n", cex.lab=1.2)
lines(x_tail, hs_cond(x_tail, 1),type="l",col=1, lwd = 2)
if(PYTHON){
  lines(x_tail, r2d2_cond(x_tail, a_pi, 1), col=2, lwd = 2)
}
lines(x_tail, ddl(x_tail, a_gamma),col=4, lwd = 2)
legend("topright", col = c(1,2,3,4), lty=1, legend = c("HS","R2D2","NG","DL"), bty = "n", lwd=2)
hc <- function(x, df, log = FALSE){
  ret <- log(2) - (log(pi)+log(df*(1+(x/df)^2)))
  if(!log) ret <- exp(ret)
  ret
}
dgg <- function(x, a, d, p, log = FALSE){
  ret <- log(p) - d*log(a) + (d-1)*log(x) -
    lgamma(d/p) - (x/a)^p
  if(!log) ret <- exp(ret)
  ret
}
shape <- b
rate <- c
q <- 1/2
d <- shape/q
p <- 1/q
a <- (1/rate)^q
x_g <- seq(.01,5,.01)
plot(0, 0, type="n", xlim=c(0,5), ylim=c(-0.015,1.2),xlab = bquote(sqrt(zeta)),
     ylab = bquote("p"~(sqrt(zeta))), lwd = 2,
     xaxs="i", yaxs="i", bty="n",cex.lab=1.2)
lines(x_g, hc(x_g, 1), col=1, lwd = 2)
lines(x_g, dgg(1/x_g,a,d,p)/x_g^2, col = 2, lwd = 2)
lines(x_g, dgg(1/x_g,a,d,p)/x_g^2, col = 3, lwd = 2, lty=2)
#lines(c(1,1),c(0,1),col=4, lwd = 2)
points(1,0, pch=16, col=4, lwd=2)
legend("topright", col = c(1,2,3,4), lty=c(1,1,1,NA), pch = c(NA,NA,NA,16) ,legend = c("HS","R2D2","NG","DL"), bty = "n", lwd=2)
dev.off()
