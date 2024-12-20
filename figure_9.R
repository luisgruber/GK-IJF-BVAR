library(scatterplot3d)

dire <- "figures/"
if(!dir.exists(dire)) dir.create(dire)

set.seed(70000)
n <- 50000
PHI <- array(NA, c(n,3,3))
for (i in seq.int(n)) {
  B <- matrix(rnorm(3*3,0,100),3,3)
  L <- diag(3)
  L[upper.tri(L)] <- rnorm(3,0,100)
  L_i <- backsolve(L, diag(3))
  PHI[i,,] <- B %*% L_i
}
png(filename = paste0(dire,"figure_9.png"),
    width = 8, height = 5, units = "in", res = 300)
# jpeg(filename = paste0(dire,"str_reduced.jpg"),
#      width = 8, height = 5, units = "in", res = 300, quality = 100)
# pdf(paste0(dire,"str_reduced.pdf"),
#     width = 8, height = 5)
par(mfrow=c(2,3), mar = c(2.7,2.7,.5,.2) + 0.1, mgp = c(1.5,.5,0))
qqplot(qnorm(ppoints(500), 0,100),PHI[,1,1], xlab = "" ,pch=".", ylab = bquote(Phi[11]))
qqplot(qnorm(ppoints(500), 0,100),PHI[,1,2], xlab = "", pch=".", ylab = bquote(Phi[12]))
qqplot(qnorm(ppoints(500), 0,100),PHI[,1,3], xlab = "", pch=".", ylab = bquote(Phi[13]))
plot(PHI[,1,1], PHI[,2,1], pch=".", xlab = bquote(Phi[11]), ylab = bquote(Phi[21]))
plot(PHI[,1,1], PHI[,1,2], pch=".", xlab = bquote(Phi[11]), ylab = bquote(Phi[12]))
plot(PHI[,1,1], PHI[,1,3], pch=".", xlab = bquote(Phi[11]), ylab = bquote(Phi[13]))
dev.off()