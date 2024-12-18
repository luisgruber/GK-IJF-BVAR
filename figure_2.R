library(scatterplot3d)

dire <- "figures/"
if(!dir.exists(dire)) dir.create(dire)
fileprefix <- paste0(dire, "figure_2")

n <- 3
r <- 100000

## HS vanilla
zeta <- abs(rt(r*3, 1))
theta <- matrix(abs(rt(r*3*n,1)), r*3,n)
phi_hs <- matrix(rnorm(r*3*n,0,theta*zeta),r*3,n)

png(filename = paste0(fileprefix, "a.png"),
    width = 5, height = 5, units = "in", res = 300)
# jpeg(filename = paste0(fileprefix, "a.jpg"),
#      width = 5, height = 5, units = "in", res = 300, quality = 100)
# pdf(filename = paste0(fileprefix, "a.pdf"),
#     width = 5, height = 5)
scatterplot3d::scatterplot3d(phi_hs, pch=".",
                             xlab = "",ylab = "",
                             zlab = "", tick.marks = FALSE,
                             color = scales::alpha(1, alpha = .3),
                             cex.main=.5,
                             box = FALSE,y.margin.add = 0 ,mar = c(.5,.5,-.1, 0) + 0.1,
                             xlim = c(-100,100),
                             ylim = c(-100,100),
                             zlim = c(-100,100))
mtext(expression(phi[1]), side = 1, adj = 0.3, line = -.5)
dims <- par("usr")
text(dims[1]+ 0.77*diff(dims[1:2]),dims[3]+ 0.15*diff(dims[3:4]),expression(phi[2]),srt=45)
text(dims[1]+ 0.015*diff(dims[1:2]), dims[3]+ 0.35*diff(dims[3:4]), expression(phi[3]))
dev.off()

## HS structured
zeta2 <- abs(rt(r*3, 1))
phi_hs2 <- phi_hs
phi_hs2[,3] <- rnorm(r*3,0,theta[,3]*zeta2)

png(filename = paste0(fileprefix, "d.png"),
    width = 5, height = 5, units = "in", res = 300)
# jpeg(filename = paste0(fileprefix, "d.jpg"),
#      width = 5, height = 5, units = "in", res = 300, quality = 100)
# pdf(paste0(fileprefix, "d.pdf"),
#     width = 5, height = 5)
scatterplot3d::scatterplot3d(phi_hs2, pch=".",
                             xlab = "",ylab = "",
                             zlab = "", tick.marks = FALSE,
                             color = scales::alpha(1, alpha = .3),
                             cex.main=.5,
                             box = FALSE,y.margin.add = 0 ,mar = c(.5,.5,-.1, 0) + 0.1,
                             xlim = c(-100,100),
                             ylim = c(-100,100),
                             zlim = c(-100,100))
mtext(expression(phi[1]), side = 1, adj = 0.3, line = -.5)
dims <- par("usr")
text(dims[1]+ 0.77*diff(dims[1:2]),dims[3]+ 0.15*diff(dims[3:4]),expression(phi[2]),srt=45)
text(dims[1]+ 0.015*diff(dims[1:2]), dims[3]+ 0.35*diff(dims[3:4]), expression(phi[3]))

dev.off()

## SSVS vanilla
p <-  rbeta(r, .1,.1)
gammas <- matrix(
  rbinom(r*n, 1, p),
  nrow = r
)
tau <- gammas
tau[which(gammas==0)] <- 0.01
tau[which(gammas==1)] <- 100
phi_ssvs <- matrix(
  rnorm(r*n, 0, as.vector(tau)),
  nrow = r
)

png(filename = paste0(fileprefix, "b.png"),
    width = 5, height = 5, units = "in", res = 300)
# jpeg(filename = paste0(fileprefix, "b.jpg"),
#      width = 5, height = 5, units = "in", res = 300, quality = 100)
# pdf(paste0(fileprefix, "b.pdf"),
#     width = 5, height = 5)
scatterplot3d::scatterplot3d(phi_ssvs, pch=".",
                             xlab = "",ylab = "",
                             zlab = "", tick.marks = FALSE,
                             color = scales::alpha(1, alpha = 0.3),
                             cex.main=.5,
                             box = FALSE,y.margin.add = 0 ,
                             mar = c(.5,.5,-.1, 0) + 0.1,
                             xlim = c(-200,200),
                             ylim = c(-200,200),
                             zlim = c(-200,200))
mtext(expression(phi[1]), side = 1, adj = 0.3, line = -.5)
dims <- par("usr")
text(dims[1]+ 0.75*diff(dims[1:2]),dims[3]+ 0.15*diff(dims[3:4]),expression(phi[2]),srt=45)
text(dims[1]+ 0.015*diff(dims[1:2]), dims[3]+ 0.35*diff(dims[3:4]), expression(phi[3]))
dev.off()

## SSVS structured
p2 <-  rbeta(r, .1,.1)
phi_ssvs2 <- phi_ssvs
gammas2 <- rbinom(r,1,p2)
tau2 <- gammas2*100
tau2[gammas2==0] <- .01
phi_ssvs2[,3] <- rnorm(r, 0, tau2)

png(filename = paste0(fileprefix, "e.png"),
    width = 5, height = 5, units = "in", res = 300)
# jpeg(filename = paste0(fileprefix, "e.jpg"),
#      width = 5, height = 5, units = "in", res = 300, quality = 100)
# pdf(paste0(fileprefix, "e.pdf"),
#     width = 5, height = 5)
scatterplot3d::scatterplot3d(phi_ssvs2, pch=".",
                             xlab = "",ylab = "",
                             zlab = "", tick.marks = FALSE,
                             color = scales::alpha(1, alpha = 0.3),
                             cex.main=.5,
                             box = FALSE,y.margin.add = 0 ,
                             mar = c(.5,.5,-.1, 0) + 0.1,
                             xlim = c(-200,200),
                             ylim = c(-200,200),
                             zlim = c(-200,200))
mtext(expression(phi[1]), side = 1, adj = 0.3, line = -.5)
dims <- par("usr")
text(dims[1]+ 0.75*diff(dims[1:2]),dims[3]+ 0.15*diff(dims[3:4]),expression(phi[2]),srt=45)
text(dims[1]+ 0.015*diff(dims[1:2]), dims[3]+ 0.35*diff(dims[3:4]), expression(phi[3]))
dev.off()

## SHM within
lambda1 <- rgamma(r, .01,.01)
lambda2 <- rgamma(r, .01,.01)
phi1 <- matrix(
  rnorm(n*r,0,sqrt(lambda1)),
  r,n
)

png(filename = paste0(fileprefix, "c.png"),
    width = 5, height = 5, units = "in", res = 300)
# jpeg(filename = paste0(fileprefix, "c.jpg"),
#      width = 5, height = 5, units = "in", res = 300, quality = 100)
# pdf(paste0(fileprefix, "c.pdf"),
#     width = 5, height = 5)
scatterplot3d::scatterplot3d(phi1, pch=".",
                             xlab = "",ylab = "",
                             zlab = "", tick.marks = FALSE,
                             color = scales::alpha(1, alpha = 0.2),
                             cex.main=.5,
                             box = FALSE,y.margin.add = 0 ,
                             mar = c(.5,.5,-.1, 0) + 0.1,
                             xlim = c(-2,2),
                             ylim = c(-2,2),
                             zlim = c(-2,2))
mtext(expression(phi[1]), side = 1, adj = 0.3, line = -.5)
dims <- par("usr")
text(dims[1]+ 0.77*diff(dims[1:2]),dims[3]+ 0.15*diff(dims[3:4]),expression(phi[2]),srt=45)
text(dims[1]+ 0.015*diff(dims[1:2]), dims[3]+ 0.35*diff(dims[3:4]), expression(phi[3]))
dev.off()

## SHM cross- vs ownlag
phi2 <- matrix(
  rnorm(n*r,0,sqrt(lambda2)),
  r,n
)

png(filename = paste0(fileprefix, "f.png"),
    width = 5, height = 5, units = "in", res = 300)
# jpeg(filename = paste0(fileprefix, "f.jpg"),
#      width = 5, height = 5, units = "in", res = 300, quality = 100)
# pdf(paste0(fileprefix, "f.pdf"),
#     width = 5, height = 5)
scatterplot3d::scatterplot3d(phi1[,1],phi1[,2],phi2[,1], pch=".",
                             xlab = "",ylab = "",
                             zlab = "", tick.marks = FALSE,
                             color = scales::alpha(1, alpha = 0.2),
                             cex.main=.5,
                             box = FALSE,y.margin.add = 0,
                             mar = c(.5,.5,-.1, 0) + 0.1,
                             xlim = c(-2,2),
                             ylim = c(-2,2),
                             zlim = c(-2,2))
mtext(expression(phi[1]), side = 1, adj = 0.3, line = -.5)
dims <- par("usr")
text(dims[1]+ 0.77*diff(dims[1:2]),dims[3]+ 0.15*diff(dims[3:4]),expression(phi[2]),srt=45)
text(dims[1]+ 0.015*diff(dims[1:2]), dims[3]+ 0.35*diff(dims[3:4]), expression(phi[3]))

dev.off()
