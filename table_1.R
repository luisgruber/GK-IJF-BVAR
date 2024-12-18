library(knitr)
library(readr)

dire <- "tables/"
if(!dir.exists(dire)) dir.create(dire)

set.seed(6000)
hoyer <- function(x, n){

  sparseness <- (sqrt(n) - sum(abs(x))/sqrt(sum(x^2))) / (sqrt(n) -1)
  return(sparseness)

}
sim_gl <- function(n, p, ind_vec = NULL ,prior,
                   dl_a,
                   r2d2_a, r2d2_b,
                   ng_a, ng_b, ng_c,
                   ssvs_c1, ssvs_c2, ssvs_s1, ssvs_s2, ssvs_f_p){
  if(is.null(ind_vec)){
    ind_vec <- rep(1, p)
  }
  groups <- unique(ind_vec)

  sim <- matrix(0, n, p)
  for(g in seq_along(groups)){
    n_g <- length(which(ind_vec==groups[g]))
    if(prior=="hs"){
      zeta <- abs(rt(n,1))
      theta <- matrix(abs(rt(n*n_g,1)), n, n_g)
      V_p_sqrt <- as.vector(theta*zeta)
    }else if(prior=="dl"){
      theta <- matrix(rgamma(n*n_g, dl_a, 1/2), n, n_g)
      psi <- matrix(rexp(n*n_g, 1/2), n, n_g)
      V_p_sqrt <- as.vector(sqrt(psi)*theta)
    }else if(prior=="r2d2"){
      psi <- matrix(rexp(n*n_g, 1/2), n, n_g)
      theta <- matrix(rgamma(n*n_g, r2d2_a, r2d2_a/2), n, n_g)
      zeta <- 1/rgamma(n, r2d2_b, r2d2_a/2)

      V_p_sqrt <- as.vector(sqrt(psi*theta*zeta/2))
    }else if(prior=="ng"){
      theta <- matrix(rgamma(n*n_g, ng_a, ng_a/2), n, n_g)
      zeta <- 1/rgamma(n, ng_b, ng_c)

      V_p_sqrt <- as.vector(sqrt(theta*zeta))
    }else if(prior=="ssvs"){
      pr <- rbeta(n, ssvs_s1, ssvs_s2)
      taus <- matrix(rbinom(n*n_g, 1, pr), n, n_g)
      taus[taus==1L] <- ssvs_c2
      taus[taus==0L] <- ssvs_c1

      V_p_sqrt <- as.vector(taus)
    }else if(prior=="ssvs_f"){
      taus <- matrix(rbinom(n*n_g, 1, ssvs_f_p), n, n_g)
      taus[taus==1L] <- ssvs_c2
      taus[taus==0L] <- ssvs_c1

      V_p_sqrt <- as.vector(taus)
    }

    sim[,ind_vec==groups[g]] <- matrix(
      rnorm(n*n_g, 0, V_p_sqrt),
      n, n_g
    )
  }
  sim
}

table_1_data <- matrix(as.numeric(NA), 2, 6)
colnames(table_1_data) <- c("MP_LIT/SHM", "HS", "DL", "R2D2", "NG", "SSVS")
rownames(table_1_data) <- c("A", "B")

p <- 1000
n <- 1e4

myas <- c(0.001, 1)

dl1 <- sim_gl(n, p, prior = "dl", dl_a = myas[1])
dl1_hoyer <- apply(dl1, 1, hoyer, p)
table_1_data["B", "DL"] <- mean(dl1_hoyer)

dl2 <- sim_gl(n, p, prior = "dl", dl_a = myas[2])
dl2_hoyer <- apply(dl2, 1, hoyer, p)
table_1_data["A", "DL"] <- mean(dl2_hoyer)

r2d21 <- sim_gl(n, p, prior = "r2d2", r2d2_a = myas[1]/2, r2d2_b = .5)
r2d21_hoyer <- apply(r2d21, 1, hoyer, p)
table_1_data["B", "R2D2"] <- mean(r2d21_hoyer)

r2d22 <- sim_gl(n, p, prior = "r2d2", r2d2_a = myas[2]/2, r2d2_b = .5)
r2d22_hoyer <- apply(r2d22, 1, hoyer, p)
table_1_data["A", "R2D2"] <- mean(r2d22_hoyer)

ng1 <- sim_gl(n, p, prior = "ng", ng_a = myas[1]/2, ng_b = .5, ng_c = myas[1]/4)
ng1_hoyer <- apply(ng1, 1, hoyer, p)
table_1_data["B", "NG"] <- mean(ng1_hoyer)

ng2 <- sim_gl(n, p, prior = "ng", ng_a = myas[2]/2, ng_b = .5, ng_c = myas[2]/4)
ng2_hoyer <- apply(ng2, 1, hoyer, p)
table_1_data["A", "NG"] <- mean(ng2_hoyer)

hs <- sim_gl(n, p, prior = "hs")
hs_hoyer <- apply(hs, 1, hoyer, p)
table_1_data["A", "HS"] <- mean(hs_hoyer)

ssvs_a <- sim_gl(n, p, prior = "ssvs_f", ssvs_c1 = .010, ssvs_c2 = 100, ssvs_f_p = 0.5)
ssvs_a_hoyer <- apply(ssvs_a, 1, hoyer , p)
table_1_data["A", "SSVS"] <- mean(ssvs_a_hoyer)

ssvs2 <- sim_gl(n, p, prior = "ssvs_f", ssvs_c1 = .010, ssvs_c2 = 100, ssvs_f_p = 0.01)
ssvs2_hoyer <- apply(ssvs2, 1, hoyer , p)
table_1_data["B", "SSVS"] <- mean(ssvs2_hoyer)

mp_lit <- matrix(rnorm(p*n),n,p)
mp_lit_hoyer <- apply(mp_lit, 1, hoyer, p)
MP_LIT <- mean(mp_lit_hoyer)

shm <- matrix(rnorm(p*n, 0, rgamma(n,.1,.1)),n,p)
shm_hoyer <- apply(shm, 1, hoyer, p)
SHM <- mean(shm_hoyer)

# check that MP_LIT and SHM are approx. equal
check <- (round(MP_LIT, 2) == round(SHM, 2))

if(check){
  table_1_data["A", "MP_LIT/SHM"] <- MP_LIT
}

knitr::kable(table_1_data, digits = 2) |>
  write_lines(paste0(dire, "table_1.txt"))
