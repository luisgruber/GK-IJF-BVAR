library(dplyr)
library(tidyr)
library(tibble)

dire <- "figures/"
if(!dir.exists(dire)) dir.create(dire)

data <- readRDS("results/figure_5.rds")

r2d2_zeta <- (data$R2D2_a_star$a / (2*data$R2D2_a_star$xi)) |> t() |> as_tibble()|>
  rename_with(~ gsub("a","zeta",.x)) |>
  pivot_longer(everything(), names_to = "zeta", values_to="value")

ng_zeta <- (data$NG_a_star$a / (2*data$NG_a_star$xi)) |> t() |> as_tibble()|>
  rename_with(~ gsub("a","zeta",.x)) |>
  pivot_longer(everything(), names_to = "zeta", values_to="value")

hs_zeta <- (data$HS$zeta) |> t() |> as_tibble() |>
  pivot_longer(everything(), names_to = "zeta")

zetas <- (r2d2_zeta |> mutate(prior="R2D2") |>
            bind_rows((hs_zeta |> mutate(prior="HS")))) |>
  bind_rows(ng_zeta |> mutate(prior="NG")) |>
  mutate(prior=factor(prior))
zetas$prior <- factor(zetas$prior, levels = c("HS","NG","R2D2"))

pdf(paste0(dire, "figure_5a.pdf"), height = 16/1.3, width = 10/1.3)
modelsquote <- c(bquote(HS^"*"),
                 bquote(NG[a]^"*"),
                 bquote(R2D2[a]^"*"))
m <- 3
p <- 4
laymat <- matrix(c(1,sort(rep(2:(p+1),8)),p+2),ncol=1)
layout(laymat)
par(mar=c(0,3,0,1)+.1)
plot.new()
axis(side = 1, at = (1:m)/m - 1/(2*m), labels = parse(text=modelsquote), lwd=0,
     line = -2, cex.axis =1.5)
par( mar=c(0,3,0,1)+.1, mgp=c(1.8,0.5,0))
for(i in seq.int(p)){

  #par(mgp=c(1,0.5,0))
  #if(i==p) par(mar=c(3,3,0,1)+.1)
  boxplot(I(sqrt(value))~zeta+prior,
          data = zetas |> filter(zeta%in%paste0("zeta",c((i*2-1),i*2))),
          log="y", horizontal = FALSE,
          xlab = "", ylab=paste0("Lag: ", i), cex.lab=1.5,
          xaxt=ifelse(i!=p,"n","s"),
          names=  rep("",2*m)
  )
  #mtext(paste0("Lag: ", i))
  #par(mgp=c(3,1.5,0))
  #axis(side = 1, at = c(1.5,3.5), labels = c("HS","R2D2"), lwd=0)
  abline(v=seq(2.5,(2*(m-1))+.5,2))

}
par( mar=c(.5,3,0,1)+.1)
plot.new()
axis(side = 1, at = (1:(2*m))/(2*m) - 1/(2*2*m), labels = rep(c("ol","cl"),m), lwd=0,
     line = -2, cex.axis=1.5)
dev.off()

dl_a <- data$DL_a_star$a |> t() |> as_tibble() |> pivot_longer(everything(), names_to = "a") |> mutate(prior = "DL")
ng_a <- data$NG_a_star$a |> t() |> as_tibble() |> pivot_longer(everything(), names_to = "a") |> mutate(prior = "NG")
r2d2_a <- data$R2D2_a_star$a |> t() |> as_tibble() |> pivot_longer(everything(), names_to = "a") |> mutate(prior = "R2D2")
as <- bind_rows(dl_a, ng_a, r2d2_a)

as$prior <- factor(as$prior, levels = c("DL","NG","R2D2"))

pdf(paste0(dire, "figure_5b.pdf"), height = 16/1.3, width = 10/1.3)
modelsquote <- c(bquote(DL[a]^"*"),
                 bquote(NG[a]^"*"),
                 bquote(R2D2[a]^"*"))
m <- 3
p <- 4
laymat <- matrix(c(1,sort(rep(2:(p+1),8)),p+2),ncol=1)
layout(laymat)
par(mar=c(0,3,0,1)+.1)
plot.new()
axis(side = 1, at = (1:m)/m - 1/(2*m), labels = parse(text=modelsquote), lwd=0,
     line = -2, cex.axis =1.5)
par( mar=c(0,3,0,1)+.1, mgp=c(1.8,0.5,0))
for(i in seq.int(p)){

  #par(mgp=c(1,0.5,0))
  #if(i==p) par(mar=c(3,3,0,1)+.1)
  boxplot(value~a+prior,
          data = as |> filter(a%in%paste0("a",c((i*2-1),i*2))),
          horizontal = FALSE, log="y",
          xlab = "", ylab=paste0("Lag: ", i), cex.lab=1.5,
          xaxt=ifelse(i!=p,"n","s"),
          names=  rep("",2*m)
  )
  #mtext(paste0("Lag: ", i))
  #par(mgp=c(3,1.5,0))
  #axis(side = 1, at = c(1.5,3.5), labels = c("HS","R2D2"), lwd=0)
  abline(v=seq(2.5,(2*(m-1))+.5,2))

}
par( mar=c(.5,3,0,1)+.1)
plot.new()
axis(side = 1, at = (1:(2*m))/(2*m) - 1/(2*2*m), labels = rep(c("ol","cl"),m), lwd=0,
     line = -2, cex.axis=1.5)
dev.off()
