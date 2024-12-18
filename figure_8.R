library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(lubridate)

dire <- "figures/"
if(!dir.exists(dire)) dir.create(dire)

data <- readRDS("results/table_4.rds")

avg_LPL <- function(LPLs){
  numericalnormalizer <- max(LPLs) - 700
  log(mean(exp(LPLs - numericalnormalizer))) + numericalnormalizer
}

normalize_stable <- function(logx){
  logx_max <- max(logx)
  x_rel <- exp(logx - logx_max)
  x_rel/sum(x_rel)
}
dma <- function(LPL, alpha, selection=FALSE){
  names <- colnames(LPL)
  LPL <- as.matrix(LPL)
  T <- nrow(LPL)
  M <- ncol(LPL)
  # weight storage
  omega_t_pred <- matrix(as.numeric(NA), T, M)
  colnames(omega_t_pred) <- names
  # initialize weights
  omega_t_pred_temp <- rep(1/M,M)

  for(t in seq.int(T)){
    # discount and normalize latest weights
    omega_t_pred[t,] <- normalize_stable(alpha*log(omega_t_pred_temp))
    # update weights
    omega_t_pred_temp <- normalize_stable(log(omega_t_pred[t,])+LPL[t,] - max(LPL[t,]))
  }

  if(selection){
    DMA_LPL <- vapply(seq.int(T),
                      FUN = function(i){
                        LPL[i,which.max(omega_t_pred[i,])]
                      },
                      FUN.VALUE = numeric(1L)
    )
  }else{
    numerical_normalizer <- apply(LPL, 1, max) - 700
    DMA_LPL <- log(rowSums(omega_t_pred*exp(LPL-numerical_normalizer))) + numerical_normalizer
  }

  return(list(weights=omega_t_pred,
              LPL = LPL,
              DMA_LPL = DMA_LPL))
}

ldf <- function(LPL, alphas, layers = 2){
  n <- nrow(LPL)
  k <- ncol(LPL)
  m <- length(alphas)
  ldf_LPL <- array(as.numeric(NA), c(n,m,layers))
  LPL_tmp <- LPL
  for(l in seq.int(layers)){
    for(a in seq_along(alphas)){
      ldf_LPL[,a,l] <- dma(LPL = LPL_tmp, alpha = alphas[a], selection = FALSE)$DMA_LPL
    }
    LPL_tmp <- ldf_LPL[,,l]
  }
  ldf_LPL
}

scores_aggr <- data |> group_by(model, priorU, est_end, p, score, ahead, fore_date) |>
  summarise(LPL = avg_LPL(LPL)) |> ungroup()|>
  mutate(model = case_when(model %in% "HM"~"SHM",
                           model %in%"NG_star"~"NG*",
                           model %in%"R2D2_star"~"R2D2*",
                           model %in% "HS_star"~"HS*",
                           model %in%"DL_a_star"~"DL*_a",
                           model %in%"NG_a_star"~"NG*_a",
                           model %in%"R2D2_a_star"~"R2D2*_a",
                           model %in%"SSVS2_star"~"SSVS*_p",
                           model %in%"SSVS2_f"~"SSVS",
                           model %in%"SSVS2"~"SSVS_p",
                           .default = model),
         priorU = case_when(priorU %in% "HM"~"GS",
                            priorU %in% "SSVS2_f"~"AT",
                            priorU %in% "R2D2"~"AT",
                            priorU %in% "NG"~"AT",
                            priorU %in% "DL"~"AT",
                            .default = priorU))

# nber recessions
# nber_rec <- tis::nberDates() |> as_tibble() |>
#   mutate(Start = lubridate::ymd(Start),
#          End = lubridate::ymd(End)) |>
#   filter(Start >= "1980-06-01",
#          Start <= "2020-03-01")
nber_rec <- readRDS("data/nber_recessions.rds")

# zero lower bound
zlb <- data.frame(Start = lubridate::ymd("2009-03-01"), End = lubridate::ymd("2015-12-01"))

myLPLs <- c("LPL", "LPL_VoI")#"LPL_VoI" #"LPL"
for(myLPL in myLPLs){
  LPL_selection <- scores_aggr |>
    filter(ahead==1L, model%in%c("DL*_a","HS*","NG*_a","R2D2*_a", "SHM"), #p==2, #
           fore_date<lubridate::ymd("2020-06-01"),
           priorU=="HS",
           score == {{myLPL}}) |>
    mutate(model_p=paste0(model,"/",p)) |>
    select(fore_date, model_p, LPL) |>
    pivot_wider(names_from = model_p, values_from = LPL) |>
    arrange(fore_date)
  myDMA <- dma(LPL_selection[,-1], 0.99, selection = FALSE)
  alphas <- c(seq(.1,.9,.1),.95,.99,1)
  myldf <- ldf(LPL_selection[,-1], alphas = alphas, layers = 10L)
  # # colSums(myldf[,,4])
  #
  DMA_data <- bind_cols(date=LPL_selection$fore_date,myDMA$weights[,]) |> #-nrow(myDMA$weights)
    as_tibble() |>
    pivot_longer(-date, names_to = "model_p", values_to = "weights~(alpha==0.99)") |>
    mutate(model = str_split(model_p,"/") |> map_chr(~.x[1]),
           p = str_split(model_p,"/") |> map_chr(~.x[2])) |>
    select(date, model, p, `weights~(alpha==0.99)`) |>
    group_by(model, p)# |>
  #mutate(p=as_factor(p))

  LPL_data <- bind_cols(date=LPL_selection$fore_date,myDMA$LPL,"LDF/0"=myldf[,1,4]) |>
    as_tibble() |>  bind_cols("DMA/0"=myDMA$DMA_LPL) |>
    mutate(across(-date, ~.x-`DMA/0`)) |>
    pivot_longer(-date, names_to = "model_p", values_to = "LPL") |>
    mutate(model = str_split(model_p,"/") |> map_chr(~.x[1]),
           p = str_split(model_p,"/") |> map_chr(~.x[2])) |>
    select(date, model, p, LPL) |>
    group_by(model, p) |>
    mutate(#p=as_factor(p),
      "log~predictive~likelihoods"=cumsum(LPL))
  DMA_LPL_merge <- DMA_data |> full_join(LPL_data) |>
    pivot_longer(`weights~(alpha==0.99)`:`log~predictive~likelihoods`, values_to = "score", names_to = "DMA") |>
    mutate(score, DMA = factor(DMA, levels = c("weights~(alpha==0.99)","log~predictive~likelihoods")))



  ggplot(DMA_LPL_merge |> filter(model!="DMA", DMA!="LPL", model != "LDF")) +
    geom_line(aes(x = date, y = score, color = model, linetype=p), lwd=0.8) +
    facet_wrap(.~DMA, scales = "free_y", labeller = "label_parsed") +
    geom_hline(yintercept = 0) +
    geom_line(data = DMA_LPL_merge |> filter(model == "LDF", DMA=="log~predictive~likelihoods"),
              aes(date, score), color = "black", lwd=.8) +
    geom_text(aes(x = x,y = y,label = label), hjust = "inward",
              data.frame(x=lubridate::ymd("2020-03-01"), y = if(myLPL=="LPL") 15 else 5, label = "LDF", DMA= factor("log~predictive~likelihoods"))) +
    geom_rect(data = nber_rec, aes(xmin = Start, xmax = End,
                                   ymin = -Inf, ymax = Inf),
              fill='grey80', alpha=0.7) +
    #("2009-03-01"), ("2015-12-01")
    geom_rect(data = zlb, aes(xmin = Start, xmax = End,
                              ymin = -Inf, ymax = Inf),
              fill='lightgreen', alpha=0.5)+
    scale_x_date(#minor_breaks = seq(ymd("1980-01-01"), ymd("2020-01-01"), by="10 year"),
      breaks = seq(lubridate::ymd("1980-01-01"), lubridate::ymd("2020-01-01"), by="5 year"),
      #date_labels = "%Y",
      labels = c(as.vector(sapply(seq(1985,2015,10), \(x) c("",x))),"")) +
    theme_bw() +
    theme(plot.background = element_blank(),
          #         #panel.background = element_rect(fill="white", color = "black"),
          #         axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(),
          legend.key.width = unit(.75,"in"))

  ind <- if(myLPL=="LPL") "a" else "b"
  ggsave(paste0("figures/figure_8", ind,".pdf"), width = 12, height = 5, dpi = 300)
}

