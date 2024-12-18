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

scores_BF <- scores_aggr |> mutate(m_p = paste0(model,"-", p, "-", priorU)) |>
  select(!any_of(c("model", "p", "priorU"))) |>
  pivot_wider(names_from = m_p, values_from = LPL) |>
  mutate(across(-c(1:4), ~.x-`MP_LIT-4-HS`)) |>
  pivot_longer(-c(1:4), names_to = "m_p", values_to = "BF") |>
  separate_wider_delim(m_p, delim="-", names = c("model", "p", "priorU")) |>
  group_by(model, priorU, p, score, ahead) |>
  mutate(cum_BF = cumsum(BF)) |>
  ungroup()

myLPLs <- c("LPL", "LPL_VoI")
for(myLPL in myLPLs){
  for(h in c(1L,4L)){
    ggplot(scores_BF |> filter(fore_date <= lubridate::ymd("2020-03-01"),
                               model %in% c("HS*", "R2D2*_a", "DL*_a","NG*_a", "SHM"),
                               ahead == h, score == {{myLPL}})) +
      geom_line(aes(x = fore_date, y = cum_BF, color = model)) +
      geom_hline(yintercept = 0) +
      facet_wrap(.~p, nrow = 1, labeller = "label_both") +
      geom_rect(data = nber_rec, aes(xmin = Start, xmax = End,
                                     ymin = -Inf, ymax = Inf),
                fill='grey80', alpha=0.7) +
      #("2009-03-01"), ("2015-12-01")
      geom_rect(data = zlb, aes(xmin = Start, xmax = End,
                                ymin = -Inf, ymax = Inf),
                fill='lightgreen', alpha=0.5) +
      theme_bw()+
      theme(plot.background = element_blank(),
            #panel.background = element_rect(fill="white", color = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "bottom",
            legend.key.width = unit(1,"in"),
            legend.direction = "horizontal") +
      scale_x_date(#minor_breaks = seq(ymd("1980-01-01"), ymd("2020-01-01"), by="10 year"),
        breaks = seq(ymd("1980-01-01"), ymd("2020-01-01"), by="5 year"),
        #date_labels = "%Y",
        labels = c(as.vector(sapply(seq(1985,2015,10), \(x) c("",x))),"")) +
      guides(color=guide_legend(nrow = 1))

    num <- if(myLPL == "LPL") 6L else 7L
    ind <- if(h==1L) "a" else "b"
    ggsave(paste0(dire, "figure_", num,ind, ".pdf"), width = 20/1.5, height = 5/1.5, dpi = 300)
  }
}

