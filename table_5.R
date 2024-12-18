library(dplyr)
library(tidyr)
library(readr)
library(knitr)

dire <- "tables/"
if(!dir.exists(dire)) dir.create(dire)

avg_LPL <- function(LPLs){
  numericalnormalizer <- max(LPLs) - 700
  log(mean(exp(LPLs - numericalnormalizer))) + numericalnormalizer
}

scores_main <- readRDS("results/table_4.rds")
scores_robustness <- readRDS("results/table_5.rds")
scores <- bind_rows(scores_main, scores_robustness)

scores_aggr <- scores |> group_by(model, priorU, est_end, p, score, ahead, fore_date) |>
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

scores_BF <- scores_aggr |> mutate(m_p = paste0(model,"-", p, "-", priorU)) |>
  select(!any_of(c("model", "p", "priorU"))) |>
  pivot_wider(names_from = m_p, values_from = LPL) |>
  mutate(across(-c(1:4), ~.x-`MP_LIT-4-HS`)) |>
  pivot_longer(-c(1:4), names_to = "m_p", values_to = "BF") |>
  separate_wider_delim(m_p, delim="-", names = c("model", "p", "priorU")) |>
  group_by(model, priorU, p, score, ahead) |>
  mutate(cum_BF = cumsum(BF)) |>
  ungroup()
HS_AT <- scores_BF |> filter(model %in% c("HS", "HS*"), priorU == "HS") |>
  mutate(priorU = "AT")
scores_BF <- bind_rows(scores_BF, HS_AT)


lags <- 2:4
myorder <- c("MP_LIT", "SHM" , "HS", "DL",  "NG", "R2D2", "SSVS", "HS*", "DL*_a", "NG*_a", "R2D2*_a")
captions <- "21"
options(knitr.kable.NA = "-")
scores_dat <- scores_BF |> filter(ahead==1L, fore_date == "2020-03-01",
                                  score=="LPL", p %in% lags,
                                  model %in% myorder) |>
  select(model, p, priorU, cum_BF) |> arrange(desc(priorU), p) |>
  mutate(p=paste0("p=",p)) |>
  pivot_wider(names_from = p, values_from = cum_BF) |>
  pivot_wider(names_from = priorU, values_from = 3:5) |>
  arrange(match(model, myorder))
knitr::kable((scores_dat), digits = 0,
             caption = "Sum of relative one-step-ahead log predictive likelihoods for several VARs equipped with different priors on u.") |>
  readr::write_lines(file = paste0(dire, "table_5.txt"))
