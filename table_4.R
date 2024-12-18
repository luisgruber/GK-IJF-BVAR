library(dplyr)
library(tidyr)
library(knitr)
library(readr)

data <- readRDS("results/table_4.rds")

dire <- "tables/"
if(!dir.exists(dire)) dir.create(dire)

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

scores_BF <- scores_aggr |> mutate(m_p = paste0(model,"-", p, "-", priorU)) |>
  select(!any_of(c("model", "p", "priorU"))) |>
  pivot_wider(names_from = m_p, values_from = LPL) |>
  mutate(across(-c(1:4), ~.x-`MP_LIT-4-HS`)) |>
  pivot_longer(-c(1:4), names_to = "m_p", values_to = "BF") |>
  separate_wider_delim(m_p, delim="-", names = c("model", "p", "priorU")) |>
  group_by(model, priorU, p, score, ahead) |>
  mutate(cum_BF = cumsum(BF)) |>
  ungroup()

myorder <- c("MP_LIT", "SHM" , "HS", "DL",  "NG", "R2D2", "DL_a", "NG_a",
             "R2D2_a", "NG*", "R2D2*", "HS*", "DL*_a", "NG*_a", "R2D2*_a",
             "SSVS", "SSVS_p", "SSVS*_p")
options(knitr.kable.NA = "-")
aheads <- c(1L,4L)
scores_kable <- vector("list", length(aheads))
names(scores_kable) <- paste0("h_", aheads)
for(h in seq_along(aheads)){
  scores_kable[[h]] <- scores_BF |> filter(ahead==aheads[h], fore_date == "2020-03-01",
                                           score=="LPL") |>
    select(model, p, priorU, cum_BF) |>
    mutate(p=paste0("p=",p)) |>
    pivot_wider(names_from = p, values_from = cum_BF) |>
    arrange(match(model, myorder)) |>
    select(-priorU)
}

readr::write_lines(
  knitr::kable(cbind(scores_kable$h_1,scores_kable$h_4[,-1]),linesep = "", digits = 0,
               caption = "Sum of one and four-steps-ahead log predictive likelihoods relative to MP LIT(4)."),
  file = paste0(dire, "table_4.txt")
)

