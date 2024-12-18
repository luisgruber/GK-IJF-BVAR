library(dplyr)
library(tidyr)
library(knitr)
library(readr)

dire <- "tables/"
if(!dir.exists(dire)) dir.create(dire)

# Read data
scores <- readRDS("results/table_2.rds")

# Order of models
myorder <- c("MP_LIT", "SHM" ,
             "DL", "DL_a","DL*_a",
             "HS","HS*",
             "NG", "NG_a", "NG*", "NG*_a",
             "R2D2",  "R2D2_a",  "R2D2*", "R2D2*_a",
             "SSVS", "SSVS_p", "SSVS*_p")

scores |> as_tibble() |>
  group_by(model, size, scenario, observations) |>
  summarise(RMSE = 100*median(score)) |> ungroup() |>
  mutate(so = paste0(observations, "_", scenario)) |>
  select(-c(scenario, observations)) |>
  pivot_wider(names_from = so, values_from = RMSE) |>
  mutate(model = as.character(model)) %>% arrange(model)|>
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
                           .default = model)) |>
  arrange(match(model, myorder)) |> select(-size) |>
  knitr::kable(digits = 2, caption = "100xRMSE") |> 
  readr::write_lines(paste0(dire, "table_2.txt"))
