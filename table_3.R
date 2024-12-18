library(readr)
library(knitr)
library(dplyr)
library(tidyr)

dire <- "tables/"
if(!dir.exists(dire)) dir.create(dire)

data <- readRDS("results/table_3.rds")

mymodels <- c("DL_a", "DL*_a", "SHM", "HS", "HS*", "MP_LIT", "NG_a",  "NG*_a",
              "R2D2_a", "R2D2*_a",  "SSVS", "SSVS*_p")

filen <- paste0(dire, "table_3.txt")
file.create(filen)
options(knitr.kable.NA = "-")
for(mod in mymodels){
  df <- data |> filter(model=={{mod}}) |>
    mutate(l = paste0("l=",l)) |> select(-a) |>
    pivot_wider(names_from = l, values_from = h) |>
    select(starts_with("l="))

    readr::write_lines(knitr::kable(df, digits = 2, caption = mod), filen, append = TRUE) #format = "latex", vline = "", toprule="\\hline\\hline", midrule="\\hline",

}
