library(tidyverse)

all_csvs <- list.files("/mnt/epistop/", full.names = TRUE, pattern = ".csv")

model_stats <-
  pbapply::pblapply(all_csvs, read_csv,
                    col_types = cols(
    features = col_character(),
    n = col_double(),
    thr = col_double(),
    mcc.test.mean = col_double(),
    mmce.test.mean = col_double(),
    npv.test.mean = col_double(),
    ppv.test.mean = col_double()
  )) %>% 
  bind_rows() %>% 
  arrange(desc(mcc.test.mean))

good_models <-
  model_stats %>% 
  filter(n > 48) %>% 
  top_n(1000, wt = mcc.test.mean)

model_stats %>% 
  filter(!str_detect(features, "\\.y\\.n")) %>%
  arrange(desc(mcc.test.mean))

selected_models <-
  good_models %>% 
  top_n(100, )
