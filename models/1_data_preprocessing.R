# Libraries ---------------------------------------------------------------
library(mlr)

# Data --------------------------------------------------------------------

combined_list <-
  list(
    TSC = TSC,
    proteome = proteome_first_tp,
    miRNA = miRNA,
    RNA = RNA,
    metabolite = metabolite,
    HMGB = HMGB,
    isomiR = isomiR,
    WGS_SNP = WGS_SNP
  )

# Feature selection -------------------------------------------------------

# target column must be a factor
patients <-
  patients %>%
  mutate(seizures_any = as.factor(seizures_any))


# Completed dataset -------------------------------------------------------

all_variables <-
reduce(combined_list, full_join, by = c( "code", "sample_code")) %>% 
  left_join(patients %>% select(code, seizures_any, sex), by = "code") %>% 
  select(-code, -sample_code)

all_completed <-
  all_variables %>% 
  na.omit()

save(all_variables, all_completed, file = here::here("models/full_dataset.RData"), compress = "gzip")

# joined datasets ---------------------------------------------------------

# joins target variable
combined_list <-
  combined_list %>%
  lapply(function(x) {
    x %>%
      left_join(patients %>%
                  select(code, seizures_any),
                by = "code") 
  })

trainData <-
  combined_list %>% 
  lapply(function(x){
    x %>% 
      select(-code, -sample_code)
  })

trainTasks <-
  trainData %>%
  lapply(function(x) {
    makeClassifTask(data = x, target = "seizures_any")
  })

filter_values <-
  trainTasks %>%
  lapply(function(x){
     generateFilterValuesData(x, method = "auc")
  })

filter_values <-
  filter_values %>%
  lapply(function(x){
    data.frame("feature_id" = x$data$name, "importance" = x$data$value) %>% 
      arrange(desc(importance))
  })

normalized_top <- mapply(function(x, y){
  top_cols <- 
    y %>% 
    filter(importance > 0.1) %>%
    top_n(30) %>% 
    pull(feature_id) %>% 
    as.character()
  x %>% 
    select(top_cols, seizures_any, code, sample_code)
}, combined_list, filter_values, SIMPLIFY = FALSE, USE.NAMES = TRUE)

normalized_top <- normalized_top[-which(names(normalized_top) == "HMGB")]

normalized_top$proteome <- 
  normalized_top$proteome[!(names(normalized_top$proteome) %in%
  c("P13645", "P35908", "P0DOY3", "Q92520", "P06312", "Q28107", "O95897", "P02654", "P01834"))]

features_in_analysis <-
  reduce(normalized_top, full_join, by = c("seizures_any", "code", "sample_code")) %>% 
  full_join(eCRF, by = c("code", "sample_code")) %>% ncol -2

no_of_combinations <-
  length(c(
    combn(features_in_analysis, 1, simplify = FALSE),
    combn(features_in_analysis, 2, simplify = FALSE),
    combn(features_in_analysis, 3, simplify = FALSE)
  ))

message("There will be ", no_of_combinations, " combinations")

# Combine datasets --------------------------------------------------------

combinations <-
  c(
    combn(length(normalized_top), 1, simplify = FALSE),
    combn(length(normalized_top), 2, simplify = FALSE),
    combn(length(normalized_top), 3, simplify = FALSE)
  )

dataset <-
  full_join(eCRF,
            reduce(normalized_top, full_join,
                   by = c("seizures_any", "code", "sample_code")),
            by = c("code", "sample_code")) %>% 
  select(-code, -sample_code) %>% 
  select(-seizures_any, everything())

raw_values <- raw_values %>% left_join(patients %>%  select(code, seizures_any))
selected_variables <- names(dataset)
colnames(raw_values) <- make.names(names(raw_values)) 
raw_values_2 <- raw_values[, c(selected_variables, "code", "sample_code", "age_weeks")]

save(dataset, raw_values,  file = here::here("models/seizures_normalized_c3.RData"), compress = "gzip")
