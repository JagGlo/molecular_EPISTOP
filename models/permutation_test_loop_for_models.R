source("~/interim_report/reports/data_preprocessing.R") # preprocessed data & additional files

n <- 100
results <- list()

for (i in 1:n) {

  print(patients$seizures_any)
  
  patients$seizures_any <- sample(patients$seizures_any)
  
  sapply(paste0("/mnt/epistop/", dir(path="/mnt/epistop/", pattern="chunk_")), unlink)
  
  tryCatch({
    source(here::here("models/1_data_preprocessing.R"))
 
    results[[length(results)+1]] <- no_of_combinations
  
  print(filter_values$TSC)
  print(dataset)
  
  if(no_of_combinations < 1500000) {

    source(here::here("models/2_regression.R"))
  
    source(here::here("models/3_read_chunks.R"))

    write_csv(model_stats, paste0("~/interim_report/models/iteration_", i, ".csv"))
  }
  },
  error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  }
  )
  
}

write.csv(results %>% unlist(), file = here::here("models/Results.csv"), append = T)

