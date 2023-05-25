library(tidyverse)
options(tidyverse.quiet = TRUE)
library(readxl)
library("xlsx")

# Read data ---------------------------------------------------------------

raw_data <- readRDS(here::here("univariate_analysis/data/raw_data.RDS"))

# Prepare annotations -----------------------------------------------------

# unified table with sample data
sample_data <-
  raw_data$annotations %>% 
  filter(qualified_for_analysis %in% c("yes", "dev_epi")) %>% 
  mutate(
    sample_age_w = as.numeric(sample_age_w),
    seizures_any = if_else(is.na(seizures_any), "No", seizures_any),
    Dev_Age_Grouping_new = as.factor(Dev_Age_Grouping_new) %>% fct_relevel("0-10 weeks", "11-40 weeks", "> 40 weeks")
  )

comparisons <-
  sample_data %>%
  dplyr::select(sample_anon, sample_age_w, contains("FT05"),contains("FT06"),contains("FT07"),contains("FT08"),contains("FT09"),contains("FT10"), status, VGB_treatment, vgb_pretreatment)

comparisons <-
  comparisons %>%
  gather(comparison, condition, -sample_anon, -sample_age_w, -status, -VGB_treatment) %>%
  filter(!is.na(condition), !is.na(sample_anon)) %>%
  group_by(comparison) %>%
  filter(!duplicated(sample_anon)) %>%
  mutate(status = if_else(is.na(status), "never seiz - never VGB", status),
         VGB_treatment = if_else(is.na(VGB_treatment), "No", VGB_treatment),
         comparison_id = str_extract_all(comparison, "^\\(FT.{1,2}\\)"))  %>% 
  ungroup()

# Prepare data ------------------------------------------------------------

convertGeneIDs <- read_csv(here::here("univariate_analysis/data/convertGeneIDs.csv"))

to_save <- 
  pull(read_csv(here::here("univariate_analysis/data/names_of_selected_RNAs.csv")))

RNA_raw <- 
  raw_data$rna_log2_CPM_TMM_ZS_age_corr

RNA_raw <- RNA_raw %>% select(sample_anon, to_save)

# Script ------------------------------------------------------------------
data <- RNA_raw

final_results <-
  comparisons %>%
  split(., .$comparison) %>%
  lapply(function(x) {
   
    x$condition <- factor(x$condition)
    
    tryCatch({
      annotation <-
        x %>%
        # filter(VGB_treatment == "No") %>% 
        select(sample_anon, condition)
      
      df <-
        annotation %>% 
        left_join(data) %>% 
        mutate(condition = as.factor(condition))
      
      to_save <-
        df[-1] %>% 
        gather(key, value, -condition) %>% 
        group_by(key, condition) %>% 
        summarise(valid = sum(value > 1, na.rm = T),
                  n_half = n()/2) %>% 
        filter(valid > n_half) %>% 
        distinct(key)
      
      
      df <- df[, c("sample_anon", "condition", to_save[[1]])]
      cl <- parallel::makeForkCluster(10)
      results <-
        df %>% 
        gather(key, value, -sample_anon, -condition) %>% 
        split(., .$key) %>%
        parallel::parLapplyLB(cl, ., function(x) {
  
          results <-
            kruskal.test(value ~ condition, data = x) %>% 
            broom::tidy() %>% 
            transmute(`KW chi-sq` = statistic, p.value.KW = p.value) 
          browser()
          medians <-
            x %>% 
            group_by(condition) %>% 
            summarise(median = median(value, na.rm = T),
                      n = sum(!is.na(value)))
          
          ns <- paste(medians$n, collapse = " + ")
          
          dunntestres <-FSA::dunnTest(x$value, x$condition, method="bh")
          
          results <-
            bind_cols(
              results,
              valid = ns,
              dunntestres$res
            ) %>% 
            add_column(molecule_name = x$key[1], .before = 1)
          
          results_med <- data.frame(Comparison = c(paste(medians$condition[1], medians$condition[2], sep = " - "),
                                                   paste(medians$condition[1], medians$condition[3], sep = " - "),
                                                   paste(medians$condition[2], medians$condition[3], sep = " - ")),
                                    compared_n = c(paste(medians$n[1], medians$n[2], sep = " + "),
                                                   paste(medians$n[1], medians$n[3], sep = " + "),
                                                   paste(medians$n[2], medians$n[3], sep = " + ")),
                                    medians = c(paste(round(medians$median[1],1), round(medians$median[2],1), sep = " vs "),
                                                paste(round(medians$median[1],1), round(medians$median[3],1), sep = " vs "),
                                                paste(round(medians$median[2],1), round(medians$median[3],1), sep = " vs ")),
                                    log2FoldChange = NA,
                                    fold = NA)
          results_med$log2FoldChange[1] <- medians$median[1] - medians$median[2]
          results_med$log2FoldChange[2] <- medians$median[1] - medians$median[3]
          results_med$log2FoldChange[3] <- medians$median[2] - medians$median[3]
          results_med$fold <- 2 ^ results_med$log2FoldChange
          
          results <-
            results %>% left_join(results_med, by = "Comparison")
          
          return(results)
          
        }) %>% 
        bind_rows()
      
      parallel::stopCluster(cl)
      rm(cl)
     
      results <-
        results %>% 
        left_join(
          results %>%
            distinct(molecule_name, p.value.KW) %>% 
            mutate(p.value.fdr = p.adjust(p.value.KW, method = "BH"))
        ) %>% 
        select(molecule_name, p.value.KW, p.value.fdr, everything())
      
      name_comp <- str_replace(x$comparison[1], "/", "") 
      
      KW_results <-
        results[1:5] %>%
        distinct(molecule_name, p.value.KW, .keep_all = T)
      
      # write.xlsx(KW_results, paste0("publication_ready/RNASeq/KW_ZScore_cor_", name_comp, ".xlsx"), sheetName= "KW results", append = F)
      write_csv(KW_results, paste0("publication_ready/RNASeq/KW_results/KW_ZScore_cor_", name_comp, ".csv"))
      
      # write.xlsx(results, paste0("publication_ready/RNASeq/KW_ZScore_cor_", name_comp, "_Dunn.xlsx"), append = F)
      write_csv(results, paste0("publication_ready/RNASeq/KW_results/KW_ZScore_cor_", name_comp, "_Dunn.csv"))
      fold_sig <- subset(results, results$p.value.fdr < 0.05 & results$P.adj < 0.05 & (results$fold >= 1.5 | results$fold <= 0.66))
      
      fold_sig$comparison <- name_comp
      
      return(fold_sig)
      
    }, error = function(error_condition) {
      return(NULL)
    })
    
  })

results_to_save <-
  final_results %>% 
  bind_rows()

results_to_save <- 
  results_to_save %>%  
  left_join(convertGeneIDs %>% dplyr::select(ensembl_gene_id, hgnc_symbol), by = c("molecule_name" = "ensembl_gene_id")) %>% 
  mutate(gene = if_else(!is.na(hgnc_symbol), hgnc_symbol, molecule_name)) 

results_to_save %>% write_csv("publication_ready/RNASeq/KW_results/ZScore_results.csv")

results_to_save %>% 
  distinct(comparison, molecule_name) %>% 
  count(comparison) %>% 
  write_csv("publication_ready/RNASeq/KW_results/ZScore_summary.csv")
