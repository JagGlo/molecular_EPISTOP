library(tidyverse)
options(tidyverse.quiet = TRUE)
library(readxl)
library("xlsx")

# Prepare annotations -----------------------------------------------------

# unified table with sample data
# unified table with sample data
sample_data <-
  read_excel("publication_ready/data/EPISTOP-Annotations for final Tests (210221 - Age Filter).xlsx", 
             sheet = "TestGroups(prep)", na = "NA"
  ) %>%
  dplyr::rename(samples = `Matabolomics Sample #`) %>%
  filter(!is.na(samples)) %>% 
  filter(qualified_for_analysis %in% c("yes", "dev_epi")) %>% 
  mutate(
    sample_age_w = as.numeric(sample_age_w),
    seizures_any = if_else(is.na(seizures_any), "No", seizures_any),
    Dev_Age_Grouping = as.factor(Dev_Age_Grouping) %>% fct_relevel("0-5 weeks", "6-11 weeks", "12-21 weeks", "22-34 weeks", "35-102 weeks", "> 103 weeks"),
    Dev_Age_Grouping_new = as.factor(Dev_Age_Grouping_new) %>% fct_relevel("0-10 weeks", "11-40 weeks", "> 40 weeks")
  )

comparisons <-
  sample_data %>%
  dplyr::select(samples, sample_age_w, contains("FT05"),contains("FT06"),contains("FT07"),contains("FT08"),contains("FT09"),contains("FT10"), status, VGB_treatment, vgb_pretreatment)

comparisons <-
  comparisons %>%
  gather(comparison, condition, -samples, -sample_age_w, -status, -VGB_treatment, -vgb_pretreatment) %>%
  filter(!is.na(condition), !is.na(samples)) %>%
  group_by(comparison) %>%
  filter(!duplicated(samples)) %>%
  mutate(status = if_else(is.na(status), "never seiz - never VGB", status),
         VGB_treatment = if_else(is.na(VGB_treatment), "No", VGB_treatment),
         vgb_pretreatment = if_else(is.na(vgb_pretreatment), "No", vgb_pretreatment),
         comparison_id = str_extract_all(comparison, "^\\(FT.{1,2}\\)"))  %>% 
  ungroup()

# Prepare data ------------------------------------------------------------

# ZScore
data <-
  read_excel("publication_ready/data/Metabolomics Data before and after Correction.xlsx",
             sheet = "Batch VGB Age corrected") %>%
  gather('samples', 'value', -Metabolite) %>%
  spread(Metabolite, value) %>%
  mutate_if(is.character, as.numeric) %>%
  arrange(samples) %>% 
  left_join(sample_data %>% transmute(samples,`Proteomics File Name` = if_else(!is.na(`Proteomics File Name`), `Proteomics File Name`, `Metabolomics updated Sample ID`))) %>% 
  group_by(`Proteomics File Name`) %>% 
  summarise_all(mean, na.rm = T) 

data$samples[data$`Proteomics File Name` == "Ol-L_female_6m_Crtl"] <- 46
data$samples[data$`Proteomics File Name` == "S-A_female_24m_Crtl"] <- 48
data$samples[data$`Proteomics File Name` == "P-M_female_21m_1_Crtl"] <- 49
data$samples[data$`Proteomics File Name` == "Wo-An_1_Crtl"] <- 50

data <- data %>% select(-`Proteomics File Name`)
# from 332 to 306
# MM
# data <-
#   read_csv("publication_ready/data/Metabolomics_MM_corrected_data.csv")


# Script ------------------------------------------------------------------

final_results <-
  comparisons %>%
  split(., .$comparison) %>%
  lapply(function(x) {
    
    x$condition <- factor(x$condition)
   
    tryCatch({
      annotation <-
        x %>%
        # filter(VGB_treatment == "No") %>% 
        select(samples, condition)
      
      df <-
        annotation %>% 
        left_join(data) %>% 
        mutate(condition = as.factor(condition))
      
      results <-
        df %>% 
        gather(key, value, -samples, -condition) %>% 
        split(., .$key) %>%
        lapply(function(x) {
          
          x <- 
            x %>% 
            filter(!is.na(value), !is.na(condition))
          
          results <-
            kruskal.test(value ~ condition, data = x) %>% 
            broom::tidy() %>% 
            transmute(`KW chi-sq` = statistic, p.value.KW = p.value) 
          
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
      
      write.xlsx(KW_results, paste0("publication_ready/Metabolomics/KW_results/KW_ZScore_cor_", name_comp, ".xlsx"), sheetName= "KW results", append = F)
      
      write.xlsx(results, paste0("publication_ready/Metabolomics/KW_results/KW_ZScore_cor_", name_comp, ".xlsx"), sheetName= "Dunn results", append = T)
      
      
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

results_to_save %>% write.xlsx("publication_ready/Metabolomics/KW_results/ZScore_results.xlsx")

results_to_save %>% 
  distinct(comparison, molecule_name) %>% 
  count(comparison) %>% 
  write.xlsx("publication_ready/Metabolomics/KW_results/ZScore_summary.xlsx")

# Presenting the results --------------------------------------------------

# present_results <- function(results_csv, before_correction, data_after_correction, annotation_set){
#   
#   results_csv %>% 
#     separate(feature, "feature", sep = ";") %>% 
#     split(., .$comparison) %>% 
#     lapply(function(x) {
#       
#       comparison_id <-  str_extract(x$comparison[1], "^\\(FT\\d{1,2}\\)")
#       signif_molecules <- x$feature
#       
#       for(i in 1:length(signif_molecules)){
#         
#         molecule_name <- signif_molecules[i]
#         
#         `1_before_correction` <- unlist(data_before_correction[data_before_correction[,1] == molecule_name,])
#         `1_before_correction` <- data.frame(`1_before_correction`) %>% rownames_to_column("samples")
#         `2_after_correction` = unlist(data_after_correction[data_after_correction[,1] == molecule_name,])
#         `2_after_correction` <- data.frame(`2_after_correction`) %>% rownames_to_column("samples")
#         
#         data <- 
#           `1_before_correction` %>% 
#           left_join(
#             `2_after_correction`
#           ) %>% 
#           left_join(
#             annotation_set %>% 
#               mutate(samples = as.character(samples)) %>% 
#               select(samples, sample_age_w, qualified_for_analysis, TSC_control = Condition, comparison = starts_with(comparison_id))
#           ) %>% 
#           filter(qualified_for_analysis %in% c("yes", "dev_epi")) %>% 
#           filter(!is.na(comparison)) %>%
#           distinct(samples, .keep_all = T) %>% 
#           gather(stage, value, c(`X1_before_correction`, `X2_after_correction`)) %>% 
#           mutate(value = as.numeric(value),
#                  sample_age_w = as.numeric(sample_age_w),
#                  stage = str_remove(stage, "^X")) %>% 
#           filter(!is.na(value))
#         
#         scatter_plot <-
#           data %>% 
#           ggplot(aes(x = sample_age_w, y = value)) + 
#           geom_point(aes(color = comparison, shape = TSC_control)) +
#           theme_bw() +
#           theme( legend.position = "bottom", legend.direction = "vertical") +
#           geom_smooth(aes(color = comparison), method = lm) +
#           ggpubr::stat_cor(aes(color = comparison), method = "spearman", label.x.npc = "middle", label.y.npc = "top") +
#           facet_wrap(~stage) +
#           xlab("age (weeks)") 
#         
#         box_plot <-
#           data %>% 
#           ggplot(aes(x = comparison, y = value)) + 
#           geom_boxplot(aes(color = comparison)) +
#           geom_jitter(aes(color = comparison)) +
#           theme_bw() +
#           theme(legend.position = "none") +
#           scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
#           facet_wrap(~stage)
#         
#         plot <-
#           ggpubr::ggarrange(scatter_plot, box_plot, nrow = 2) %>% 
#           ggpubr::annotate_figure(top = ggpubr::text_grob(paste(comparison_id, molecule_name)))
#         
#         ggsave(paste0("~/interim_report/additional/Proteomics/modification_03_03_age_groups_filter/results/", paste(comparison_id, molecule_name), ".png"), width = 9, height = 10, plot = plot)
#         
#       }
#     })
# }
# # Prepare data ------------------------------------------------------------
# 
# data_before_correction <- 
#   read_excel("additional/Proteomics/modification_21_01/Proteomics Results for Testing_21_01.xlsx", 
#              sheet = "before Correction") %>% 
#   separate(ProteinGroup, "Proteins", sep = ";")
# 
# data_after_correction <-
#   read_excel("additional/Proteomics/modification_21_01/Proteomics Results for Testing_21_01.xlsx",
#              sheet = "after Correction") %>% 
#   separate(ProteinGroup, "Proteins", sep = ";")
# 
# 
# # Run function ------------------------------------------------------------
# 
# present_results(results_to_save, data_before_correction, data_after_correction, sample_data)
