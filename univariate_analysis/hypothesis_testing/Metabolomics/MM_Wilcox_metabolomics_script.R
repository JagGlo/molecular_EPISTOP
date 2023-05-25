library(tidyverse)
library(readxl)
library(openxlsx)

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
  dplyr::select(sample_anon, sample_age_w, starts_with("(FT"), status, VGB_treatment,
                -c(contains("FT06"),contains("FT07"),contains("FT08"),contains("FT09"),contains("FT10"))) %>% 
  mutate(`(FT05) VGB V24 (seizures vs no seizures) - modified` = case_when(
    `(FT05) VGB V24 (seizures vs no seizures (any) - central corr)` ==  "treated - Seizures developed" ~  "Seizures developed", 
    `(FT05) VGB V24 (seizures vs no seizures (any) - central corr)` == "not treated - no Seizures developed" ~ "no Seizures developed",
    `(FT05) VGB V24 (seizures vs no seizures (any) - central corr)` == "treated - no Seizures developed" ~ "no Seizures developed"
  )) %>% 
  select(-`(FT05) VGB V24 (seizures vs no seizures (any) - central corr)`)

comparisons <-
  comparisons %>%
  gather(comparison, condition, -sample_anon, -sample_age_w, -status, -VGB_treatment) %>%
  filter(!is.na(condition), !is.na(sample_anon)) %>%
  group_by(comparison) %>%
  filter(!duplicated(sample_anon)) %>%
  mutate(status = if_else(is.na(status), "never seiz - never VGB", status),
         VGB_treatment = if_else(is.na(VGB_treatment), "No", VGB_treatment),
         comparison_id = str_extract_all(comparison, "^\\(FT.{1,2}\\)")) %>% 
  ungroup()

# Prepare data ------------------------------------------------------------

data <-
  raw_data$met_batch_VGB_MM_age_corr

# Statistics comparison function -----------------------------------------

compare_means <- function(df, hypothesis = "two.sided", allow_parametric = TRUE) {
  split_factor <-
    df[[1]]
  assertthat::assert_that(is.factor(split_factor))
  assertthat::assert_that(length(levels(split_factor)) == 2)
  assertthat::noNA(split_factor)
  numerics <-
    df %>%
    select(-1)
  colnames <- names(numerics)
  normal_cols <- # first assumption of t-tests
    numerics %>%
    sapply(function(i) {
      testA <- tryCatch({
        shapiro.test(i[split_factor == levels(split_factor)[1]])
      }, warning = function(w) {
        warning(w)
        return(list(p.value = 1))
      }, error = function(e) {
        warning(e)
        return(list(p.value = 1))
      })
      
      testB <- tryCatch({
        shapiro.test(i[split_factor == levels(split_factor)[2]])
      }, warning = function(w) {
        warning(w)
        return(list(p.value = 1))
      }, error = function(e) {
        warning(e)
        return(list(p.value = 1))
      })
      
      testA$p.value >= 0.05 & testB$p.value >= 0.05
    })
  equal_vars <- # second assuption of t-tests
    numerics %>%
    sapply(function(i) {
      test <- var.test(
        i[split_factor == levels(split_factor)[1]],
        i[split_factor == levels(split_factor)[2]]
      )
      test$p.value >= 0.05
    })
  t_eligible <- normal_cols & equal_vars
  KW_tests <-
    numerics %>%
    lapply(function(i) wilcox.test(
      x = i[split_factor == levels(split_factor)[1]],
      y = i[split_factor == levels(split_factor)[2]],
      alternative = hypothesis
    ))
  t_tests <-
    numerics %>%
    lapply(function(i) t.test(
      x = i[split_factor == levels(split_factor)[1]],
      y = i[split_factor == levels(split_factor)[2]],
      alternative = hypothesis
    ))
  lapply(seq_along(t_eligible), function(x) {
    if (t_eligible[x] == TRUE & allow_parametric & as.numeric(table(split_factor)[1]) >= 30 & as.numeric(table(split_factor)[2]) >= 30) return(t_tests[[x]] %>% broom::tidy()) else return(KW_tests[[x]] %>% broom::tidy())
  }) %>%
    bind_rows() %>%
    select(statistic, p.value, method, alternative) %>%
    mutate(alternative = case_when(
      alternative == "less" ~ paste(levels(split_factor)[1], "<", levels(split_factor)[2]),
      alternative == "greater" ~ paste(levels(split_factor)[1], ">", levels(split_factor)[2]),
      alternative == "two.sided" ~ paste(levels(split_factor)[1], "!=", levels(split_factor)[2])
    )) %>%
    add_column(feature = colnames, .before = 1) %>%
    left_join(
      df %>%
        gather("feature", "val", -!!sym(names(df)[1])) %>%
        group_by(feature, !!sym(names(df)[1])) %>%
        summarise(valid = sum(!is.na(val))) %>%
        group_by(feature) %>%
        summarise(valid = paste(valid, collapse = " + "))
    ) %>%
    mutate(p.value.fdr = p.adjust(p.value, method = "BH")) %>%
    select(feature, valid, statistic, p.value, p.value.fdr, method, alternative)
}

# Volcano plot function ---------------------------------------------------

plot_volcano <- function(results, data, sampling_data, comparison_name) {
  m1 <- as.data.frame(results)
  m1$neglog <- -log10(m1$p.value.fdr)
  m1.sig <- subset(m1, m1$p.value.fdr < 0.05 & (m1$fold >= 1.5 | m1$fold <= 0.66))
  m1.sig.up <- subset(m1.sig, m1.sig$fold >= 1.5)
  m1.sig.down <- subset(m1.sig, m1.sig$fold <= 0.66)
  
  if (nrow(m1.sig) > 0) {
    p <- ggplot(
      data = m1,
      aes(x = log2FoldChange, y = neglog)
    ) +
      geom_point(alpha = 0.4, size = 1.75) +
      xlab("log2 fold change") + ylab("-log10 p-value") +
      theme_bw() +
      theme(legend.position = "bottom", legend.direction = "vertical")
    
    if (nrow(m1.sig.up) > 0) {
      p <-
        p + geom_point(aes(color = "1"), alpha = 0.4, size = 1.75, data = m1.sig.up)
      
      # top3_up <-
      #   m1.sig.up %>%
      #   arrange(desc(fold), p.value.fdr) %>%
      #   head(3) 
      # 
      # data_up <-
      #   data %>% 
      #   select(sample_anon, top3_up$feature) %>%
      #   gather(metabolite, value, -sample_anon) %>%
      #   left_join(sampling_data %>% transmute(age_weeks = sample_age_w, sample_anon, condition)) 
      # 
      # p_scatter_up <-
      #   data_up %>% 
      #   ggplot(aes(x = age_weeks, y = value, color = condition)) +
      #   geom_jitter(height  = 0) +
      #   facet_wrap(~metabolite, scales = "free_y") +
      #   theme_bw() +
      #   theme(legend.position = "bottom")
      # 
      # p_box_up <-
      #   data_up %>%
      #   ggplot(aes(x = condition, y = value, color = condition)) +
      #   geom_boxplot() +
      #   facet_wrap(~metabolite, scales = "free_y") +
      #   theme_bw() +
      #   theme(legend.position = "none")
      # 
      # p_up <- ggpubr::ggarrange(p_scatter_up, p_box_up, nrow = 2) 
      # 
      # ggsave(paste0("~/interim_report/additional/Metabolomics/", comparison_name, "_up_reg.png"), width = 9, height = 10, plot = p_up)
    } else {p <- p + geom_point(aes(color = "1"), alpha = 0)}
    
    if (nrow(m1.sig.down) > 0) {
      p <- p + geom_point(aes(color = "2"), alpha = 0.4, size = 1.75, data = m1.sig.down) 
      
      # top3_down <-
      #   m1.sig.down %>%
      #   arrange(fold, p.value.fdr) %>%
      #   head(3) 
      # 
      # data_down <-
      # data %>% 
      #   select(sample_anon, top3_down$feature) %>%
      #   gather(metabolite, value, -sample_anon) %>%
      #   left_join(sampling_data %>% transmute(age_weeks = sample_age_w, sample_anon, condition)) 
      # 
      # p_scatter_down <-
      #   data_down %>% 
      #   ggplot(aes(x = age_weeks, y = value, color = condition)) +
      #   geom_jitter(height  = 0) +
      #   facet_wrap(~metabolite, scales = "free_y") +
      #   theme_bw() +
      #   theme(legend.position = "bottom")
      # 
      # p_box_down <-
      #   data_down %>% 
      #   ggplot(aes(x = condition, y = value, color = condition)) +
      #   geom_boxplot() +
      #   facet_wrap(~metabolite, scales = "free_y") +
      #   theme_bw() +
      #   theme(legend.position = "none")
      # 
      # p_down <- ggpubr::ggarrange(p_scatter_down, p_box_down, nrow = 2) 
      # 
      # ggsave(paste0("~/interim_report/additional/Metabolomics/", comparison_name, "_down_reg.png"), width = 9, height = 10, plot = p_down)
      
    } else {p <- p + geom_point(aes(color = "2"), alpha = 0)}
    
    p <- p + theme(legend.title = element_blank()) + 
      scale_color_manual(
        labels = c(
          "1" = paste0(nrow(m1.sig.up), " over-expressed metabolites in ", levels(sampling_data$condition)[1], " (n = ", sum(sampling_data$condition == levels(sampling_data$condition)[1]), ")"),
          "2" = paste0(nrow(m1.sig.down), " over-expressed metabolites in ", levels(sampling_data$condition)[2], " (n = ", sum(sampling_data$condition == levels(sampling_data$condition)[2]), ")")
        ),
        values = c("1" = "red3","2" = "blue3"),
        drop = FALSE
      )
    
    ggsave(paste0("~/interim_report/additional/Metabolomics/modification_03_03_age_groups_filter/", comparison_name, ".png"), width = 9, height = 6, plot = p)
  }
}

# Script ------------------------------------------------------------------

#data <- data[,!(colnames(data) %in% excluded_metabolites[[1]])]

final_results <-
  comparisons %>%
  split(., .$comparison) %>%
  lapply(function(x) {
    
    x$condition <- factor(x$condition)
      
      annotation <-
        x %>%
        # filter(VGB_treatment == "No") %>% 
        select(sample_anon, condition)
      
      df <-
        annotation %>% 
        left_join(data) %>% 
        mutate(condition = as.factor(condition))
    tryCatch({
      results <-
        compare_means(df[-1]) 
      
      group_1 <- as.matrix(df[df$condition == levels(df$condition)[1],-c(1:2)])
      group_2 <- as.matrix(df[df$condition == levels(df$condition)[2],-c(1:2)])
      results$group1_med <- matrixStats::colMedians(group_1, na.rm = T)
      results$group2_med <- matrixStats::colMedians(group_2, na.rm = T)
      results$log2FoldChange <-  results$group1_med - results$group2_med
      results$fold <- 2 ^ results$log2FoldChange
      
      name_comp <- str_replace(x$comparison[1], "/", "") 
      
      openxlsx::write.xlsx(results, paste0(here::here("publication_ready/Metabolomics/Wilcox_results/MM_"), name_comp, ".xlsx"))
      
      # plot_volcano(results, df, x, name_comp)
      
      fold_sig <- subset(results, results$p.value.fdr < 0.05 & (results$fold >= 1.5 | results$fold <= 0.66))
      
      fold_sig$comparison <- name_comp
      
      return(fold_sig)
  }, error = function(error_condition) {
    return(NULL)
  })
      
  })

results_to_save <-
  final_results %>% 
  bind_rows()

results_to_save %>% write.xlsx("publication_ready/Metabolomics/Wilcox_results/MM_results.xlsx")

results_to_save %>% 
  count(comparison) %>% 
  write.xlsx("publication_ready/Metabolomics/Wilcox_results/MM_summary.xlsx")

# Checking the results ----------------------------------------------------

present_results <- function(results_csv, before_correction, data_after_correction, annotation_set){
  
  results_csv %>% 
    split(., .$comparison) %>% 
    lapply(function(x) {
      
      comparison_id <-  str_extract(x$comparison[1], "^\\(FT\\d{1,2}\\)")
      signif_molecules <- x$feature
      
      for(i in 1:length(signif_molecules)){
        
        molecule_name <- signif_molecules[i]
        
        `1_before_correction` <- unlist(data_before_correction[data_before_correction[,1] == molecule_name,])
        `1_before_correction` <- data.frame(`1_before_correction`) %>% rownames_to_column("sample_anon")
        `2_after_correction` = unlist(data_after_correction[data_after_correction[,1] == molecule_name,])
        `2_after_correction` <- data.frame(`2_after_correction`) %>% rownames_to_column("sample_anon")
        
        data <- 
          `1_before_correction` %>% 
          left_join(
            `2_after_correction`
          ) %>% 
          left_join(
            annotation_set %>% 
              mutate(sample_anon = as.character(sample_anon)) %>% 
              select(sample_anon, sample_age_w, qualified_for_analysis, TSC_control = Condition, comparison = starts_with(comparison_id))
          ) %>% 
          filter(qualified_for_analysis %in% c("yes", "dev_epi")) %>% 
          filter(!is.na(comparison)) %>%
          distinct(sample_anon, .keep_all = T) %>% 
          gather(stage, value, c(`X1_before_correction`, `X2_after_correction`)) %>% 
          mutate(value = as.numeric(value),
                 sample_age_w = as.numeric(sample_age_w),
                 stage = str_remove(stage, "^X")) %>% 
          filter(!is.na(value))
        
        scatter_plot <-
          data %>% 
          ggplot(aes(x = sample_age_w, y = value)) + 
          geom_point(aes(color = comparison, shape = TSC_control)) +
          theme_bw() +
          theme( legend.position = "bottom", legend.direction = "vertical") +
          geom_smooth(aes(color = comparison), method = lm) +
          ggpubr::stat_cor(aes(color = comparison), method = "spearman", label.x.npc = "middle", label.y.npc = "top") +
          facet_wrap(~stage) +
          xlab("age (weeks)") 
        
        box_plot <-
          data %>% 
          ggplot(aes(x = comparison, y = value)) + 
          geom_boxplot(aes(color = comparison)) +
          geom_jitter(aes(color = comparison)) +
          theme_bw() +
          theme(legend.position = "none") +
          scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
          facet_wrap(~stage)
        
        plot <-
          ggpubr::ggarrange(scatter_plot, box_plot, nrow = 2) %>% 
          ggpubr::annotate_figure(top = ggpubr::text_grob(paste(comparison_id, molecule_name)))
        
        ggsave(paste0("~/interim_report/additional/Metabolomics/modification_03_03_age_groups_filter/results/", paste(comparison_id, molecule_name), ".png"), width = 9, height = 10, plot = plot)
        
      }
    })
}

data_before_correction <- 
  read_excel("additional/Metabolomics/Metabolomics Data before and after Correction.xlsx", 
             sheet = "Batch corrected")

data_after_correction <-
  read_excel("additional/Metabolomics/Metabolomics Data before and after Correction.xlsx", 
             sheet = "Batch VGB Age corrected")



# Run function ------------------------------------------------------------

present_results(results_to_save, data_before_correction, data_after_correction, sample_data)
