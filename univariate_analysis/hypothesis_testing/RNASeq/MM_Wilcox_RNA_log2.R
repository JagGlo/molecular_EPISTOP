library(tidyverse)
library(readxl)
  
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

# Valid values ------------------------------------------------------------
  
RNA_raw <- raw_data$rna_log2_CPM_TMM_MM_age_corr
  
  
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
        #   
        #   top3_up <-
        #     m1.sig.up %>%
        #     filter(fold != 0, fold != as.numeric(Inf)) %>% 
        #     arrange(desc(log2FoldChange), p.value.fdr) %>%
        #     head(3) 
        #   
        #   data_up <-
        #     data %>% 
        #     select(sample_anon, top3_up$feature) %>%
        #     gather(gene, value, -sample_anon) %>%
        #     left_join(sampling_data %>% transmute(age_weeks = sample_age_w, sample_anon, condition)) %>% 
        #     left_join(convertGeneIDs %>% select(ensembl_gene_id, hgnc_symbol), by = c("gene" = "ensembl_gene_id")) %>% 
        #     mutate(gene = if_else(!is.na(hgnc_symbol), hgnc_symbol, gene))
        #   
        #   p_scatter_up <-
        #     data_up %>% 
        #     ggplot(aes(x = age_weeks, y = value, color = condition)) +
        #     geom_jitter() +
        #     facet_wrap(~gene, scales = "free_y") +
        #     theme_bw() +
        #     theme(legend.position = "bottom")
        #   
        #   p_box_up <-
        #     data_up %>%
        #     ggplot(aes(x = condition, y = log2(value), color = condition)) +
        #     geom_boxplot() +
        #     facet_wrap(~gene, scales = "free_y") +
        #     theme_bw() +
        #     theme(legend.position = "none")
        #   
        #   p_up <- ggpubr::ggarrange(p_scatter_up, p_box_up, nrow = 2) 
        #   
        #   ggsave(paste0("~/interim_report/additional/RNASeq/modification_03_03_age_groups_filter/", comparison_name, "_up_reg.png"), width = 9, height = 10, plot = p_up)
      } else {p <- p + geom_point(aes(color = "1"), alpha = 0)}
      # 
      if (nrow(m1.sig.down) > 0) {
        p <- p + geom_point(aes(color = "2"), alpha = 0.4, size = 1.75, data = m1.sig.down)
        
        #   top3_down <-
        #     m1.sig.down %>%
        #     filter(fold != 0, fold != as.numeric(Inf)) %>% 
        #     arrange(log2FoldChange, p.value.fdr) %>%
        #     head(3) 
        #   
        #   data_down <-
        #     data %>% 
        #     select(sample_anon, top3_down$feature) %>%
        #     gather(gene, value, -sample_anon) %>%
        #     left_join(sampling_data %>% transmute(age_weeks = sample_age_w, sample_anon, condition)) %>% 
        #     left_join(convertGeneIDs %>% select(ensembl_gene_id, hgnc_symbol), by = c("gene" = "ensembl_gene_id")) %>% 
        #     mutate(gene = if_else(!is.na(hgnc_symbol), hgnc_symbol, gene))
        #   
        #   p_scatter_down <-
        #     data_down %>% 
        #     ggplot(aes(x = age_weeks, y = value, color = condition)) +
        #     geom_jitter() +
        #     facet_wrap(~gene, scales = "free_y") +
        #     theme_bw() +
        #     theme(legend.position = "bottom")
        #   
        #   p_box_down <-
        #     data_down %>% 
        #     ggplot(aes(x = condition, y = log2(value), color = condition)) +
        #     geom_boxplot() +
        #     facet_wrap(~gene, scales = "free_y") +
        #     theme_bw() +
        #     theme(legend.position = "none")
        #   
        #   p_down <- ggpubr::ggarrange(p_scatter_down, p_box_down, nrow = 2) 
        #   
        #   ggsave(paste0("~/interim_report/additional/RNASeq/modification_03_03_age_groups_filter/", comparison_name, "_down_reg.png"), width = 9, height = 10, plot = p_down)
        #   
      } else {p <- p + geom_point(aes(color = "2"), alpha = 0)}
      
      p <- p + theme(legend.title = element_blank()) + 
        scale_color_manual(
          labels = c(
            "1" = paste0(nrow(m1.sig.up), " over-expressed genes in ", levels(sampling_data$condition)[1], " (n = ", sum(sampling_data$condition == levels(sampling_data$condition)[1]), ")"),
            "2" = paste0(nrow(m1.sig.down), " over-expressed genes in ", levels(sampling_data$condition)[2], " (n = ", sum(sampling_data$condition == levels(sampling_data$condition)[2]), ")")
          ),
          values = c("1" = "red3","2" = "blue3"),
          drop = FALSE
        )
      
      ggsave(paste0("~/interim_report/additional/RNASeq/modification_03_03_age_groups_filter/", comparison_name, ".png"), width = 9, height = 6, plot = p)
    }
  }
  
  
  
  # Compare means -----------------------------------------------------------
  
  compare_means <- function(df, hypothesis = "two.sided") {
    split_factor <-
      df[[1]]
    assertthat::assert_that(is.factor(split_factor))
    assertthat::assert_that(length(levels(split_factor)) == 2)
    assertthat::noNA(split_factor)
    numerics <-
      df %>%
      dplyr::select(-1)
    colnames <- names(numerics)
    KW_tests <-
      numerics %>%
      lapply(function(i) wilcox.test(
        x = i[split_factor == levels(split_factor)[1]],
        y = i[split_factor == levels(split_factor)[2]],
        alternative = hypothesis
      ))
    
    lapply(KW_tests, function(x) {
      x %>% broom::tidy()
    }) %>%
      bind_rows() %>%
      dplyr::select(statistic, p.value, method, alternative) %>%
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
      dplyr::select(feature, valid, statistic, p.value, p.value.fdr, method, alternative)
  }
  
  # Pathway analysis / Gene ontology ----------------------------------------
  
  folds_pa_go <- function(fold_sig, name_comp) {
    
    fold_con <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "gene_biotype"), filters = "ensembl_gene_id", values = fold_sig$feature, mart = ensembl)
    
    fold_pathways <- enrichPathway(gene = fold_con$entrezgene_id, pvalueCutoff = 0.05, readable = T, organism = "human", pAdjustMethod = "BH")
    
    fold_pathways@result %>% write_csv(paste0("~/interim_report/additional/RNASeq/final_filtered/", name_comp, "_pa.csv"))
    
    genes_found_PA <- sum(fold_pathways@result$p.adjust <= fold_pathways@pvalueCutoff)
    
    if (genes_found_PA > 0) {
      fold_entrezgene <- fold_sig$log2FoldChange
      names(fold_entrezgene) <- fold_con$entrezgene_id
      
      dotplot(fold_pathways, showCategory = 20)
      ggsave(paste0("~/interim_report/additional/RNASeq/final_filtered/", name_comp, "_pa_dot.png"), width = 10, height = 15)
      cnetplot(fold_pathways, showCategory = 10, foldChange = fold_entrezgene)
      ggsave(paste0("~/interim_report/additional/RNASeq/final_filtered/", name_comp, "_pa_cnet.png"), width = 16, height = 12)
    }
    
    fold_GO <- enrichGO(gene = fold_con$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    
    fold_GO@result %>% write_csv(paste0("~/interim_report/additional/RNASeq/final_filtered/", name_comp, "_go.csv"))
    
    genes_found_GO <- sum(fold_GO@result$p.adjust <= fold_GO@pvalueCutoff)
    
    if (genes_found_GO > 0) {
      fold_entrezgene <- fold_sig$log2FoldChange
      names(fold_entrezgene) <- fold_con$entrezgene_id
      dotplot(fold_GO, showCategory = 20)
      ggsave(paste0("~/interim_report/additional/RNASeq/final_filtered/", name_comp, "_go_dot.png"), width = 10, height = 15)
      cnetplot(fold_GO, showCategory = 10, foldChange = fold_entrezgene)
      ggsave(paste0("~/interim_report/additional/RNASeq/final_filtered/", name_comp, "_go_cnet.png"), width = 16, height = 12)
    }
  }
  
  
  # Script ------------------------------------------------------------------
  data <- RNA_raw
  
  final_results <-
    comparisons %>%
    split(., .$comparison) %>%
    lapply(function(x) {
      
      x$condition <- factor(x$condition)
      
      RNA <- RNA_raw[, colnames(RNA_raw) %in% x$sample_anon]
      
      annotation <-
        x %>%
        filter(sample_anon %in% colnames(RNA)) %>% # not always sample_anon are found in data
        dplyr::select(-comparison, -sample_age_w)
      # more than 50% of valid values inside groups
      group_1 <- RNA[,colnames(RNA) %in% annotation$sample_anon[annotation$condition == levels(annotation$condition)[1]]]
      group_1 <- group_1[rowSums(group_1 > 1) > ncol(group_1)/2,]
      group_2 <- RNA[,colnames(RNA) %in% annotation$sample_anon[annotation$condition == levels(annotation$condition)[2]]]
      group_2 <- group_2[rowSums(group_2 > 1) > ncol(group_2)/2,]
      
      RNA <-
        RNA[na.omit(unique(c(rownames(group_1), rownames(group_2)))),] %>% 
        t(.) %>% 
        as.data.frame() %>% 
        rownames_to_column() %>%
        dplyr::rename(sample_anon = rowname) %>% 
        gather(gene_id, value, -sample_anon) 
      
      RNA <-
        RNA %>% 
        # anti_join(
        #   sample_anon_to_filter %>%
        #     select(gene_id, sample_anon)) %>%
        spread(gene_id, value)
      
      df <- 
        annotation %>% 
        select(sample_anon, condition) %>% 
        left_join(RNA, by = c("sample_anon")) 
      
      results <-
        compare_means(df[-1]) 
      
      group_1 <- as.matrix(df[df$condition == levels(df$condition)[1],-c(1:2)])
      group_2 <- as.matrix(df[df$condition == levels(df$condition)[2],-c(1:2)])
      results$group1_med <- matrixStats::colMedians(group_1, na.rm = T)
      results$group2_med <- matrixStats::colMedians(group_2, na.rm = T)
      results$log2FoldChange <-  results$group1_med - results$group2_med
      results$fold <- 2 ^ results$log2FoldChange
      
      name_comp <- str_replace(x$comparison[1], "/", "") 
      
      xlsx::write.xlsx(results, paste0(here::here("publication_ready/RNASeq/Wilcox_results/MM_"), name_comp, ".xlsx"))
      
      # plot_volcano(results, df, x, name_comp)
      
      fold_sig <- subset(results, results$p.value.fdr < 0.05 & (results$fold >= 1.5 | results$fold <= 0.66))
      
      #if (nrow(fold_sig) > 50) {
      #   #folds_pa_go(fold_sig, name_comp)
      #}
      
      fold_sig$comparison <- name_comp
      
      return(fold_sig)
      
    })
  
  results_to_save <-
    final_results %>% 
    bind_rows()
  
  results_to_save <- 
    results_to_save %>%  
    left_join(convertGeneIDs %>% dplyr::select(ensembl_gene_id, hgnc_symbol), by = c("feature" = "ensembl_gene_id")) %>% 
    mutate(gene = if_else(!is.na(hgnc_symbol), hgnc_symbol, feature)) 
  
  results_to_save %>% xlsx::write.xlsx("publication_ready/RNASeq/Wilcox_results/MM_results.xlsx")
  
  results_to_save %>% 
    count(comparison) %>% 
    xlsx::write.xlsx("publication_ready/RNASeq/Wilcox_results/MM_summary.xlsx")
  
  #RNA_raw %>% as.data.frame() %>% write_csv("~/interim_report/additional/RNASeq/modification_03_03_age_groups_filter/log2_after.csv")
  
  results_to_save %>% count(comparison) %>% View()
  
  results_to_save %>% mutate(gene = if_else(!is.na(hgnc_symbol), hgnc_symbol, feature)) %>% count(gene) %>% arrange(desc(n)) %>% View()
  
  
  # Checking the results ----------------------------------------------------
  
  present_results <- function(results_csv, before_correction, data_after_correction, annotation_set){
    
    results_csv %>% 
      split(., .$comparison) %>% 
      lapply(function(x) {
        
        comparison_id <-  str_extract(x$comparison[1], "^\\(FT\\d{1,2}\\)")
        signif_molecules <- x$feature
        names <- x$gene
        
        for(i in 1:length(signif_molecules)){
          #rename(!!sym(column_id) := common_dictionary_row_id, !!sym(new_name) := value)
          
          molecule_name <- signif_molecules[i]
          pretty_name <- names[i]
          
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
            ggpubr::annotate_figure(top = ggpubr::text_grob(paste(comparison_id, pretty_name)))
          
          ggsave(paste0("~/interim_report/additional/RNASeq/modification_03_03_age_groups_filter/results/", paste(comparison_id, pretty_name), ".png"), width = 9, height = 10, plot = plot)
          
        }
      })
  }
  
  data_before_correction <- #RNA_raw %>% as.data.frame() %>% rownames_to_column(var = "gene") 
    #   read_excel("additional/RNASeq/log2_RNAseq_CPM_TMM (before and after age correction).xlsx", 
    #              sheet = "before Correction") 
    
    data_after_correction <- RNA_raw %>% as.data.frame() %>% rownames_to_column(var = "gene") 
  
  
  present_results(results_to_save, data_before_correction, data_after_correction, sample_data)
  