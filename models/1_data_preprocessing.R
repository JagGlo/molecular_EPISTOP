# Libraries ---------------------------------------------------------------
library(mlr)

# Data --------------------------------------------------------------------

source("~/interim_report/reports/data_preprocessing.R") # preprocessed data & additional files
#patients$seizures_any <- sample(patients$seizures_any, replace = F)
#patients$seizures_any <- sample(patients$seizures_any, replace = F)

# Eligible blood samples --------------------------------------------------

eligible_samples <-
  sampling_data %>% 
  filter(correct_name %in% c("V1", "V1/VaEEG")) %>% 
  filter(!(code %in% c("01-025", "01-017", "01-059", "01-032", "04-001", "09-005", "09-006"))) %>% # excluded
  transmute(code, 
            sample_code = correct_name,
            sample_date) %>% 
  group_by(code, sample_code) %>% 
  slice(1) %>% # one patient has two V1 samples looking good
  ungroup() %>% 
  left_join(patients %>% select(code, first_SS, first_AEDrug, seizures_any)) %>% 
  filter(first_SS > sample_date | is.na(first_SS)) %>% # remove "subclinical seizures" samples
  mutate(status = case_when(
    is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - never VGB",
    seizures_any == TRUE ~ "ever seizure",
    !is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - VGB treated"
    ),
    seizures_any = seizures_any %>% as.factor()) %>% 
  filter(status != "never seiz - VGB treated") %>% 
  select(-status, -first_SS, -first_AEDrug, -seizures_any)

raw_values <- eligible_samples

eligible_samples <- eligible_samples %>% select(-sample_date)

# eligible_samples %>% 
#   count(status)

# status                         n
# 1 ever seizure                54
# 2 never seiz - never VGB      14
# 3 never seiz - VGB treated    14

# Exclude SS samples ------------------------------------------------------

# eligible_samples %>% 
#   filter(first_SS <= sample_date)


# MRI ---------------------------------------------------------------------

MRI <-
  eligible_samples %>% 
  inner_join(
    mri_baseline %>%
      mutate(
        code = str_extract(`Pt ID`, "\\d{2}-\\d{3}"),
        TBV = TBV %>% as.numeric()
      ) %>%
      filter(!(code %in% c("01-025", "01-017", "01-059","01-032", "04-001", "09-005", "09-006"))) %>% # excluded patients
      mutate_at(vars(`Tuber y/n`:Calcification), funs(as.numeric(case_when(. == "no" ~ 0, . == "yes" ~ 1)))) %>%
      select(-c(`Pt ID`:`Timepoint MRI`)),
    by = "code"
  ) %>% 
  mutate(sample_code = if_else(is.na(sample_code), "V1", sample_code)) %>%
  mlr::removeConstantFeatures() %>%
  mlr::normalizeFeatures(method = "standardize")

colnames(MRI) <- make.names(names(MRI)) # convert names to R standard

assert_that(n_distinct(MRI %>% select(code, sample_code)) == nrow(MRI))
assert_that(any(duplicated(MRI %>% select(code, sample_code))) == FALSE)

# multiple missings in this dataset;  "yes/no" questions -> 93/93  completed, some variables -> over 90% completed (87/93)

MRI_over_90pcent <-
  MRI %>% 
  select(-c(TBV:Lesion_Brain))

MRI_only_completed <-
  MRI %>% 
  select(-c(Tuber_vol:Lesion_Brain))

# proteom -----------------------------------------------------------------

proteome <-
  eligible_samples %>%
  inner_join(
    protein %>%
      left_join(protein_meta %>% select(protein_ID, protein_id)) %>%
      select(-protein_ID) %>%
      spread(protein_id, intensity) %>%
      left_join(protein_sample_meta %>% transmute(sample, sample_code = new_sample_type, code)) %>%
      mutate(code = case_when(
        is.na(code) ~ sample,
        TRUE ~ code
      )) %>%
      select(-sample)
  ) %>%
  mlr::removeConstantFeatures() %>% 
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(proteome %>% select(code, sample_code)) == nrow(proteome))
assert_that(any(duplicated(proteome %>% select(code, sample_code))) == FALSE)

# proteom_V1_V1/VaEEG -----------------------------------------------------------------

proteome_first_tp <-
  eligible_samples %>%
  inner_join(
    proteome_V1_V1VaEEG
  ) %>%
  mlr::removeConstantFeatures()

raw_values <- raw_values %>% left_join(proteome_first_tp)

proteome_first_tp <-
  proteome_first_tp %>% 
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(proteome_first_tp %>% select(code, sample_code)) == nrow(proteome_first_tp))
assert_that(any(duplicated(proteome_first_tp %>% select(code, sample_code))) == FALSE)

# proteom_raw -----------------------------------------------------------------

proteome_raw <-
  eligible_samples %>%
  inner_join(
    proteome_raw
  ) %>%
  mlr::removeConstantFeatures() %>% 
  left_join(patients %>% select(code, seizures_any), by = "code")

proteom_seiz_T <- as.matrix(proteome_raw %>% filter(seizures_any == T) %>% select(-code, -sample_code, -seizures_any))
proteom_seiz_F <- as.matrix(proteome_raw %>% filter(seizures_any == F) %>% select(-code, -sample_code, -seizures_any))

# more than 3 reads in at least 75% of samples within groups

proteom_T_treshold <- proteom_seiz_T > 0.05
proteom_T_treshold <-
  as.data.frame(colSums(proteom_T_treshold)) %>%
  rownames_to_column() %>%
  filter(as.numeric(`colSums(proteom_T_treshold)`) >= 0.75 * nrow(proteom_seiz_T))

proteom_F_treshold <- proteom_seiz_F > 0.05
proteom_F_treshold <-
  as.data.frame(colSums(proteom_F_treshold)) %>%
  rownames_to_column() %>%
  filter(as.numeric(`colSums(proteom_F_treshold)`) >= 0.75 * nrow(proteom_seiz_F))

selected_genes <-
  bind_rows(
    proteom_F_treshold %>% transmute(seizure = rowname),
    proteom_T_treshold %>% transmute(seizure = rowname),
    .id = "id"
  ) %>%
  mutate(id = ifelse(id == 1, "seiz_F", "seiz_T"))

selected_genes <- unique(selected_genes$seizure)

proteom_filtered <- data.frame(seizures_any = proteome_raw$seizures_any,
                               code = proteome_raw$code, 
                               sample_code = proteome_raw$sample_code, 
                               proteome_raw[, selected_genes])

proteome_seiz_T_filter <- as.matrix(proteom_filtered %>% filter(seizures_any == T) %>% 
                                 select(-code, -sample_code, -seizures_any))

proteome_seiz_F_filter <- as.matrix(proteom_filtered %>% filter(seizures_any == F) %>% 
                                      select(-code, -sample_code, -seizures_any))

# fold changes
fold <- matrixStats::colMedians(proteome_seiz_T_filter) / matrixStats::colMedians(proteome_seiz_F_filter)
fold <- as.data.frame(fold)
rownames(fold) <- colnames(proteome_seiz_T_filter)
fold <-
  fold %>%
  rownames_to_column() %>% 
  filter(fold >= 2 | fold <= 0.5) 

# raw_values <- raw_values %>% left_join(proteome_raw)

proteome_raw <-
  proteom_filtered %>% 
  select(-seizures_any) %>% 
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(proteome_raw %>% select(code, sample_code)) == nrow(proteome_raw))
assert_that(any(duplicated(proteome_raw %>% select(code, sample_code))) == FALSE)

rm(proteom_seiz_T, proteom_seiz_F, proteom_T_treshold, proteom_F_treshold, 
   proteom_filtered, proteome_seiz_T_filter, proteome_seiz_F_filter)

# RNA ---------------------------------------------------------------------

RNA <-
  eligible_samples %>%
  inner_join(
    RNA_raw %>%
      gather(sample_id, expression, -X1) %>%
      spread(X1, expression) %>%
      left_join(RNA_sample_meta %>% transmute(sample_id = GSID, code, sample_code = correct_name)) %>%
      select(-sample_id)
  ) %>%
  group_by(code, sample_code) %>%
  filter(!(duplicated(sample_code))) %>%
  ungroup() %>%
  mlr::removeConstantFeatures() %>% 
  left_join(patients %>% select(code, seizures_any), by = "code")

# filtering
RNA_seiz_T <- as.matrix(RNA %>% filter(seizures_any == T) %>% select(-code, -sample_code, -seizures_any))
RNA_seiz_F <- as.matrix(RNA %>% filter(seizures_any == F) %>% select(-code, -sample_code, -seizures_any))

# more than 3 reads in at least 75% of samples within groups

RNA_T_treshold <- RNA_seiz_T > 0.05
RNA_T_treshold <-
  as.data.frame(colSums(RNA_T_treshold)) %>%
  rownames_to_column() %>%
  filter(as.numeric(`colSums(RNA_T_treshold)`) >= 0.75 * nrow(RNA_seiz_T))

RNA_F_treshold <- RNA_seiz_F > 0.05
RNA_F_treshold <-
  as.data.frame(colSums(RNA_F_treshold)) %>%
  rownames_to_column() %>%
  filter(as.numeric(`colSums(RNA_F_treshold)`) >= 0.75 * nrow(RNA_seiz_F))

selected_genes <-
  bind_rows(
    RNA_F_treshold %>% transmute(seizure = rowname),
    RNA_T_treshold %>% transmute(seizure = rowname),
    .id = "id"
  ) %>%
  mutate(id = ifelse(id == 1, "seiz_F", "seiz_T"))

selected_genes <- unique(selected_genes$seizure)

RNA_filtered <- data.frame(RNA$code, RNA$sample_code, RNA$seizures_any, RNA[, selected_genes])

RNA_seiz_T_filter <- as.matrix(RNA_filtered %>% filter(RNA.seizures_any == T) %>% 
                                 select(-RNA.code, -RNA.sample_code, -RNA.seizures_any))

RNA_seiz_F_filter <- as.matrix(RNA_filtered %>% filter(RNA.seizures_any == F) %>% 
                                 select(-RNA.code, -RNA.sample_code, -RNA.seizures_any))

# fold changes
fold <- matrixStats::colMedians(RNA_seiz_T_filter) / matrixStats::colMedians(RNA_seiz_F_filter)
fold <- as.data.frame(fold)
rownames(fold) <- colnames(RNA_seiz_F_filter)
fold <-
  fold %>%
  rownames_to_column() %>% 
  filter(fold >= 2 | fold <= 0.5) %>%
  filter(fold != 0)

filtered_RNA <- c("code", "sample_code", fold$rowname)

RNA <-
  RNA[, filtered_RNA]

raw_values <- raw_values %>% left_join(RNA)

RNA <-
  RNA %>% 
  mlr::normalizeFeatures(method = "standardize")

colnames(RNA) <- make.names(names(RNA))

assert_that(n_distinct(RNA %>% select(code, sample_code)) == nrow(RNA))
assert_that(any(duplicated(RNA %>% select(code, sample_code))) == FALSE)

rm(RNA_seiz_T, RNA_seiz_F, RNA_T_treshold, RNA_F_treshold, selected_genes, RNA_filtered, 
   RNA_seiz_T_filter, RNA_seiz_F_filter, filtered_RNA)


# isomiRs -----------------------------------------------------------------

isomiR_meta <- 
  read_csv("~/interim_report/input/WP3/WP3_isomiR_annotation_norm.csv") %>% 
  rename(GSID = X1) %>% 
  filter(Notes == "none") %>% 
  mutate(code = str_extract(PatientID, "\\d{2}-\\d{3}")) %>% 
  left_join(sampling_data %>% select(code, sample_code, correct_name), by = c("code", "Time_Point_2" = "sample_code")) %>% 
  mutate(correct_name = if_else(!is.na(code) & is.na(correct_name), Time_Point_2, correct_name),
         correct_name = case_when(
           Condition == "control" & Time_Point_1 == "V1" ~ "Early_Control",
           Condition == "control" & Time_Point_1 == "V24" ~ "Late_Control",
           TRUE ~ correct_name)) %>% 
  filter(!duplicated(GSID)) #removes duplicated sample o1-015 V1 & 16-007 V24

isomiR_data <- 
  read_csv("~/interim_report/input/WP3/WP3_isomiR_normCountCPM.csv") %>% 
  rename(isomiR = X1) %>% 
  gather(sample, value, -isomiR) %>% 
  spread(isomiR, value) %>% 
  left_join(isomiR_meta %>% transmute(GSID, code, sample_code = correct_name), by = c("sample" = "GSID")) 

isomiR <-
  eligible_samples %>%
  inner_join(
    isomiR_data) %>%
  left_join(
    patients %>% select(code, seizures_any)
  ) %>% 
  select(-sample)

colnames(isomiR) <- make.names(names(isomiR)) # convert names to R standard

# filtering
isomiR_seiz_T <- as.matrix(isomiR %>% filter(seizures_any == T) %>% select(-code, -sample_code, -seizures_any)) 
isomiR_seiz_F <- as.matrix(isomiR %>% filter(seizures_any == F) %>% select(-code, -sample_code, -seizures_any))

# more than 3 reads in at least 75% of samples within groups

isomiR_T_treshold <- isomiR_seiz_T > 0.625
isomiR_T_treshold <-
  as.data.frame(colSums(isomiR_T_treshold)) %>%
  rownames_to_column() %>%
  filter(as.numeric(`colSums(isomiR_T_treshold)`) >= 0.75 * nrow(isomiR_seiz_T) & `colSums(isomiR_T_treshold)` > 0)

isomiR_F_treshold <- isomiR_seiz_F > 0.625
isomiR_F_treshold <-
  as.data.frame(colSums(isomiR_F_treshold)) %>%
  rownames_to_column() %>%
  filter(as.numeric(`colSums(isomiR_F_treshold)`) >= 0.75 * nrow(isomiR_seiz_F) & `colSums(isomiR_F_treshold)` > 0)

selected_genes <-
  bind_rows(
    isomiR_F_treshold %>% transmute(seizure = rowname),
    isomiR_T_treshold %>% transmute(seizure = rowname),
    .id = "id"
  ) %>%
  mutate(id = ifelse(id == 1, "seiz_F", "seiz_T"))

selected_genes <- unique(selected_genes$seizure)

isomiR_filtered <- data.frame(isomiR$code, isomiR$sample_code, isomiR$seizures_any, isomiR[, selected_genes])

isomiR_seiz_T_filter <- as.matrix(isomiR_filtered %>% filter(isomiR.seizures_any == T) %>% 
                                 select(-isomiR.code, -isomiR.sample_code, -isomiR.seizures_any))

isomiR_seiz_F_filter <- as.matrix(isomiR_filtered %>% filter(isomiR.seizures_any == F) %>% 
                                    select(-isomiR.code, -isomiR.sample_code, -isomiR.seizures_any))

# fold changes
fold <- if (length(isomiR_seiz_F_filter) > 0 & length(isomiR_seiz_T_filter) > 0) matrixStats::colMedians(isomiR_seiz_T_filter) / 
  matrixStats::colMedians(isomiR_seiz_F_filter) else matrixStats::colMedians(isomiR_seiz_T_filter) / 
  length(isomiR_seiz_F_filter)

fold <- as.data.frame(fold)
rownames(fold) <- colnames(isomiR_seiz_F_filter)
fold <-
  fold %>%
  rownames_to_column() %>% 
  filter(fold >= 2 | fold <= 0.5) %>%
  filter(fold != 0)

filtered_isomiR <- c("code", "sample_code", fold$rowname)

isomiR <-
  isomiR[, filtered_isomiR]

raw_values <- raw_values %>% left_join(isomiR)

isomiR <-
  isomiR %>% 
  mlr::removeConstantFeatures() %>% 
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(isomiR %>% select(code, sample_code)) == nrow(isomiR))
assert_that(any(duplicated(isomiR %>% select(code, sample_code))) == FALSE)

rm(isomiR_seiz_T, isomiR_seiz_F, isomiR_T_treshold, isomiR_F_treshold, selected_genes, isomiR_filtered, 
   isomiR_seiz_T_filter, isomiR_seiz_F_filter, filtered_isomiR, fold)

# miRNA -------------------------------------------------------------------

miRNA <-
  eligible_samples %>%
  inner_join(
    miRNA_raw %>%
      gather(sample_id, expression, -X1) %>%
      spread(X1, expression) %>%
      mutate(sample_id = sample_id %>% as.numeric()) %>%
      left_join(
        miRNA_sample_meta %>%
          transmute(sample_id = `Sample no. AMC`, code, sample_code = correct_name)
      ) %>%
      select(-sample_id, -`<NA>`)
  ) %>%
  mlr::removeConstantFeatures()

miRNA <-
  miRNA[!duplicated(miRNA[,c("code", "sample_code")]),]

raw_values <- raw_values %>% left_join(miRNA)

miRNA <-
  miRNA %>% 
  mlr::normalizeFeatures(method = "standardize")

colnames(miRNA) <- make.names(names(miRNA)) # convert names to R standard

assert_that(n_distinct(miRNA %>% select(code, sample_code)) == nrow(miRNA))
assert_that(any(duplicated(miRNA %>% select(code, sample_code))) == FALSE)

# metabolite --------------------------------------------------------------

metabolite <-
  eligible_samples %>%
  inner_join(
    metabolite_data_raw %>%
      mutate(
        sample_code = case_when(
          str_detect(Timepoint, "OUT") ~ "out",
          TRUE ~ Timepoint
        ),
        vgb = as.numeric(ifelse(vgb == "Yes", 1, 0))
      ) %>%
      select(-c(Batch:`VGB (0=no, 1=yes)`))
  ) %>% 
  replace(is.na(.), 0) %>%
  mlr::removeConstantFeatures()

raw_values <- raw_values %>% left_join(metabolite)

metabolite <-
  metabolite %>% 
  mlr::normalizeFeatures(method = "standardize")

colnames(metabolite) <- make.names(names(metabolite))

assert_that(n_distinct(metabolite %>% select(code, sample_code)) == nrow(metabolite))
assert_that(any(duplicated(metabolite %>% select(code, sample_code))) == FALSE)

# TSC ---------------------------------------------------------------------

TSC <-
  eligible_samples %>%
  right_join(
    mutation_data %>%
      select(code, gene, mutant_allele_freq) %>% 
      transmute(code = code, 
                TSC = case_when(
        gene == "TSC2" & mutant_allele_freq < 0.5 ~ "TSC2_m",
        gene == "TSC2" ~ "TSC2",
        gene == "TSC1" ~ "TSC1",
        gene == "NMI" ~ "NMI"
      ) %>% as.factor)
  ) %>%
  mutate(sample_code = if_else(is.na(sample_code), "V1", sample_code), ) %>%
  filter(!(code %in% c("01-025", "01-017", "01-059", "01-032", "04-001", "09-005", "09-006"))) %>% # excluded patients
  mutate_if(is.logical, as.numeric) %>%
  mutate_if(is.factor, fct_drop) %>%
  inner_join(
      patients %>% select(code, first_SS, first_AEDrug, seizures_any) %>% 
        mutate(status = case_when(
          is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - never VGB",
          seizures_any == TRUE ~ "ever seizure",
          !is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - VGB treated"
        ),
        seizures_any = seizures_any %>% as.factor()) %>% 
        filter(status != "never seiz - VGB treated") %>% 
        select(code)) %>%  
  mlr::createDummyFeatures(method = "reference") %>%
  mlr::removeConstantFeatures()

raw_values <- raw_values %>% left_join(TSC)

TSC <-
  TSC %>% 
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(TSC %>% select(code, sample_code)) == nrow(TSC))
assert_that(any(duplicated(TSC %>% select(code, sample_code))) == FALSE)


# WGS ---------------------------------------------------------------------

WGS_SNP <-
  read_excel("~/interim_report/input/WP3/WP3_WGS_SNP_genotypes.xlsx",
             na = "NA"
  )

names(WGS_SNP) <- paste(names(WGS_SNP), WGS_SNP[1, ], sep = "_")
WGS_SNP <-
  eligible_samples %>%
  right_join(
  WGS_SNP[-1, -1] %>%
    rename(patient_id = `...2_patient id`, seizure_status = `...3_seizure status`) %>% 
  separate(patient_id, c("wgs_code", "patient_id", "additional"), sep = "_") %>% 
  mutate(code = str_extract(`patient_id`, "\\d{2}-\\d{3}"),
         code = case_when(
           wgs_code == "7042-011-001-021" ~ "06-002",
           code == "08-003" ~ "09-003",
           TRUE ~ code
         )) %>% 
    filter(!is.na(code)) %>% 
    filter(!(code %in% c("01-025", "01-017", "01-059","01-032", "04-001", "09-005", "09-006"))),  # excluded patients,
  by = "code"
  ) %>% 
  mutate(sample_code = if_else(is.na(sample_code), "V1", sample_code)) %>% 
  select(-wgs_code, -patient_id, -additional, -seizure_status) %>% 
  inner_join(
    patients %>% select(code, first_SS, first_AEDrug, seizures_any) %>% 
      mutate(status = case_when(
        is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - never VGB",
        seizures_any == TRUE ~ "ever seizure",
        !is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - VGB treated"
      ),
      seizures_any = seizures_any %>% as.factor()) %>% 
      filter(status != "never seiz - VGB treated") %>% 
      select(code)) %>% 
  mutate_at(vars(`rs346291 (C/T)_chr6:80564836`:`rs1130183 (G/A)_chr1:160011512`), as.numeric) %>%  
  mlr::removeConstantFeatures()

raw_values <- raw_values %>% left_join(WGS_SNP)

WGS_SNP <-
  WGS_SNP %>% 
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(WGS_SNP %>% select(code, sample_code)) == nrow(WGS_SNP))
assert_that(any(duplicated(WGS_SNP %>% select(code, sample_code))) == FALSE)

colnames(WGS_SNP) <- make.names(names(WGS_SNP)) # convert names to R standard

# eCRF --------------------------------------------------------------------

eCRF <-
  eligible_samples %>%
  right_join(
    patients %>% transmute(code, sexMale = if_else(sex == "Male", 1, 0), dob)
  ) %>%
  mutate(sample_code = if_else(is.na(sample_code), "V1", sample_code), 
  is_aEEG = if_else(sample_code == "V1/VaEEG", 1, 0)) %>% 
  filter(!(code %in% c("01-025","01-032", "01-017", "01-059", "04-001", "09-005", "09-006"))) %>% # excluded patients
  inner_join(
    patients %>% select(code, first_SS, first_AEDrug, seizures_any) %>% 
      mutate(status = case_when(
        is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - never VGB",
        seizures_any == TRUE ~ "ever seizure",
        !is.na(first_AEDrug) & seizures_any == FALSE ~ "never seiz - VGB treated"
      ),
      seizures_any = seizures_any %>% as.factor()) %>% 
      filter(status != "never seiz - VGB treated") %>% 
      select(code)) %>%
  mutate_if(is.logical, as.numeric) %>%
  mutate_if(is.factor, fct_drop) %>%
  mlr::removeConstantFeatures()

raw_values <- 
  raw_values %>% 
  left_join(eCRF) %>% 
  mutate(age_weeks = round(as.numeric(sample_date - dob)/7,1)) 

eCRF <-
  eCRF %>% 
  select(-dob) %>% 
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(eCRF %>% select(code, sample_code)) == nrow(eCRF))
assert_that(any(duplicated(eCRF %>% select(code, sample_code))) == FALSE)


# HMGB --------------------------------------------------------------------

HMGB <-
  eligible_samples %>% 
  inner_join(
    inflammatory_data %>% 
      select(code, sample_code, hmgb_level)
  ) %>% 
  mlr::removeConstantFeatures()

HMGB <-
  HMGB[!duplicated(HMGB[,c("code", "sample_code")]),] %>%
  mlr::normalizeFeatures(method = "standardize")

assert_that(n_distinct(HMGB %>% select(code, sample_code)) == nrow(HMGB))
assert_that(any(duplicated(HMGB %>% select(code, sample_code))) == FALSE)


# Making list -------------------------------------------------------------

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

# seizures_any_clinical      n
# 1 FALSE                    38
# 2 TRUE                     59
  
all_completed <-
  all_variables %>% 
  na.omit()

save(all_variables, all_completed, file = here::here("models/full_dataset.RData"), compress = "gzip")

# n = 67 
# seizures_any_clinical      n
# 1 FALSE                    28
# 2 TRUE                     39


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

source(here::here("models/isomiRs_PCA_clust.R"))

filter_values$isomiR <- filter_values$isomiR[filter_values$isomiR$feature_id %in% selected_isomiRs, ]

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
