rm(list = ls())

library(magrittr)
library(tidyverse)
library(functionjunction)
library(RColorBrewer)
library(alluvial)
library(viridis)
library(MLmetrics)
source("./helper_functions.R")

### Read in data
tr_te_dat <- read_rds("../paper_er_icu_sepsis/data_overview/tr_te_dat.RDS")
universe <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_221020.RDS")$universe

meta <- tr_te_dat$test$meta 
expr <- tr_te_dat$test$expr$raw 

meta_no_hc <- tr_te_dat$test$meta %>% filter(!condition == "healthy_control")
expr_no_hc <- tr_te_dat$test$expr$raw %>% dplyr::select(one_of("ensembl_gene_id", meta_no_hc$sample_identifier))

### Create outcomes
meta_no_hc %<>% 
  dplyr::select(one_of("sample_identifier", "condition", "sequencing_month_year", 
                       "first_at_ed_sofa", "worst_within_72_sofa","outcome_mortality", 
                       "sofa_24_only_sev_group","micro_blood_culture_pathogen")) %>% 
  dplyr::rename(sofa = "first_at_ed_sofa", survive = "outcome_mortality") %>% 
  mutate(sofa_progress = ifelse(worst_within_72_sofa - sofa > 0, "worse", 
                                ifelse(worst_within_72_sofa == sofa, "same", "improve"))) %>% 
  mutate(survive = ifelse(survive == 0, "survived", "dead")) %>% 
  mutate(sofa_24_only_sev_group = ifelse(is.na(sofa_24_only_sev_group), NA, .$sofa_24_only_sev_group)) %>% 
  mutate(sofa_24_only_sev_group = factor(sofa_24_only_sev_group, levels = c("low", "intermediate", "high")))

meta_no_hc %<>% 
  mutate(High_Low = ifelse(sofa_24_only_sev_group == "high", "High",
                           ifelse(sofa_24_only_sev_group == "low", "Low", NA))) %>%
  mutate(High_Low = factor(High_Low, levels = c("Low", "High"))) %>% 
  
  mutate(High_Int = ifelse(sofa_24_only_sev_group == "high", "High",
                           ifelse(sofa_24_only_sev_group == "intermediate", "Intermediate", NA))) %>%
  mutate(High_Int = factor(High_Int, levels = c("Intermediate", "High"))) %>% 
  
  mutate(HighInt_Low = ifelse(is.na(sofa_24_only_sev_group) , NA, 
                              ifelse(sofa_24_only_sev_group %in% c("high", "intermediate"), "HighInt", "Low"))) %>% 
  mutate(HighInt_Low = factor(HighInt_Low, levels = c("Low", "HighInt"))) %>% 
  mutate(High_IntLow = ifelse(is.na(sofa_24_only_sev_group) , NA, 
                              ifelse(sofa_24_only_sev_group %in% c("high"), "High", "IntLow"))) %>% 
  mutate(High_IntLow = factor(High_IntLow, levels = c("IntLow", "High"))) %>%
  
  mutate(Worse_Improve = ifelse(sofa_progress == "worse", "Worse", 
                                ifelse(sofa_progress == "improve", "Improve", NA))) %>% 
  mutate(Worse_Improve = factor(Worse_Improve, levels = c("Improve", "Worse"))) %>% 
  mutate(High_Low_BC = ifelse(is.na(sofa_24_only_sev_group), NA,
                              ifelse(sofa_24_only_sev_group == "high" & micro_blood_culture_pathogen == 1, "High_BC",
                                     ifelse(sofa_24_only_sev_group == "low" & micro_blood_culture_pathogen == 0, "Low_BC",
                                            ifelse(sofa_24_only_sev_group == "intermediate", NA, NA))))) %>% 
  mutate(High_Low_BC = factor(High_Low_BC, levels = c("Low_BC", "High_BC"))) %>% 
  mutate(survive = factor(survive, levels = c("survived", "dead")))

# Read in classifiers 
classifiers <- read_rds("./results/classifiers_res_tr_icu.rds")
genes_to_use <- read_rds("./results/severity_signatures.rds")

# Predict and batch correct based on outcome
cols_to_keep <- c("sequencing_month_year")
#outcomes <- c("survive", "Worse_Improve","High_Low","High_Int", "HighInt_Low", "High_IntLow")
outcomes <- c("High_Low", "HighInt_Low")

# Create AUC function
get_AUC_sens_spec <- function(df, pred.prob.name) {
  res <- data.frame(ROC = AUC(df[[pred.prob.name]], df$comparison_1_0 ),
                    Sens = Sensitivity(df$comparison, df$pred_name ),
                    Spec = Specificity(df$comparison, df$pred_name)
                    # F1_S = F1_Score(df$comparison, df$pred_name)
  )
  return(res)
}

pred_res <- list()
for (i in outcomes) {
  # Create metadata & normalize data & batch correct
  met <- meta_no_hc %>% 
    dplyr::select(one_of("sample_identifier", i, cols_to_keep)) %>% 
    dplyr::rename(comparison = i) %>% 
    filter(!is.na(comparison)) 
  exp <- expr_no_hc  %>% 
    dplyr::select(one_of("ensembl_gene_id", met$sample_identifier)) %>% 
    column_to_rownames(var = "ensembl_gene_id") %>% 
    as.matrix() %>% 
    DESeq2::varianceStabilizingTransformation() %>%
    sva::ComBat(batch = met$sequencing_month_year,  mod = model.matrix(~comparison, data=met)) %>%
    t() %>% 
    as.data.frame()
  
  right_order <-  all(rownames(exp) == met$sample_identifier)
  if (!right_order) {
    cat("Things are not in the right order")
    stop()
  }
  
  # Predict, get AUCs
  levs <- levels(met$comparison)
  
  print(i)
  for (gene_input in names(classifiers[[i]])) {
    for (mod in names(classifiers[[i]][[gene_input]]) ) {
      # Filt 
      exp_filt <- exp %>% 
        rownames_to_column(var = "sample_identifier") %>% 
        dplyr::select(one_of("sample_identifier", genes_to_use[[gene_input]])) %>% 
        column_to_rownames(var = "sample_identifier")  %>% 
        scale(center = TRUE, scale = TRUE)
      
      # Predict  
      pred_df <- classifiers[[i]][[gene_input]][[mod]] %>% predict(exp_filt, type = "prob")
      pred_df <- data.frame(sample_identifier = rownames(exp_filt), pred = pred_df, model_used = mod, genes_used = gene_input) %>%
        left_join(met) %>%
        mutate(pred_name = ifelse(!!sym(paste0("pred.", levs[2])) >= 0.5, levs[2], levs[1]),
               comparison_1_0 = ifelse(comparison == levs[2], 1, 0))
      
      if(length(unique(pred_df$pred_name)) == 1) {
        next()}
      
      pred_df <- pred_df %>% get_AUC_sens_spec(pred.prob.name = paste0("pred.", levs[2]))
      
      pred_res[[i]][[gene_input]][[mod]] <- pred_df
      
      
    }
  }
}

combine_all_tr <- classifiers %>% 
  map(~map(.x, ~map(.x, ~data.frame(ROC = mean(.x$resample$ROC), Sens = mean(.x$resample$Sens), Spec = mean(.x$resample$Spec)) ))) %>% 
  map(~map(.x, ~bind_rows(.x, .id = "algo"))) %>% 
  map(~bind_rows(.x, .id = "gene_input")) %>% 
  bind_rows(.id = "outcome") %>% 
  mutate(tr_te = "train")

combine_all_te <- pred_res %>% 
  map(~map(.x, ~bind_rows(.x, .id = "algo"))) %>% 
  map(~bind_rows(.x, .id = "gene_input")) %>% 
  bind_rows(.id = "outcome") %>% 
  mutate(tr_te = "test")

combine_all_tr %>% 
  bind_rows(combine_all_te) %>% 
  mutate(ROC = round(ROC, 2), Sens = round(Sens, 2), Spec = round(Spec, 2)) %>% 
  unite(col = "results", ROC, Sens, Spec, sep = "; ", remove = TRUE) %>% 
  pivot_wider(names_from = tr_te, values_from = results) %>% 
  #filter(outcome == "High_Low") %>% 
  # filter(algo %in% c("lm", "elastic", "lasso")) %>% 
  View("hi")
