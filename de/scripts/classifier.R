rm(list = ls())

library(tidyverse)
library(magrittr)
library(functionjunction)
library(DESeq2)
library(viridis)
library(ggrepel)
library(UpSetR)
library(caret)
library(doParallel)
source("./helper_functions.R")

### Read in data
data <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_261120.RDS")
universe <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_261120.RDS")$universe
tr_te_dat <- read_rds("../paper_er_icu_sepsis/data_overview/tr_te_dat.RDS")
icu_dat <- read_rds("../paper_er_icu_sepsis/data_overview/icu_fir_sec_dat.RDS")

### Severity measures: "sofa_24_sev_group", "qsofa_sev_group", "sofa_72_sev_group"
sev_model <- data.frame(
  model = c("severity", "severity", "severity", "severity"),
  group = rev(c("high", "intermediate", "low", "healthy_control")),
  color = rev(c("red", "darkorange1" ,"burlywood4", "black")),
  proper_name = rev(c("High", "Int", "Low", "HC"))
)

### Set up data
expr <- tr_te_dat$train$expr$raw %>% 
  left_join(icu_dat$first$expr$raw)

keep <- c("sample_identifier", "condition", "sequencing_month_year",  "outcome_icu_survival",  "outcome_mortality",
          "first_at_ed_sofa","worst_within_72_sofa", "outcome_sofa_second", "outcome_sofa_first", "outcome_sofa_third",
          "age", "gender", "micro_blood_culture_pathogen")
meta <- icu_dat$first$meta %>% 
  dplyr::select(one_of(keep)) %>% 
  dplyr::rename(sofa = "outcome_sofa_second", survive = "outcome_icu_survival") %>% 
  mutate(sofa_progress = ifelse(outcome_sofa_third - outcome_sofa_first > 0, "worse",
                                ifelse(outcome_sofa_third == outcome_sofa_first, "same", "improve"))) %>% 
  mutate(survive = ifelse(survive == 1, "survived", "dead"))
meta <- tr_te_dat$train$meta %>% 
  dplyr::select(one_of(keep)) %>% 
  dplyr::rename(sofa = "first_at_ed_sofa", survive = "outcome_mortality") %>% 
  mutate(sofa_progress = ifelse(worst_within_72_sofa - sofa > 0, "worse", 
                                ifelse(worst_within_72_sofa == sofa, "same", "improve"))) %>% 
  mutate(survive = ifelse(survive == 0, "survived", "dead")) %>% 
  bind_rows(meta) 

### Remove controls
meta_no_hc <- meta %>% filter(!condition == "healthy_control")
expr_no_hc <- expr %>% dplyr::select(one_of("ensembl_gene_id", meta_no_hc$sample_identifier))
all(meta_no_hc$sample_identifier == colnames(expr_no_hc)[-1])

### Get severity groups
meta_no_hc %<>% 
  mutate(sofa_24_only_sev_group = 
           ifelse(condition == "healthy_control", "healthy_control", 
                  ifelse( sofa >= 5, "high", 
                          ifelse(sofa >= 2, "intermediate","low")))) %>% 
  mutate(sofa_24_only_sev_group = ifelse(is.na(sofa_24_only_sev_group), NA, .$sofa_24_only_sev_group)) %>% 
  mutate(sofa_24_only_sev_group = factor(sofa_24_only_sev_group, levels = c("low", "intermediate", "high")))

### Now make outcomes to model for classifiers. HighInt vs Low, High vs IntLow, Died vs Survive.
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

all(colnames(expr_no_hc)[-1] == meta_no_hc$sample_identifier)

### Read in signatures
genes_to_use <- read_rds("./results/severity_signatures.rds")

#### Start Classifiers
cols_to_keep <- c("sequencing_month_year")
#outcomes <- c("survive", "Worse_Improve","High_Low","High_Int", "HighInt_Low", "High_IntLow")
outcomes <- c("High_Low", "HighInt_Low")

cl <- makePSOCKcluster(16)
registerDoParallel(cl)
classifiers <- list() 

for (i in outcomes) {
  cat(paste0("Outcome: ", i, "\n"))
  
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
    sva::ComBat(batch = met$sequencing_month_year,  mod = model.matrix(~as.factor(comparison), data=met)) %>% 
    as.data.frame() 
  
  # Only select gene sets which matter to the specific comparison
  levs <- levels(unique(met$comparison))
  # geneset_to_try <- c(tolower(paste0(levs[2], "_", levs[1])), paste0(tolower(paste0(levs[2], "_", levs[1])), "_POS"),
  #                     paste0(levs[2], "_", levs[1]), paste0(levs[2], "_", levs[1], "_POS"),
  #                     "cr_de", "cr_combine, "endotox", "inflam", "endotox31", "patent_1",  "patent_2",  "patent_3",  "patent_4",  "patent_5")
  # geneset_to_try <- c(tolower(paste0(levs[2], "_", levs[1])), paste0(tolower(paste0(levs[2], "_", levs[1])), "_POS"),
  #                     paste0(levs[2], "_", levs[1]), paste0(levs[2], "_", levs[1], "_POS"))
  # geneset_to_try <- names(genes_to_use) %>% str_subset(paste(tolower(levs), collapse = "|"))
  # geneset_to_try <- c(geneset_to_try, c("cr_combine", "cr_de", "endotox31", "multinom_sig"))
  geneset_to_try <- names(genes_to_use)
  
  for (genes in names(genes_to_use)   ) {
    
    
    if (genes %in% geneset_to_try ) {
      
      cat(paste0(genes, "\n"))
      
      # Filter genes 
      exp_filt <- exp %>% 
        t() %>% 
        as.data.frame() %>%         
        rownames_to_column(var = "sample_identifier") %>% 
        dplyr::select(one_of("sample_identifier", genes_to_use[[genes]] )) %>% 
        column_to_rownames(var = "sample_identifier") %>% 
        scale(center = TRUE, scale = TRUE) %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "sample_identifier")
      
      # Check if expression and metadata are in the right order.
      right_order <-  all(exp_filt$sample_identifier == met$sample_identifier)
      if (!right_order) {
        cat("Things are not in the right order")
        stop()
      }
      
      # Create input matrix
      input <- met %>% dplyr::select(one_of("sample_identifier", "comparison")) %>% 
        left_join(exp_filt) %>% 
        column_to_rownames(var = "sample_identifier")
      
      # Create fitControl
      set.seed(1)
      fitControl <- trainControl(
        method = "repeatedcv",
        savePredictions = "all",
        number = 10,
        repeats = 10,
        classProbs = TRUE,
        summaryFunction =  twoClassSummary,
        sampling = "down",
        allowParallel = TRUE
      )
      
      # Set up classifiers 
      
      # print("sparseLDA")
      # set.seed(1)
      # classifiers[[i]][[genes]][["lasso"]] <- train(comparison ~ ., data = input ,
      #                                               method = "sparseLDA",
      #                                               trControl = fitControl,
      #                                               preProc = c("center", "scale"),
      #                                               metric= "ROC",
      #                                               tuneGrid = expand.grid(
      #                                                 NumVars= c(5, 10, 25, 50, 100),
      #                                                 lambda=seq(0, 10, by = .1)))
      
      
      
      print("lasso")
      set.seed(1)
      classifiers[[i]][[genes]][["lasso"]] <- train(comparison ~ ., data = input ,
                                                    method = "glmnet",
                                                    trControl = fitControl,
                                                    preProc = c("center", "scale"),
                                                    metric= "ROC",
                                                    tuneGrid = expand.grid(
                                                      .alpha= 1,
                                                      .lambda=seq(0, 10, by = .1)))
      print("elastic")
      set.seed(1)
      classifiers[[i]][[genes]][["elastic"]] <- train(comparison ~ ., data = input ,
                                                      method = "glmnet",
                                                      trControl = fitControl,
                                                      preProc = c("center", "scale"),
                                                      metric= "ROC",
                                                      tuneGrid = expand.grid(
                                                        .alpha= seq(0,1, by = .1),
                                                        .lambda=seq(0, 10, by = .1)))
      # print("lm")
      # set.seed(1)
      # classifiers[[i]][[genes]][["lm"]] <- train(comparison ~ ., data = input ,
      #                                            method = "glmnet",
      #                                            trControl = fitControl,
      #                                            preProc = c("center", "scale"),
      #                                            metric= "ROC",
      #                                            tuneGrid = expand.grid(
      #                                              .alpha= 1,
      #                                              .lambda= 0))
      
      print("rf")
      set.seed(1)
      rf_params <- c(ceiling(sqrt(ncol(input)))) %>% unique()
      classifiers[[i]][[genes]][["rf"]]  <- train(comparison ~ ., data = input ,
                                                  method = "rf",
                                                  trControl = fitControl,
                                                  preProc = c("center", "scale"),
                                                  metric= "ROC",
                                                  tuneGrid = expand.grid(mtry = rf_params))
      print("svm")
      set.seed(1)
      classifiers[[i]][[genes]][["svm"]] <- train(comparison ~ ., data = input ,
                                                  method = "svmRadialCost",
                                                  trControl = fitControl,
                                                  preProc = c("center", "scale"),
                                                  metric= "ROC",
                                                  tuneGrid = expand.grid(C = c(0.1, 0.2, 0.5, 1, 1.5, 2, 5, 10, 20)))
      
    }
  }
}

stopCluster(cl)

# Save results
classifiers %>% write_rds(paste0("./results/classifiers_res_tr_icu.rds"))
classifiers <- read_rds("./results/classifiers_res_tr_icu.rds")

combine_all <- classifiers %>% 
  map(~map(.x, ~map(.x, ~data.frame(ROC = mean(.x$resample$ROC), Sens = mean(.x$resample$Sens), Spec = mean(.x$resample$Spec))))) %>% 
  map(~map(.x, ~bind_rows(.x, .id = "algo"))) %>% 
  map(~bind_rows(.x, .id = "gene_input")) %>% 
  bind_rows(.id = "outcome")
combine_all %>% View()

# Look at coefficients 
sev_3_gene_sig <- classifiers$HighInt_Low$high_low$lasso %>% 
  extract_coefs() %>% 
  pull(ensembl_gene_id)
sev_3_gene_sig

expr_no_hc_vst <- expr_no_hc %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  as.matrix() %>% 
  DESeq2::varianceStabilizingTransformation() %>%
  sva::ComBat(batch = meta_no_hc$sequencing_month_year) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id")

expr_no_hc_vst %>% 
  filter(ensembl_gene_id %in% sev_3_gene_sig) %>% 
  pivot_longer(-"ensembl_gene_id", "sample_identifier") %>% 
  left_join(meta_no_hc) %>% 
  filter(!is.na(sofa_24_only_sev_group)) %>% 
  left_join(universe) %>% 
  ggplot(aes(x = sofa_24_only_sev_group, y = value, color = sofa_24_only_sev_group)) + 
  geom_boxplot(outlier.shape=NA, alpha = 0.5, lwd = 0.5) + 
  geom_point(position = "jitter", alpha = 0.2) +
  facet_grid(cols = vars(hgnc_symbol), scales = "free") + ylab("") + xlab("") +
  scale_color_manual("", labels =sev_model$proper_name[-1] , values = as.character(sev_model$color[-1] ) )
  
ggsave("./de/figures/de_high_vs_low_lasso_sig.tiff", scale = 1.1, height = 4, width = 5.5)
ggsave("./de/figures/de_high_vs_low_lasso_sig.png", scale = 1.1, height = 4, width = 5.5)
classifiers$HighInt_Low$high_low$lasso %>% 
  extract_coefs() %>% 
  write_csv("./de/results/classifier_High_Low_6_gene.csv")
