#### Script to measure enrichment of ICU signatures in SS Train and Test patients.
rm(list = ls())
source("./misc/helper_functions.R")

### Which signatures are worthy of exploration
SS_ROAST_Enr <- read_rds("./validation/results/SS_ROAST_Enrich.rds")
SS_ROAST_Enr %>% 
  map(~map(.x, ~filter(.x, PValue <= 0.05))) %>% 
  map(~map(.x, ~rownames(.x))) %>% 
  map(~unlist(.x, use.names = FALSE)) %>% 
  unlist(use.names = FALSE) %>% 
  unique() %>% 
  sort()

### Read in signatures 
all_sigs <- read_rds("./validation/results/DE_ICA_sigs.RDS")
all_sigs <- all_sigs[str_detect(names(all_sigs),"ICU")]

### Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))


SS_Sig_Model <- function() {
  
    expr <- tr_te_dat$er_tr$expr
    meta <- tr_te_dat$er_tr$meta %>% 
      dplyr::select(one_of("sample_identifier", "sequencing_month_year","HighInt_Low")) %>% 
      filter(!is.na(HighInt_Low))

    expr_filt <- expr %>%
      dplyr::select(one_of("ensembl_gene_id", meta$sample_identifier)) %>%
      column_to_rownames(var = "ensembl_gene_id") %>%
      as.matrix() %>%
      varianceStabilizingTransformation() %>%
     # sva::ComBat(batch = meta$sequencing_month_year) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "ensembl_gene_id") %>% 
      filter(ensembl_gene_id %in% all_sigs$ICU_survive_DE_UP) %>% 
      column_to_rownames(var = "ensembl_gene_id") 
    
    expr_filt_t <- expr_filt %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "sample_identifier") %>% 
      left_join(meta) %>% 
      dplyr::select(-one_of("sequencing_month_year")) %>% 
      column_to_rownames(var = "sample_identifier")
    
    all(meta$sample_identifier == rownames(expr_filt_t))
  

  set.seed(1)
  fitControl <- trainControl(
    method = "repeatedcv",
    savePredictions = "all",
    number = 10,
    repeats = 10,
    classProbs = TRUE,
    #sampling = "down",
    summaryFunction =  twoClassSummary,
    allowParallel = TRUE
  )
  
  cl <- makePSOCKcluster(16)
  registerDoParallel(cl)
  classifier_res <- train(HighInt_Low ~ ., data = expr_filt_t, 
                                          method = "glmnet",
                                          trControl = fitControl,
                                          preProc = c("center", "scale"),
                                          #family = "multinomial",
                                          #type.multinomial = "grouped",
                                          #standardize = TRUE,
                                          #metric= "AUC",
                                          #metric = "ROC",
                                          tuneGrid = expand.grid(.alpha= 1,  .lambda=seq(0, 10, by = .1)))
  stopCluster(cl)
  
  
  auc_sens_spec <- function(mod) {
    y_len = length(mod$levels)
    if (y_len == 2) {
      res = list(auc = mod$resample$ROC %>% mean(),
                 #f1 = mod$resample$Mean_F1 %>% na.omit() %>% mean(),
                 sens = mod$resample$Sens %>% mean(),
                 spec = mod$resample$Spec %>% mean())
    } else {
      res = list(auc = mod$resample$AUC %>% mean(),
                 f1 = mod$resample$Mean_F1 %>% na.omit() %>% mean(),
                 sens = mod$resample$Mean_Sensitivity %>% mean(),
                 spec = mod$resample$Mean_Specificity %>% mean())
    }
    
    return(res)
  }
  
  classifier_res %>% 
    auc_sens_spec()  
  
  classifier_res %>% 
    extract_coefs()
    
}











SS_Sig_Cluster <- function(){

  expr <- tr_te_dat$er_tr$expr
  meta <- tr_te_dat$er_tr$meta %>%
    filter(!is.na(sofa_sev_24))

  expr_filt <- expr %>%
    dplyr::select(one_of("ensembl_gene_id", meta$sample_identifier)) %>%
    column_to_rownames(var = "ensembl_gene_id") %>%
    as.matrix() %>%
    varianceStabilizingTransformation() %>%
    sva::ComBat(batch = meta$sequencing_month_year) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "ensembl_gene_id") %>% 
    filter(ensembl_gene_id %in% all_sigs$ICU_survive_DE_UP) %>% 
    column_to_rownames(var = "ensembl_gene_id") %>% 
    t() %>% scale(center = TRUE, scale = TRUE) %>%  t() %>%
    as.data.frame()

  all(meta$sample_identifier == colnames(expr_filt))

  column_ha <- ComplexHeatmap::HeatmapAnnotation(Severity = meta$sofa_sev_24)

  hey <- expr_filt %>%
    ComplexHeatmap::Heatmap(column_split = meta$sofa_sev_24,
                            name = "Z-score",
                            top_annotation = column_ha,
                            column_labels = FALSE,
                            show_column_names = FALSE,
                            show_column_dend = FALSE,
                            show_row_dend = FALSE,
                            cluster_columns = TRUE,
                            cluster_rows = TRUE,
                            cluster_column_slices = FALSE,
                            border = TRUE)

  hey
}
