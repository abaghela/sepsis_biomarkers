# Build Elastic Nets for mortality/severity . 

rm(list = ls())
source("./misc/helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

ICU_Train_Mod <- function(outcome) {
  
  # Set up data 
  meta <- tr_te_dat$icu$meta %>% 
    dplyr::select(one_of("sample_identifier", outcome)) %>% 
    na.omit()
  expr <- tr_te_dat$icu$expr %>% 
    dplyr::select(one_of("ensembl_gene_id", meta$sample_identifier)) %>%
    column_to_rownames(var = "ensembl_gene_id") %>%
    as.matrix() %>%
    varianceStabilizingTransformation() %>%
    as.data.frame() %>% 
    rownames_to_column(var = "ensembl_gene_id")
  
  all(colnames(expr)[-1] == meta$sample_identifier)
  
  # Read in genesets of interest & combine !HAVE TO ADD ICA SIGS!
  ICU_DE <- read_rds("./de/results/ICU_DE_res.rds") %>% 
    keep(names(.) %in% c("survive","High_Low_10_CutOff" )) %>% 
    map(~map(.x$de, ~filter(.x, de == "de")))
  
  ICU_DE_ALL <- ICU_DE %>% 
    map(~bind_rows(.x, .id = "comparison")) %>% 
    map(~.x$ensembl_gene_id) %>% 
    set_names(paste0("ICU_",names(.), "_DE_ALL"))
  ICU_DE_ALL[["ICU_combine_DE_ALL"]] <- ICU_DE_ALL %>% flatten_chr() %>% unique()
  ICU_DE_UP <- ICU_DE %>% 
    map(~bind_rows(.x, .id = "comparison")) %>% 
    map(~filter(.x, direction == "up")) %>% 
    map(~.x$ensembl_gene_id) %>% 
    set_names(paste0("ICU_",names(.), "_DE_UP"))
  ICU_DE_UP[["ICU_combine_DE_UP"]] <- ICU_DE_UP %>% flatten_chr() %>% unique()
  
  all_sigs <- append(ICU_DE_ALL, ICU_DE_UP)
  
  # Write F1 summary function 
  f1 <- function(data, lev = NULL, model = NULL) {
    f1_val <- MLmetrics::F1_Score(y_pred = data$pred,
                                  y_true = data$obs,
                                  positive = lev[1])
    c(F1 = f1_val)
  }
  
  # Enable parallel computing 
  cl <- makePSOCKcluster(16)
  registerDoParallel(cl)
  
  # Loop through each signature
  classifier_res <- list()
  for (sig in names(all_sigs)) {
    cat(paste0("\n ", sig, "\n"))
    
    # Filter expression data to only include the genes of interest
    expr_filt <- expr %>% 
      filter(ensembl_gene_id %in% all_sigs[[sig]] ) %>% 
      column_to_rownames(var = "ensembl_gene_id") 
    
    # Combine expression and metadata with outcome of interest
    expr_filt_t <- expr_filt %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "sample_identifier") %>% 
      left_join(meta, by = "sample_identifier") %>% 
      #dplyr::select(-one_of("sequencing_month_year")) %>% 
      column_to_rownames(var = "sample_identifier")
    
    all(meta$sample_identifier == rownames(expr_filt_t))
    
    # Define training options
    set.seed(1)
    fitControl <- trainControl(
      method = "repeatedcv",
      savePredictions = "all",
      number = 10,
      repeats = 10,
      classProbs = TRUE,
      #sampling = "down",
      summaryFunction =  f1,
      allowParallel = TRUE
    )
    
    # Elastic Net 
    classifier_res[[sig]] <- train(
      
      as.formula(paste0(outcome, " ~ .")), 
      data = expr_filt_t, 
      method = "glmnet",
      trControl = fitControl,
      preProc = c("center", "scale"),
      metric = "F1",
      tuneGrid = expand.grid(.alpha= seq(0,1, by = 0.1),  .lambda=seq(0, 10, by = 0.1))
      
    )
    
    }
  
  stopCluster(cl)
  
  return(classifier_res)
}

outcomes <- c("survive", "High_IntLow")
res <- outcomes %>% 
  map(~ICU_Train_Mod(outcome = .x)) %>% 
  set_names(outcomes)
res %>% map(~map(.x, ~.x$resample$F1 %>% na.omit() %>% mean()))
res %>% write_rds("./train_test_mods/results/ICU_Train_Mods.RDS")

res$survive$ICU_combine_DE_UP %>% extract_coefs() %>% View("first")
res$survive$ICU_survive_DE_UP %>% extract_coefs() %>% View("second")
classifier_res$resample$ROC %>% mean()
classifier_res$resample$Sens %>% mean()
classifier_res$resample$Spec %>% mean()
classifier_res$pred %>% filter(alpha == 1 & lambda == 0) %>% View()


