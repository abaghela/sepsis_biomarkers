rm(list = ls())
source("./misc/helper_functions.R")

# Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

### Look at factors pre-correction and also correct data if desired
expr_pcs_lm_factors <- function(expr, meta, deconv = NULL,  incl_cell_props = FALSE, correct_expr = FALSE) {
  
  # Metadata - Keep variables of interest 
  met <- meta %>% 
    dplyr::select(one_of("sample_identifier", "sequencing_month_year", "age", "gender"))
  
  # VST data
  expr_vst <- expr %>% 
    column_to_rownames(var = "ensembl_gene_id") %>% 
    as.matrix() %>% 
    DESeq2::varianceStabilizingTransformation()
  
  # Batch correct if more than one batch
  if (length(unique(met$sequencing_month_year)) > 1) {
    cat("\n Batch Correction being performed. \n")
    expr_vst <- expr_vst %>% 
      sva::ComBat(batch = met$sequencing_month_year) %>% 
      as.data.frame()
  }
  
  # Remove low variance genes - GLM's will not fit if response is near constant 
  gene_vars <- expr_vst %>% apply(1, function(x) {var(x)})
  expr_vst <- expr_vst[-c(which(gene_vars < 0.01)),]
  
  # Perform PCA on expression data 
  expr_pre_corr_pca <- expr_vst %>% 
    t() %>% 
    perform_pca()
  
  # Perform PCA on deconvolution data 
  if (incl_cell_props) {
    deconv_pca <- deconv %>% 
      rownames_to_column(var = "sample_identifier") %>% 
      filter(sample_identifier %in% meta$sample_identifier) %>% 
      remove_zero_var_cols() %>% 
      column_to_rownames(var = "sample_identifier") %>% 
      perform_pca() 
    
    pcs_to_incl <- deconv_pca$pov %>% explain_95(perc = 90)
    cat(paste0("Cell Proportion PCs used: ", pcs_to_incl, "\n"))
    pcs_to_incl <- paste0("PC", 1:pcs_to_incl)
    
    deconv_pca_x <- deconv_pca$x %>% 
      dplyr::select(one_of("sample_identifier", pcs_to_incl)) %>% 
      dplyr::rename_at(vars(matches("PC")), ~ paste0("CiSo_", .))
    
    met <- met %>% 
      left_join(deconv_pca_x,  by = "sample_identifier")
  
  }

  # Create metadata with PCs and factors of interest
  expr_pre_corr_pca_pcs_meta <- expr_pre_corr_pca$x %>% 
    dplyr::select(one_of("sample_identifier", paste0("PC", 1:10))) %>% 
    left_join(met, by = "sample_identifier") %>% 
    column_to_rownames(var = "sample_identifier")
  
  # Regress PCs on factors 
  pca_contribs <- list()
  for (i in 1:10) {
    pca_contribs[[paste0("PC", i)]] <- 
      lm(
        formula(paste0("PC", i," ~ ", paste(colnames(met)[-c(1,2)], collapse = " +"  ) )),
        data = expr_pre_corr_pca_pcs_meta
      ) %>%  broom::tidy() 
  }
  
  # Plot
  breaksList = seq(0, 10, by = .1)
  
  htmap_pre_corr <- pca_contribs %>% 
    bind_rows(.id = "PC") %>% 
    filter(term != "(Intercept)") %>% 
    mutate(log_p = -log10(p.value)) %>% 
    dplyr::select(one_of("PC", "term", "log_p")) %>% 
    pivot_wider(names_from = "PC", values_from = "log_p") %>% 
    column_to_rownames(var = "term") %>% 
    pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,  breaks = breaksList, 
             color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100) )
  
  #return(htmap_pre_corr)
  
  if (correct_expr) {
    cat("\n Correcting expression matrix with Ridge regression. \n")
    # Create the model matrix of unwanted covariates
    model.mat <- met %>% 
      dplyr::select(-one_of("sequencing_month_year")) %>% 
      column_to_rownames(var = "sample_identifier") 
    model.mat <- model.matrix( ~ .-1, model.mat)
    
    # Write GLMs with Ridge penalty 
    set.seed(1)
    glms <- c(1:nrow(expr_vst)) %>% 
      map(~glmnet::cv.glmnet(model.mat, y = as.numeric(expr_vst[.x,]), alpha = 0, nfolds = 5, keep = FALSE))
    names(glms) <- rownames(expr_vst)
    glms_resid_mat <- c(1:nrow(expr_vst)) %>% 
      map(~as.numeric(expr_vst[.x,]) - predict(glms[[.x]], model.mat, s = 'lambda.1se')[,1]) %>% 
      map_df(~.x) %>% 
      set_rownames(rownames(expr_vst)) 
    
    # Perform PCA on new matrix 
    expr_post_corr_pca <- glms_resid_mat %>% 
      t() %>% 
      perform_pca()
    
    # Create metadata with PCs and factors of interest
    expr_post_corr_pca_pcs_meta <- expr_post_corr_pca$x %>% 
      dplyr::select(one_of("sample_identifier", paste0("PC", 1:10))) %>% 
      left_join(met, by = "sample_identifier") %>% 
      column_to_rownames(var = "sample_identifier")
    
    # Regress PCs on factors 
    post_pca_contribs <- list()
    for (i in 1:10) {
      post_pca_contribs[[paste0("PC", i)]] <- 
        lm(
          formula(paste0("PC", i," ~ ", paste(colnames(met)[-c(1,2)], collapse = " +"  ) )),
          data = expr_post_corr_pca_pcs_meta
        ) %>%  broom::tidy() 
    }
    
    # Plot
    breaksList = seq(0, 10, by = .1)
    
    htmap_post_corr <- post_pca_contribs %>% 
      bind_rows(.id = "PC") %>% 
      filter(term != "(Intercept)") %>% 
      mutate(log_p = -log10(p.value)) %>% 
      dplyr::select(one_of("PC", "term", "log_p")) %>% 
      pivot_wider(names_from = "PC", values_from = "log_p") %>% 
      column_to_rownames(var = "term") %>% 
      pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,  breaks = breaksList, 
               color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100) )
    
    
    return(list(htmap_pre_corr = htmap_pre_corr, expr_corr = glms_resid_mat, htmap_post_corr = htmap_post_corr))
  }
  
return(list(htmap_pre_corr = htmap_pre_corr))
  
}

### SS Patients
expr <- tr_te_dat$er_tr$expr
meta <- tr_te_dat$er_tr$meta 
deconv <- deconv_res$er_tr$cibersort
all(colnames(expr)[-1] == meta$sample_identifier)

SS_corr <- expr_pcs_lm_factors(expr, meta, deconv, incl_cell_props = TRUE, correct_expr = TRUE)
SS_corr %>% write_rds("./modules/results/SS_expr_corr.RDS")

### ICU Patients
expr <- tr_te_dat$icu$expr
meta <- tr_te_dat$icu$meta 
deconv <- deconv_res$icu$cibersort
all(colnames(expr)[-1] == meta$sample_identifier)

ICU_corr <- expr_pcs_lm_factors(expr, meta, deconv, incl_cell_props = TRUE, correct_expr = TRUE)
ICU_corr %>% write_rds("./modules/results/ICU_expr_corr.RDS")
