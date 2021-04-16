rm(list = ls())
source("./helper_functions.R")

# Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")

expr <- tr_te_dat$er_tr$expr %>% dplyr::select(!contains("hc"))
meta <- tr_te_dat$er_tr$meta %>% filter(condition != "healthy_control")
all(colnames(expr)[-1] == meta$sample_identifier)

### VST data
expr_vst <- expr %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  as.matrix() %>% 
  DESeq2::varianceStabilizingTransformation()  %>% 
  sva::ComBat(batch = meta$sequencing_month_year) %>% 
  as.data.frame()

### Read in deconvolution data 
deconv_res <- read_rds("./deconvolution/deconv_res.rds")$er_tr$cibersort %>% 
  mutate(cell_type = janitor::make_clean_names(cell_type)) %>% 
  pivot_longer(cols = -cell_type, names_to = "sample_identifier") %>% 
  pivot_wider(id_cols = sample_identifier,names_from = cell_type ) %>% 
  filter(sample_identifier %in% meta$sample_identifier) %>% 
  remove_zero_var_cols() 
deconv_res_pca <- deconv_res %>% 
  column_to_rownames(var = "sample_identifier") %>% 
  perform_pca() 
deconv_res_pca_x <- deconv_res_pca$x %>% 
  dplyr::select(one_of("sample_identifier", paste0("PC", 1:12))) %>% 
  dplyr::rename_at(vars(matches("PC")), ~ paste0(., "_CiSo"))
meta_vars_to_corr <- meta %>% 
  dplyr::select(one_of("sample_identifier", "age", "gender")) %>%
  mutate(gender = ifelse(gender == 1, "F", "M")) %>% 
  left_join(deconv_res_pca_x) 

### Look at factors pre-correction
expr_pcs_lm_factors <- function(expr_norm, meta_w_select_vars) {
  
  # Perform PCA on expression data 
  expr_pre_corr_pca <- expr_norm %>% 
    t() %>% 
    perform_pca()
  
  # Create metadata with PCs and factors of interest
  expr_pre_corr_pca_pcs_meta <- expr_pre_corr_pca$x %>% 
    dplyr::select(one_of("sample_identifier", paste0("PC", 1:10))) %>% 
    left_join(meta_w_select_vars, by = "sample_identifier") %>% 
    column_to_rownames(var = "sample_identifier")
  
  # Regress PCs on factors 
  pca_contribs <- list()
  for (i in 1:10) {
    pca_contribs[[paste0("PC", i)]] <- 
      lm(
        formula(paste0("PC", i," ~ ", paste(colnames(meta_w_select_vars)[-1], collapse = " +"  ) )),
        data = expr_pre_corr_pca_pcs_meta
      ) %>%  broom::tidy() 
  }
  
  # Plot
  breaksList = seq(0, 10, by = .1)
  
  htmap <- pca_contribs %>% 
    bind_rows(.id = "PC") %>% 
    filter(term != "(Intercept)") %>% 
    mutate(log_p = -log10(p.value)) %>% 
    dplyr::select(one_of("PC", "term", "log_p")) %>% 
    pivot_wider(names_from = "PC", values_from = "log_p") %>% 
    column_to_rownames(var = "term") %>% 
    pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,  breaks = breaksList, 
             color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100) )
  
  return(htmap)
  
}

expr_pcs_lm_factors(expr_vst, meta_vars_to_corr)

#### Correct data 
# Create the model matrix of unwanted covariates
model.mat <- meta_vars_to_corr %>% 
  column_to_rownames(var = "sample_identifier") 
model.mat <- model.matrix( ~ .-1, model.mat)

# Write GLMs with Ridge penalty 
glms <- c(1:nrow(expr_vst)) %>% 
  map(~glmnet::cv.glmnet(model.mat, y = as.numeric(expr_vst[.x,]), alpha = 0, keep = FALSE))
names(glms) <- rownames(expr_vst)
glms_resid_mat <- c(1:nrow(expr_vst)) %>% 
  map(~as.numeric(expr_vst[.x,]) - predict(glms[[.x]], model.mat, s = 'lambda.1se')[,1]) %>% 
  map_df(~.x) %>% 
  set_rownames(rownames(expr_vst)) 

expr_pcs_lm_factors(glms_resid_mat, meta_vars_to_corr)
expr_pcs_lm_factors(glms_resid_mat,  dplyr::select(meta, one_of("sample_identifier", "sequencing_month_year")))

write_rds(glms_resid_mat, "./modules/results/expr_corr_resid.rds")
