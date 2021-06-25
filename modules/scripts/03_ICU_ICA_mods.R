rm(list = ls())
source("./misc/helper_functions.R")
source("./modules/scripts/ModuleFuncs.R")

### Read in expression data 
# Uncorrected
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

meta <- tr_te_dat$icu$meta 
expr <- tr_te_dat$icu$expr %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  as.matrix() %>% 
  DESeq2::varianceStabilizingTransformation() 
deconv_res <- deconv_res$icu$cibersort
all(colnames(expr) == meta$sample_identifier)

# Corrected 
expr_corr <- read_rds("./modules/results/ICU_expr_corr.RDS")
expr_corr <- expr_corr$expr_corr

# Combine 
all(colnames(expr) == colnames(expr_corr))
expr_mats <- list(expr_uncorr = expr, expr_corr = expr_corr)

### ICA Pipeline
# Perform ICA 
ICA_res <- expr_mats %>% 
  map(~perform_ica(.x, kurtosis_keep = 10, qval_filt = 0.001, seed = 1))

# Plot ICA Module Pathways.
ICA_mod_plot <- ICA_res %>% 
  map(~plot_pathwy_mods(.x))

# Get Eigenvectors
ICA_mod_eigen <- ICA_res %>% 
  map2(expr_mats, ~get_mod_eigenvect(.x, .y))

# Get module eigenvector-outcome associations. 
vars <- c("High_Low","HighInt_Low", "survive")
ICA_mod_eigen_clin <- list(
  expr_uncorr = lm_clin_eig(ICA_mod_eigen$expr_uncorr, meta,  vars_of_int = vars , deconv = deconv_res, incl_cell_props = TRUE),
  expr_corr = lm_clin_eig(ICA_mod_eigen$expr_corr, meta, vars_of_int = vars , incl_cell_props = FALSE)
)

# Get module signatures 
ICA_mod_sigs <- ICA_mod_eigen_clin %>% 
  map2(ICA_res, ~get_mod_sigs(.x, .y,  vars_of_int = vars, top_n_mods = 10, sig_qualifier = "ICU") )

# Save
ICA_res %>% write_rds("./modules/results/ICU_ICA_res.rds")
ICA_mod_plot %>% write_rds("./modules/results/ICU_ICA_mod_plot.rds")
ICA_mod_eigen %>% write_rds("./modules/results/ICU_ICA_mod_eigen.rds")
ICA_mod_eigen_clin %>% write_rds("./modules/results/ICU_ICA_mod_eigen_clin.rds")
ICA_mod_sigs %>% write_rds("./modules/results/ICU_ICA_mod_sigs.rds")

# Averge module size
ICA_res %>% 
  map(~map(.x$ICA_mods_filt, ~nrow(.x))) %>% 
  map(~flatten_int(.x)) %>% 
  map(~mean(.x))

# Module similarity 
ICA_jac_clust <- ICA_res %>% 
  map(~get_jac_mat(.x$ICA_mods_filt) %>% 
        t() %>% 
        vegan::vegdist(method = "jaccard", binary = TRUE, diag = TRUE) %>% 
        hclust %>% plot())


