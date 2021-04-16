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
  DESeq2::varianceStabilizingTransformation() %>% 
  sva::ComBat(batch = meta$sequencing_month_year) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id")

### Read in Cell proportion corrected data 
expr_cor <- read_rds("./modules/results/expr_corr_resid.rds") %>% 
  rownames_to_column(var = "ensembl_gene_id")

all(colnames(expr_vst) == colnames(expr_cor))

### Combine 
expr_mods <- list(expr_vst = expr_vst, expr_cor = expr_cor)

# Read in components
ICA_res <- read_rds("./modules/results/ICA_res.rds")

# Get module eigenvectors
get_mod_eigenvect <- function(ICA_mods, expr_mat) {
  
  eig <- ICA_mods$ICA_mods_filt %>% 
    discard(~nrow(.x) == 0) %>% # For some reason there are modules with no genes in them... need to check this out. 
    map(~filter(expr_mat, ensembl_gene_id %in% .x$ensembl_gene_id)) %>% 
    map(~column_to_rownames(.x, var = "ensembl_gene_id")) %>% 
    map(~t(.x) %>% perform_pca() %>% pluck("x")) %>% 
    map(~dplyr::select(.x, one_of("sample_identifier", "PC1")))
  return(eig)

  }

ICA_mod_eigen <- ICA_res %>% map2(expr_mods, ~get_mod_eigenvect(.x, .y))

ICA_mod_eigen %>% write_rds("./modules/results/ICA_mod_eigen.rds")
