source("./helper_functions.R")
library(immunedeconv)

# Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS") 

# Specify CIBERSORT script
set_cibersort_binary("../deconvolution/CIBERSORT.R")
set_cibersort_mat("../deconvolution/LM22.txt")
lm22 <- read_tsv("../deconvolution/LM22.txt")

deconv_res <- list()
for (i in names(tr_te_dat)) {
  expr_mat <- tr_te_dat[[i]]$expr %>% 
    left_join(universe) %>% 
    dplyr::select(!one_of("ensembl_gene_id", "entrezgene_id", "description")) %>% 
    filter(!is.na(hgnc_symbol)) %>% 
    filter(!hgnc_symbol %in% c("", "POLR2J3", "POLR2J4", "TBCE")) %>% 
    column_to_rownames(var = "hgnc_symbol") %>% 
    as.matrix() %>% 
    DESeq2::varianceStabilizingTransformation() 
  
  # Batch correct if required
  
  batches = tr_te_dat[[i]]$meta$sequencing_month_year %>% unique() %>% length()
  
  if (batches >= 2) {
    expr_mat =  sva::ComBat(expr_mat, tr_te_dat[[i]]$meta$sequencing_month_year)
  }
   
  
  res <- list(
    #mcp_counter = deconvolute(expr_mat, "mcp_counter"),
    #epic = deconvolute(expr_mat, "epic"), 
    #quantiseq = deconvolute(expr_mat, "quantiseq"), 
    #xcell = deconvolute(expr_mat, "xcell"),
    cibersort = deconvolute(expr_mat, "cibersort")
    #cibersort_abs = deconvolute(expr_mat, "cibersort_abs")
  )
  
  # Make nice 
  res %<>% 
    map(~mutate(.x, cell_type = janitor::make_clean_names(cell_type) )) %>% 
    map(~remove_zero_var_cols(.x)) %>% 
    map(~column_to_rownames(.x, var = "cell_type")) %>% 
    map(~t(.x))

  deconv_res[[i]] <- res
}

deconv_res %>% write_rds("./create_tr_te/deconv_res.rds")

system("mv ./CIBERSORT-Results.txt ./create_tr_te/")
