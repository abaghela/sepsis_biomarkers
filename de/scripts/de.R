rm(list = ls())
source("./helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")

expr <- tr_te_dat$er_tr$expr %>% dplyr::select(!contains("hc"))
meta <- tr_te_dat$er_tr$meta %>% filter(condition != "healthy_control")
all(colnames(expr)[-1] == meta$sample_identifier)

### PCA
pca <- expr %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  as.matrix() %>% 
  DESeq2::varianceStabilizingTransformation()  %>% 
  sva::ComBat(batch = meta$sequencing_month_year) %>%
  t() %>% 
  perform_pca()

pca$x %>%
  left_join(meta) %>% 
  ggplot(aes(x = PC1, y = PC2,  fill = HighInt_Low_BC)) +
  geom_point(size = 4, colour="black", pch=21 ) + 
  xlab(paste0("PC1 (", pca$pov[1]*100, "%)" )) + 
  ylab(paste0("PC2 (", pca$pov[2]*100, "%)" )) + 
  theme_minimal() +
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13))

### Add in cell proportion information
deconv_res <- read_rds("./deconvolution/deconv_res.rds")

deconv_res_tr_df <- deconv_res$er_tr$cibersort %>% 
  dplyr::select(one_of("cell_type", meta$sample_identifier)) %>% 
  mutate(cell_type = janitor::make_clean_names(cell_type)) %>% 
  remove_zero_var_cols() %>% 
  column_to_rownames(var = "cell_type") %>% 
  t()
deconv_res_pca <- deconv_res_tr_df %>% perform_pca()
deconv_res_pca_x <- deconv_res_pca$x %>% dplyr::select(one_of("sample_identifier", paste0("PC", 1:explain_95(deconv_res_pca$pov))))

meta %<>% left_join(deconv_res_pca_x)

### DE
cols_to_keep <- c("sequencing_month_year", "age", "gender", "PC1")
outcomes <- c("survive","icu_adm","culture","sofa_sev_24","sofa_sev_prog", "HighInt_Low_BC", "HighInt_Low", "High_IntLow")

# Get an idea of patient numbers
for (i in outcomes){
  meta %>% group_by_at(i) %>% summarize(n=n(), .groups = 'drop') %>% print()
}

perform_de_pthwy_enr <- function(expr, meta, outcome) {
  
  # Print the outcome
  cat(paste0("Outcome: ", outcome, "\n"))
  
  # Set up meta and expression data
  met <- meta %>% 
    dplyr::select(one_of("sample_identifier", outcome, cols_to_keep)) %>% 
    dplyr::rename(comparison = outcome) %>% 
    filter(!is.na(comparison))
  
  exp <- expr %>% 
    dplyr::select(one_of("ensembl_gene_id", met$sample_identifier)) %>% 
    column_to_rownames(var = "ensembl_gene_id")
  
  # Check if samples are in the right order
  right_order <-  all(colnames(exp) == met$sample_identifier)
  
  if (!right_order) {
    cat("Things are not in the right order")
    stop()
  }
  
  # Perform differential expression
  de_res <- de(counts = exp, meta = met, FC = 1.5, 
               des = paste0("comparison + ", paste(cols_to_keep, collapse = " + ")), 
               main_covar = "comparison", filt = FALSE)
  
  # Perform pathway enrichment
  pthwy <- de_res %>% 
    map(~filter(.x, de == "de")) %>% 
    map(~pathway_enrichment(.x$entrezgene_id, p_val = 0.01, universe_list = universe$entrezgene_id))
  
  pthwy_up <- de_res %>% 
    map(~filter(.x, direction == "up")) %>% 
    map(~pathway_enrichment(.x$entrezgene_id, p_val = 0.01, universe_list = universe$entrezgene_id))
  
  pthwy_down <- de_res %>% 
    map(~filter(.x, direction == "down")) %>% 
    map(~pathway_enrichment(.x$entrezgene_id, p_val = 0.01, universe_list = universe$entrezgene_id))
  
  # Save Results
  res <- list()
  res[["de"]] <- de_res
  res[["pthwy"]] <- list(
    all = pthwy,
    up = pthwy_up,
    down = pthwy_down)
  return(res)
  }

res <- outcomes %>% 
  map(~perform_de_pthwy_enr(expr, meta, outcome = .x)) %>% 
  set_names(outcomes)

# Save results
res %>% write_rds(paste0("./de/results/de_tr_res.rds"))

# Get numbers
res %>% 
  map(~.x$de) %>% 
  map(~map(.x, ~filter(.x, de == "de"))) %>% 
  map(~map(.x, ~nrow(.x)))
res %>% 
  map(~.x$de) %>% 
  map(~map(.x, ~filter(.x, de == "de" & direction == "up"))) %>% 
  map(~map(.x, ~nrow(.x)))

# Read in results
res_de <- read_rds(paste0("./results/severity_combine_de_tr_icu_res.rds"))
res_pthwy <-  read_rds(paste0("./results/severity_combine_pthwy_enr_tr_icu_res.rds"))

### Pathway Enrichment plots
res_pthwy_df <- res_pthwy %>% 
  map(~map(.x, ~bind_rows(.x, .id = "comparison"))) %>% 
  map(~bind_rows(.x, .id = "direction")) %>% 
  bind_rows(.id = "outcome")
 
res_pthwy_df %>%
  plot_reactome_pathway()

### Volcano plot 
res_de$sofa_24_only_sev_group$high_low %>% vp(colour = "red", labs_to_incl = 10, fc_cutoff = 1.09, title = "High vs Low")
ggsave("./de/figures/de_high_vs_low.tiff", scale = 1.2)
ggsave("./de/figures/de_high_vs_low.png", scale = 1.2)
