rm(list = ls())
source("./misc/helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

expr <- tr_te_dat$er_tr$expr
meta <- tr_te_dat$er_tr$meta 
all(colnames(expr)[-1] == meta$sample_identifier)

### Add in Endotype classification
endotype_res <- read_rds("../paper_er_icu_sepsis/endotype/consensus_clust_res/kmed_endotype_all_res_clusters.rds")
endotype_res <- endotype_res$top_5_percent$k_5 %>% 
  dplyr::select(one_of("sample_identifier", "cluster")) %>% 
  mutate(endotype2class = ifelse(cluster %in% c("cluster_2", "cluster_3"), "NPS_INF", "IHD_IFN_ADA")) %>% 
  mutate(endotype2class = factor(endotype2class, levels = c("IHD_IFN_ADA", "NPS_INF")))
meta %<>% left_join(endotype_res, by = "sample_identifier")

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
  ggplot(aes(x = PC1, y = PC2,  fill = qsofa_High_Low)) +
  geom_point(size = 4, colour="black", pch=21 ) + 
  xlab(paste0("PC1 (", pca$pov[1]*100, "%)" )) + 
  ylab(paste0("PC2 (", pca$pov[2]*100, "%)" )) + 
  theme_minimal() +
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13))

### DE
<<<<<<< HEAD
SS_DE <- function(expr, meta, outcome, nuisance_vars = c("sequencing_month_year", "age", "gender"), incl_cell_props = FALSE) {
=======
SS_DE <- function(expr, meta, outcome, incl_cell_props = FALSE) {
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  
  # Print the outcome
  cat(paste0("Outcome: ", outcome, "\n"))
  
  # Set up meta and expression data
  met <- meta %>% 
<<<<<<< HEAD
    dplyr::select(one_of("sample_identifier",outcome, nuisance_vars)) %>% 
    dplyr::rename(comparison = all_of(outcome)) %>% 
    na.omit()
=======
    dplyr::rename(comparison = all_of(outcome)) %>% 
    filter(!is.na(comparison))
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  
  exp <- expr %>% 
    dplyr::select(one_of("ensembl_gene_id", met$sample_identifier)) %>% 
    column_to_rownames(var = "ensembl_gene_id")
  
<<<<<<< HEAD
  des = paste0("~comparison + ", paste(nuisance_vars, collapse = " + "))
  
=======
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  # How many observations are there 
  no_pats <- met %>% group_by(comparison) %>% summarize(n= n(), .groups = "drop")
  cat(paste0("There are ",  no_pats$n[1], " ", no_pats$comparison[1], 
             " patients and ", 
             no_pats$n[2], " ", no_pats$comparison[2], " patients.\n\n" ))
  
<<<<<<< HEAD
=======
  des = "~comparison + sequencing_month_year + age + gender"
  
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  # Include cell proportions  
  if (incl_cell_props) {
    deconv_pca <- deconv_res$er_tr$cibersort %>% 
      rownames_to_column(var = "sample_identifier") %>% 
      filter(sample_identifier %in% met$sample_identifier) %>% 
      remove_zero_var_cols() %>% 
      column_to_rownames(var = "sample_identifier") %>% 
      perform_pca() 
    
    met <- met %>% left_join(deconv_pca$x,  by = "sample_identifier")
    
    pcs_to_incl <- deconv_pca$pov %>% explain_95(perc = 90)
    cat(paste0("Cell Proportion PCs used: ", pcs_to_incl, "\n"))
    pcs_to_incl <- paste0("PC", 1:pcs_to_incl)
    
    des <- paste0(des, " + ",paste(pcs_to_incl, collapse = " + "))
    
  }
  
<<<<<<< HEAD
  cat("\nDesign formula: ", des, "\n")
  
=======
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  # Check if samples are in the right order
  right_order <-  all(colnames(exp) == met$sample_identifier)
  
  if (!right_order) {
    cat("Things are not in the right order")
    stop()
  }
  
  # Perform differential expression
  de_res <- de(counts = exp, meta = met, FC = 1.5, 
               des = des, 
               main_covar = "comparison", filt = FALSE)
  
  # Perform pathway enrichment
  pthwy <- de_res %>% 
    map(~filter(.x, de == "de")) %>% 
    map(~split(.x, .x$direction)) %>% 
    map(~map(.x, ~pathway_enrichment(.x$entrezgene_id, p_val = 0.05, universe_list = universe$entrezgene_id)))
  
  # Perform mSigDB enrichment
  msigdb <-  de_res %>% 
    map(~filter(.x, de == "de")) %>% 
    map(~split(.x, .x$direction)) %>% 
    map(~map(.x, ~go_enrichment(.x$hgnc_symbol, p_val = 0.05))) 
  
  # Save Results
  res <- list()
  res[["de"]] <- de_res
  res[["pthwy"]] <- pthwy
  res[["msigdb"]] <- msigdb
  return(res)
  }

<<<<<<< HEAD
outcomes <- c("endotype2class", "survive","icu_adm", "culture", "qsofa_High_Low", "High_Low", "HighInt_Low")

=======
outcomes <- c("endotype2class", "survive","icu_adm", "culture", "qsofa_High_Low", "High_Low")
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
res <- outcomes %>% 
  map(~SS_DE(expr, meta, outcome = .x, incl_cell_props = TRUE)) %>% 
  set_names(outcomes)

# Get DE Gene Numbers
res %>% 
  map(~.x$de) %>% 
  map(~map(.x, ~de_gene_numbers(.x)))

### Save results
res %>%  
  write_rds(paste0("./de/results/SS_DE_res.rds"))

