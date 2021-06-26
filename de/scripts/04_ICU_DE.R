rm(list = ls())
source("./misc/helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

expr <- tr_te_dat$icu$expr
meta <- tr_te_dat$icu$meta 
all(colnames(expr)[-1] == meta$sample_identifier)

<<<<<<< HEAD
# Add in new outcome only relevant to ICU patients - High Severity patients who died vs did not. 
meta %<>% 
  mutate(
    High_Survive = case_when(
      sofa_sev_24 == "high" & survive == "dead" ~ "high_dead", 
      sofa_sev_24 == "high" & survive == "survive" ~ "high_survive", 
      TRUE ~ NA_character_),
    High_Survive = factor(High_Survive, levels = c("high_survive", "high_dead")) 
    )

# Increase High severity cut off
meta %<>% 
  mutate(
    High_Low_10_Cutoff = case_when(
      sofa_24 >= 10 ~ "high", 
      sofa_24 < 2 ~ "low", 
      TRUE ~ NA_character_),
    High_Low_10_Cutoff = factor(High_Low_10_Cutoff, levels = c("low", "high")) 
    )

=======
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
### PCA
pca <- expr %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  as.matrix() %>% 
  DESeq2::varianceStabilizingTransformation()  %>% 
  #sva::ComBat(batch = meta$sequencing_month_year) %>%
  t() %>% 
  perform_pca()

pca$x %>%
  left_join(meta) %>% 
<<<<<<< HEAD
  ggplot(aes(x = PC1, y = PC2,  fill = survive )) +
=======
  ggplot(aes(x = PC1, y = PC2,  fill =High_Low)) +
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  geom_point(size = 4, colour="black", pch=21 ) + 
  xlab(paste0("PC1 (", pca$pov[1]*100, "%)" )) + 
  ylab(paste0("PC2 (", pca$pov[2]*100, "%)" )) + 
  theme_minimal() +
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13))

### DE
<<<<<<< HEAD
ICU_DE <- function(expr, meta, outcome, nuisance_vars = c("age", "gender"),  incl_cell_props = FALSE) {
=======
ICU_DE <- function(expr, meta, outcome, incl_cell_props = FALSE) {
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

    # How many observations are there 
=======
  des = "~comparison + age + gender"
  
  # How many observations are there 
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  no_pats <- met %>% group_by(comparison) %>% summarize(n= n(), .groups = "drop")
  cat(paste0("There are ",  no_pats$n[1], " ", no_pats$comparison[1], 
             " patients and ", 
             no_pats$n[2], " ", no_pats$comparison[2], " patients.\n\n" ))
  
  # Include cell proportions  
  if (incl_cell_props) {
    deconv_pca <- deconv_res$icu$cibersort %>% 
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
  cat("Design formula: ", des, "\n")
  
  # Check if samples are in the right order
  right_order <-  all(colnames(exp) == met$sample_identifier)
=======
  # Check if samples are in the right order
  right_order <-  all(colnames(exp) == met$sample_identifier)
  
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
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
outcomes <- c("High_Low", "survive", "culture", "High_Low_10_CutOff")
res <- outcomes %>% 
  map(~ICU_DE(expr, meta, outcome = .x, nuisance_vars = c("age", "gender"), incl_cell_props = TRUE)) %>% 
  set_names(outcomes)

#### New outcomes or models to build
res[["survive_sev_corr"]] <- ICU_DE(expr, meta, outcome = "survive", nuisance_vars = c("sofa_24", "age", "gender"), incl_cell_props = FALSE)
res[["High_Survive"]] <- ICU_DE(expr, meta, outcome = "High_Survive", nuisance_vars = c("age", "gender"), incl_cell_props = FALSE)
=======
outcomes <- c("High_IntLow","High_Low", "survive")
res <- outcomes %>% 
  map(~ICU_DE(expr, meta, outcome = .x, incl_cell_props = TRUE)) %>% 
  set_names(outcomes)

### Save results
res %>%  
  write_rds(paste0("./de/results/ICU_DE_res.rds"))
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff

# Get DE Gene Numbers
res %>% 
  map(~.x$de) %>% 
  map(~map(.x, ~de_gene_numbers(.x)))

<<<<<<< HEAD
res$Worse_Same_Better <- NULL
res$High_Low <- NULL

### Save results
res %>%  
  write_rds(paste0("./de/results/ICU_DE_res.rds"))
=======
### Plot pathways 
pthwy_plot <- res %>% 
  map(~.x$pthwy) %>% 
  map(~map(.x, ~bind_rows(.x, .id = "direction"))) %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  bind_rows() %>% 
  separate(BgRatio,into = c("M", "N")) %>%
  mutate(Ratio = Count/as.numeric(M) ) %>% 
  left_join(pathway_hier,  by = c("ID" = "enr_pathway")) 

pthwy_plot %>% 
  #filter(top_level_pathway_descrip == pthwy_group) %>% 
  mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down"))) %>% 
  # mutate(comparison = factor(comparison, levels =  names(comparison_names), labels =  comparison_names)) %>% 
  mutate(one_lower_level_pathway_descrip_clean_name = str_wrap(one_lower_level_pathway_descrip_clean_name, width = 20)) %>% 
  #mutate(Description = str_wrap(Description, width = 40)) %>% 
  ggplot(aes(x = direction, y = Description, fill = Ratio)) + 
  geom_tile() + 
  facet_grid(cols = vars(comparison), rows = vars(top_level_pathway_descrip, one_lower_level_pathway_descrip_clean_name), scale = "free_y", space = "free_y") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA), plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(angle = 0)) + 
  scale_fill_viridis(direction = -1) + 
  ylab("") + xlab("")




pthwy_htmaps <- function(de_list, pthwy_group ){
  
  # Set up DF
  pthwy_plot <- res %>% 
    map(~.x$pthwy) %>% 
    map(~map(.x, ~bind_rows(.x, .id = "direction"))) %>% 
    map(~bind_rows(.x, .id = "comparison")) %>% 
    bind_rows() %>% 
    separate(BgRatio,into = c("M", "N")) %>%
    mutate(Ratio = Count/as.numeric(M) ) %>% 
    left_join(pathway_hier,  by = c("ID" = "enr_pathway")) 
  
  #comparison_names = c("SS vs\n HC", "ICU vs\n HC", "ICU vs\n SS")
  #names(comparison_names) = unique(pthwy_plot$comparison)
  
  # Plot
  plt <- pthwy_plot %>% 
    filter(top_level_pathway_descrip == pthwy_group) %>% 
    mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down"))) %>% 
   # mutate(comparison = factor(comparison, levels =  names(comparison_names), labels =  comparison_names)) %>% 
    mutate(one_lower_level_pathway_descrip_clean_name = str_wrap(one_lower_level_pathway_descrip_clean_name, width = 20)) %>% 
    mutate(Description = str_wrap(Description, width = 40)) %>% 
    ggplot(aes(x = direction, y = Description, fill = Ratio)) + 
    geom_tile() + 
    facet_grid(cols = vars(comparison), rows = vars(one_lower_level_pathway_descrip_clean_name), scale = "free_y", space = "free_y") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA), plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(angle = 0)) + 
    scale_fill_viridis(direction = -1) + 
    ylab("") + xlab("") + 
    ggtitle(pthwy_group)
  return(plt)
}

res %>% 
  pthwy_htmaps(pthwy_group = "Immune System")
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
