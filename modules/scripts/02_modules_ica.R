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
  as.data.frame()

### Read in Cell proportion corrected data 
expr_cor <- read_rds("./modules/results/expr_corr_resid.rds")

all(colnames(expr_vst) == colnames(expr_cor))

### Combine 
expr_mods <- list(expr_vst = expr_vst, expr_cor = expr_cor)

### PCA 
plot_PCA <- function(expr_df, pca_title, outcome) {
  pca <- expr_df %>% t() %>% as.data.frame() %>% remove_zero_var_cols() %>% perform_pca()
  pca_fig <- pca$x %>% 
    left_join(meta) %>% 
    ggplot(aes(x = PC1, y = PC2, color = meta[[outcome]])) + geom_point() + ggtitle(pca_title)
  return(pca_fig)
}

expr_mods_pca <- expr_mods %>% 
  map2(names(.), ~plot_PCA(.x, outcome = "outcome_icu_admission", pca_title = .y))

#### Determine the number of components to learn   
no_comps <- expr_mods %>% 
  map(~t(.x) %>% as.data.frame() %>% remove_zero_var_cols() %>% perform_pca() %>% pluck("pov") %>% explain_95())

# Now perform ICA
perform_ica <- function(expr, comp_no = ncol(expr), kurtosis_keep = 10, qval_filt = 0.001, seed = 1) {
  # Perform ICA
  set.seed(seed)
  ICA <- fastICA::fastICA(expr, n.comp = comp_no)
  
  # Get the S matrix
  ICA_S <- ICA$S %>% as.data.frame() %>% set_names(paste0("Mod_", 1:ncol(.)))
  
  # Keep certain components based on kurtosis 
  comp_keep <- ICA_S %>% map_dbl(~e1071::kurtosis(.x)) %>%  keep(~.x >= kurtosis_keep) %>% names()
  ICA_S_filt <- ICA_S[,comp_keep]
  
  # Filter components of interest with FDR estimation
  ICA_mods <- ICA_S_filt %>% map(~fdrtool::fdrtool(.x, plot = FALSE, verbose = FALSE))
  
  # Filter genes in each componennt 
  ICA_mods_filt <- list()
  for (comp in names(ICA_mods)) {
    mod_df <- data.frame(qval = ICA_mods[[comp]][["qval"]], ensembl_gene_id = rownames(ICA_S))
    mod_df <- mod_df %>% dplyr::filter(qval <= qval_filt) %>% left_join(universe, by = "ensembl_gene_id")
    ICA_mods_filt[[comp]] <- mod_df
  }
  
  return(list(ICA_S_filt = ICA_S_filt, ICA_mods_filt = ICA_mods_filt))
}

ICA_res <- expr_mods %>% 
  map2(no_comps, ~perform_ica(.x, comp_no = .y, kurtosis_keep = 10, qval_filt = 0.001, seed = 1))

ICA_res %>% 
  map(~map(.x$ICA_mods_filt, ~nrow(.x))) %>% 
  map(~flatten_int(.x)) %>% 
  map(~mean(.x))

ICA_jac_clust <- ICA_res %>% 
  map(~get_jac_mat(.x$ICA_mods_filt) %>% 
        t() %>% 
        vegan::vegdist(method = "jaccard", binary = TRUE, diag = TRUE) %>% 
        hclust %>% plot())

#### Plot pathway enrichment of Modules
plot_pathwy_mods <- function(ica_obj) {
  
  # Pathway Enrichment
  ICA_mods_filt_pathway_enr <- ica_obj$ICA_mods_filt %>% 
    map(~.x$entrezgene_id) %>% 
    map(~pathway_enrichment(.x, universe$entrezgene_id, ID = "entrez", p_val = 0.01))
  
  ICA_mods_filt_pathway_enr_df <- ICA_mods_filt_pathway_enr %>% 
    bind_rows(.id = "comp")
  
  pthwy_enr <- ICA_mods_filt_pathway_enr_df %>% 
    mutate(comp = paste0("Mod ", str_split(comp, "_", simplify = TRUE)[,2])) %>% 
    mutate(comp = factor(comp)) %>% # Fix this, it still excludes pathways
    separate(BgRatio,into = c("M", "N")) %>%
    mutate(Ratio = Count/as.numeric(M) ) %>% 
    left_join(pathway_hier, by = c("ID" = "enr_pathway")) %>% 
    group_by(comp) %>% 
    top_n(2, Ratio) %>% 
    ungroup() %>% 
    mutate(Description = factor(Description, levels = rev(unique(.$Description)))) %>% 
    ggplot(aes(x = comp, y = Description, fill = Ratio, size = -log10(p.adjust) )) + 
    geom_point(color = "black", pch = 22) +
    # facet_grid(rows = vars(one_lower_level_pathway_descrip_clean_name), scales = "free_y", space = "free_y") + 
    scale_fill_viridis(option = "magma", direction = -1) +
    scale_x_discrete(drop=FALSE) +
    theme_minimal() + ylab("") + xlab("") +
    theme( axis.text.x = element_text(angle = 90),  strip.text.y = element_text(angle = 0, color = "black" )) 
  
  # Perform mSigDB enrichment
  ICA_mods_filt_msigdb_enr <-  ica_obj$ICA_mods_filt %>% 
    map(~pull(.x, hgnc_symbol)) %>%
    map(~go_enrichment(.x), p_val = 0.01)
  
  ICA_mods_filt_msigdb_enr_df <- ICA_mods_filt_msigdb_enr %>% 
    map(~keep(.x, names(.) %in% c("MSigDB_Hallmark_2020"))) %>% 
    map(~bind_rows(.x, .id = "database")) %>% 
    bind_rows(.id = "comp")
  
  msigdb_enr <- ICA_mods_filt_msigdb_enr_df %>% 
    mutate(comp = paste0("Mod ", str_split(comp, "_", simplify = TRUE)[,2])) %>% 
    mutate(comp = factor(comp)) %>% 
    separate(Overlap,into = c("M", "N")) %>%
    mutate(Ratio = as.numeric(M)/as.numeric(N) ) %>% 
    group_by(comp) %>% 
    top_n(2, Ratio) %>% 
    ungroup() %>% 
    mutate(Term = factor(Term, levels = rev(unique(.$Term)))) %>% 
    ggplot(aes(x = comp, y = Term, fill = Ratio, size = -log10(Adjusted.P.value) )) + 
    geom_point(color = "black", pch = 22) +
    scale_fill_viridis(option = "magma", direction = -1) +
    scale_x_discrete(drop=FALSE) +
    theme_minimal() + ylab("") + xlab("") +
    theme( axis.text.x = element_text(angle = 90),  strip.text.y = element_text(angle = 0, color = "black" )) 
  
  return(list(reactome = pthwy_enr, msigdb = msigdb_enr ))
  }

ICA_mod_plot <- ICA_res %>% map(~plot_pathwy_mods(.x))

# Save
ICA_mod_plot %>% write_rds("./modules/results/ICA_mod_plot.rds")
ICA_res %>% write_rds("./modules/results/ICA_res.rds")
