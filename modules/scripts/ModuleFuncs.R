#### Perform ICA
perform_ica <- function(expr, kurtosis_keep = 10, qval_filt = 0.001, seed = 1) {
  
  # Determine the number of components to learn   
  comp_no <- expr %>% 
    t() %>% as.data.frame() %>% remove_zero_var_cols() %>% perform_pca() %>% pluck("pov") %>% explain_95(perc = 95)
  cat(paste0("Performing ICA with ", comp_no, " components.\n"))
  
  # Perform ICA
  set.seed(seed)
  ICA <- fastICA::fastICA(expr, n.comp = comp_no)
  
  # Get the S matrix
  ICA_S <- ICA$S %>% as.data.frame() %>% 
    set_names(paste0("Mod_",   str_pad(1:ncol(.), width = 3, side = "left", pad = "0")    ))
  
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

#### Plot pathway enrichment of Modules
plot_pathwy_mods <- function(ica_obj) {
  
  # Pathway Enrichment
  ICA_mods_filt_pathway_enr <- ica_obj$ICA_mods_filt %>% 
    map(~.x$entrezgene_id) %>% 
    map(~pathway_enrichment(.x, universe$entrezgene_id, ID = "entrez", p_val = 0.01))
  
  ICA_mods_filt_pathway_enr_df <- ICA_mods_filt_pathway_enr %>% 
    bind_rows(.id = "comp") %>% 
    separate(BgRatio,into = c("M", "N")) %>%
    mutate(Ratio = Count/as.numeric(M) ) %>% 
    left_join(pathway_hier, by = c("ID" = "enr_pathway")) %>% 
    group_by(comp) %>% 
    top_n(1, Ratio) %>% 
    ungroup()
  
  ICA_mods_filt_pathway_enr_df_plot <- ICA_mods_filt_pathway_enr_df %>% 
    #mutate(comp = factor(comp, levels = names(ICA_mods_filt_pathway_enr))) %>% 
    ggplot(aes(y = Description, x = comp, fill = Ratio )) + 
    geom_tile() + 
    facet_grid(rows = vars(top_level_pathway_descrip), scales = "free_y", space = "free_y") + 
    scale_fill_viridis("magma", direction = -1) +
    scale_x_discrete(drop=FALSE) +
    theme_minimal() + ylab("") + xlab("") +
    theme( axis.text.x = element_text(angle = 90),  strip.text.y = element_text(angle = 0, color = "black" )) 
  
  # Perform mSigDB enrichment
  ICA_mods_filt_msigdb_enr <-  ica_obj$ICA_mods_filt %>%
    map(~pull(.x, hgnc_symbol)) %>%
    map(~go_enrichment(.x,  p_val = 0.01))

  ICA_mods_filt_msigdb_enr_df <- ICA_mods_filt_msigdb_enr %>%
    map(~keep(.x, names(.) %in% c("MSigDB_Hallmark_2020"))) %>%
    map(~bind_rows(.x, .id = "database")) %>%
    bind_rows(.id = "comp") %>%
    separate(Overlap,into = c("M", "N")) %>%
    mutate(Ratio = as.numeric(M)/as.numeric(N) ) %>%
    group_by(comp) %>%
    top_n(1, Ratio) %>%
    ungroup()

  ICA_mods_filt_msigdb_enr_df_plot <- ICA_mods_filt_msigdb_enr_df %>%
    mutate(Term = factor(Term, levels = rev(unique(.$Term)))) %>%
    ggplot(aes(x = comp, y = Term, fill = Ratio )) +
    geom_tile() +
    scale_fill_viridis(option = "magma", direction = -1) +
    scale_x_discrete(drop=FALSE) +
    theme_minimal() + ylab("") + xlab("") +
    theme( axis.text.x = element_text(angle = 90),  strip.text.y = element_text(angle = 0, color = "black" ))
  
  return(list(reactome = ICA_mods_filt_pathway_enr_df_plot, msigdb = ICA_mods_filt_msigdb_enr_df_plot ))

}

get_mod_eigenvect <- function(ICA_mods, expr_mat) {
  
  eig <- ICA_mods$ICA_mods_filt %>% 
    discard(~nrow(.x) == 0) %>% # For some reason there are modules with no genes in them... need to check this out. 
    map(~as.matrix(expr_mat)[.x$ensembl_gene_id, ]) %>% 
    map(~t(.x) %>% perform_pca() %>% pluck("x")) %>% 
    map(~dplyr::select(.x, one_of("sample_identifier", "PC1")))
  
  return(eig)
  
}

#### Associate modules to clinical vars_of_int
lm_clin_eig <- function(ICA_eig, meta, vars_of_int, deconv = NULL, incl_cell_props = FALSE) {
  #vars_of_int = c("survive","survive","icu_adm", "culture", "qsofa_High_Low", "High_Low")
  met <- meta
  
  # Check things are in the right order 
  right_order <- ICA_eig %>% 
    map(~all(.x$sample_identifier == met$sample_identifier)) %>% 
    unlist() %>% 
    all()
  
  if (isFALSE(right_order)) {
    return(c("Things are not in the right order"))
  }
  
  # Design formula 
  des = "PC1 ~ "
  
  # Include cell proportions  
  if (incl_cell_props) { 
    
    if (is.null(deconv)) {
      stop()
    }
    
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
    
    if (length(unique(met$sequencing_month_year)) > 1 ) {
      paste0(des, " sequencing_month_year + ")
    }
    
    des <- paste0(des, "gender + age + ", paste(paste0("CiSo_",pcs_to_incl), collapse = " + "), " + ")
  }
  
  # Initiate loop - loop through variables, and then module
  res <- list()
  for (var in  vars_of_int) {
    for (mod in names(ICA_eig)) {
      # Create a new df with the var and eigenvector the module
      lm_df <- ICA_eig[[mod]] %>% 
        left_join(met, by = "sample_identifier")
      
      # GLM - PC1 ~ var
      res[[var]][[mod]] <- glm(formula = formula(paste0(des,  var)), data = lm_df ) %>% 
        broom::tidy()
    }
  }
  
  # Format res and plot 
  res_df <- res %>%
    map(~bind_rows(.x, .id = "comp")) %>%
    bind_rows(.id = "variable") %>%
    filter(term != "(Intercept)" ) %>%
    mutate(adj_p_value = p.adjust(p.value)) %>% 
    filter(str_detect(term,  paste(vars_of_int, collapse = "|") )) %>% 
    mutate(log10_p.value = -1*log10(p.value)) 
  
  # Plot res
  # paletteLength <- 50
  # myColor <- colorRampPalette(c("white", "white", "blue"))(paletteLength)
  # res_df %>% 
  #   dplyr::select(one_of("variable", "comp", "log10_p.value")) %>%
  #   pivot_wider(names_from = "variable", values_from = "log10_p.value") %>%
  #   column_to_rownames(var = "comp") %>% 
  #   pheatmap(color = colorRampPalette(c( "white", "blue"))(50))
  # 
  res_df_plt <- res_df %>% 
    mutate(variable = factor(variable, levels = c(vars_of_int))) %>% 
    ggplot(aes(x = variable, y = comp, fill = log10_p.value)) +
    geom_tile() +
    scale_fill_viridis(option = "cividis", direction = -1) + 
    theme_minimal() + ylab("") + xlab("") +
    theme( axis.text.x = element_text(angle = 90),  strip.text.y = element_text(angle = 0, color = "black" ))
  
  return(list(res = res, res_df_plt = res_df_plt))
} 

get_mod_sigs <- function(ICA_eig_clin, ICA_mods, vars_of_int, top_n_mods = 5,  sig_qualifier = NULL) {
  # Get the top 5 signatures (wrt p value) for each variable/outcome
  mod_sig_top <- ICA_eig_clin$res %>% 
    map(~bind_rows(.x, .id = "comp")) %>%
    bind_rows(.id = "variable") %>%
    filter(term != "(Intercept)" ) %>%
    mutate(adj_p_value = p.adjust(p.value)) %>% 
    filter(str_detect(term,  paste(vars_of_int, collapse = "|") )) %>% 
    mutate(log10_p.value = -1*log10(p.value)) %>% 
    group_by(variable) %>% 
    top_n(top_n_mods, log10_p.value) %>% 
    ungroup()
  
  # Create Signature names - Add sig qualifier if provided 
  if (is.null(sig_qualifier)) {
    mod_sig_top %<>% mutate(sig_name = paste0(variable, "_", comp )) 
    mod_sig_top %<>% mutate(sig_name_no_clin = comp) 
    
  } else {
    mod_sig_top %<>% mutate(sig_name = paste0(sig_qualifier, "_", variable, "_", comp )) 
    mod_sig_top %<>% mutate(sig_name_no_clin = paste0(sig_qualifier, "_", comp)) 
  }
  # Get module genes and add to the DF
  mod_genes <- list()
  for (i in 1:nrow(mod_sig_top)) {
    comp <- mod_sig_top$comp[i]
    mod_genes[[i]] <- ICA_mods$ICA_mods_filt[[comp]]$ensembl_gene_id
  }
  mod_sig_top[["mod_genes"]] <- mod_genes

  # Return 
  mod_sig_top %<>% dplyr::select(one_of("variable","comp", "sig_name", "sig_name_no_clin", "mod_genes"))
  return(mod_sig_top)
}
