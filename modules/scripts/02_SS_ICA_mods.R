rm(list = ls())
source("./misc/helper_functions.R")
source("./modules/scripts/ModuleFuncs.R")

### Read in expression data 
# Uncorrected
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
<<<<<<< HEAD
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))
=======
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff

meta <- tr_te_dat$er_tr$meta 
expr <- tr_te_dat$er_tr$expr %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  as.matrix() %>% 
  DESeq2::varianceStabilizingTransformation() %>% 
  sva::ComBat(batch = meta$sequencing_month_year) 
<<<<<<< HEAD
deconv_res <- deconv_res$er_tr$cibersort
=======
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
all(colnames(expr) == meta$sample_identifier)

# Corrected 
expr_corr <- read_rds("./modules/results/SS_expr_corr.RDS")
expr_corr <- expr_corr$expr_corr

# Combine 
all(colnames(expr) == colnames(expr_corr))
expr_mats <- list(expr_uncorr = expr, expr_corr = expr_corr)

### ICA Pipeline
# Perform ICA 
ICA_res <- expr_mats %>% 
  map(~perform_ica(.x, kurtosis_keep = 10, qval_filt = 0.001, seed = 1))

<<<<<<< HEAD
# Plot ICA Module Pathways.
ICA_mod_plot <- ICA_res %>% 
  map(~plot_pathwy_mods(.x))

# Get Eigenvectors
ICA_mod_eigen <- ICA_res %>% 
  map2(expr_mats, ~get_mod_eigenvect(.x, .y))

# Get module eigenvector-outcome associations. 
vars <- c("at_ed_body_temperature", "at_ed_heart_rate", "at_ed_systolic", "at_ed_diastolic", "at_ed_altered_mental_state")
ICA_mod_eigen_symp <- list(
  expr_uncorr = lm_clin_eig(ICA_mod_eigen$expr_uncorr, meta,  vars_of_int = vars , deconv = deconv_res, incl_cell_props = TRUE),
  expr_corr = lm_clin_eig(ICA_mod_eigen$expr_corr, meta, vars_of_int = vars , incl_cell_props = FALSE)
)

vars <- c("culture", "qsofa_High_Low", "High_Low","HighInt_Low", "icu_adm", "survive")
ICA_mod_eigen_clin <- list(
  expr_uncorr = lm_clin_eig(ICA_mod_eigen$expr_uncorr, meta,  vars_of_int = vars , deconv = deconv_res, incl_cell_props = TRUE),
  expr_corr = lm_clin_eig(ICA_mod_eigen$expr_corr, meta, vars_of_int = vars , incl_cell_props = FALSE)
)

# Get module signatures 
ICA_mod_sigs <- ICA_mod_eigen_clin %>% 
  map2(ICA_res, ~get_mod_sigs(.x, .y, vars_of_int = vars, top_n_mods = 10, sig_qualifier = "SS") )
=======
# Plot ICA Results
ICA_mod_plot <- ICA_res %>% map(~plot_pathwy_mods(.x))

# Get Eigenvectors
ICA_mod_eigen <- ICA_res %>% map2(expr_mats, ~get_mod_eigenvect(.x, .y))

# Get module eigenvector-outcome associations. 
lm_clin_eig <- function(ICA_eig, meta) {
  
  meta_vars_of_int <- meta %>% 
    dplyr::select(one_of("sample_identifier","survive","survive","icu_adm", "culture", "qsofa_High_Low", "High_Low"))
  ICA_eig <- ICA_mod_eigen$expr_uncorr
  
  # Check things are in the right order 
  right_order <- ICA_eig %>% 
    map(~all(.x$sample_identifier == meta_vars_of_int$sample_identifier)) %>% 
    unlist() %>% 
    all()
  
  if (isFALSE(right_order)) {
    return(c("Things are not in the right order"))
  }
  
  # What variables do we have?
  lm_vars <- colnames(meta_vars_of_int[,-1])
  
  # Initiate loop - loop through variables, and then module
  res <- list()
  for (var in lm_vars) {
    for (mod in names(ICA_eig)) {
      # Create a new df with the var and eigenvector the module
      lm_df <- ICA_eig[[mod]] %>% 
        left_join(meta_vars_of_int, by = "sample_identifier") %>% 
        dplyr::select(one_of("sample_identifier","PC1", var)) %>% 
        na.omit()
      
      # GLM - PC1 ~ var
      res[[var]][[mod]] <- glm(formula = formula(paste0("PC1 ~ ", var)), data = lm_df ) %>% 
          broom::tidy()
    }
  }
  
  # Format res
  res_df <- res %>% 
    map(~bind_rows(.x, .id = "comp")) %>% 
    bind_rows(.id = "variable") %>% 
    filter(term != "(Intercept)" ) %>% 
    #mutate(adj_p_value = p.adjust(p.value)) %>% 
    mutate(log10_p.value = -log10(p.value))
  
  # Plot res
  paletteLength <- 50
  myColor <- colorRampPalette(c("white", "whit", "blue"))(paletteLength)
  res_df %>% 
    dplyr::select(one_of("variable", "comp", "log10_p.value")) %>% 
    pivot_wider(names_from = "variable", values_from = "log10_p.value") %>% 
    column_to_rownames(var = "comp") %>% 
    pheatmap(color = colorRampPalette(c( "white", "blue"))(50))
    # ggplot(aes(x = variable, y = comp, fill = log10_p.value)) +
    # geom_tile() + 
    # scale_fill_viridis(option = "magma", direction = -1) +
    # scale_x_discrete(drop=FALSE) +
    # theme_minimal() + ylab("") + xlab("") +
    # theme( axis.text.x = element_text(angle = 90),  strip.text.y = element_text(angle = 0, color = "black" )) 
  hey %>% View()
  
  return(res)
} 


clinical_assoc <- ICA_mod_eigen %>% map(~lm_clin_eig(.x, meta))
clinical_assoc$expr_uncorr$survive

>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff

# Save
ICA_res %>% write_rds("./modules/results/SS_ICA_res.rds")
ICA_mod_plot %>% write_rds("./modules/results/SS_ICA_mod_plot.rds")
ICA_mod_eigen %>% write_rds("./modules/results/SS_ICA_mod_eigen.rds")
<<<<<<< HEAD
ICA_mod_eigen_clin %>% write_rds("./modules/results/SS_ICA_mod_eigen_clin.rds")
ICA_mod_eigen_symp %>% write_rds("./modules/results/SS_ICA_mod_eigen_symp.rds")
ICA_mod_sigs %>% write_rds("./modules/results/SS_ICA_mod_sigs.rds")
=======
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff

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
