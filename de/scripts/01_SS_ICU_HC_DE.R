rm(list = ls())
source("./misc/helper_functions.R")

### Read in data
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

### Set up expr
expr <- tr_te_dat$er_tr$expr %>% 
  left_join(tr_te_dat$icu$expr) %>% 
  left_join(tr_te_dat$hc$expr) %>% 
  column_to_rownames(var = "ensembl_gene_id")

meta <- tr_te_dat %>% 
  map(~.x$meta) %>% 
  keep(!names(.) %in% c("er_te")) %>% 
  map(~dplyr::select(.x, one_of("sample_identifier", "condition", "sequencing_month_year", "age", "gender"))) %>% 
  bind_rows()
all(colnames(expr) == meta$sample_identifier)

deconv_res <- deconv_res$er_tr$cibersort %>% 
  bind_rows(deconv_res$icu$cibersort) %>% 
  bind_rows(deconv_res$hc$cibersort) %>% 
  rownames_to_column(var = "sample_identifier")
all(meta$sample_identifier == deconv_res$sample_identifier)

### Plot cell proportions 
keep <- c("monocyte","neutrophil", "macrophage_m0","b_cell_plasma", "t_cell_cd4_naive", "t_cell_regulatory_tregs", "t_cell_cd8")
deconv_res_df <- deconv_res %>% 
  pivot_longer(cols = -c("sample_identifier"), "cell_type") %>% 
  left_join(meta) %>% 
  filter(cell_type %in% keep)
  
deconv_res_df %>% 
  mutate(condition = factor(condition, levels = c("healthy_control", "suspected_sepsis", "icu"), labels = c("HC", "SS", "ICU"))) %>% 
  ggplot(aes(x = condition, y = value, fill = condition)) + 
  geom_boxplot() + 
  facet_wrap(~cell_type, nrow = 1, scale = "free") +
  xlab("") + ylab("")
ggsave("./de/figures/ss_icu_hc_cell_prop.png", width = 14, height = 3, scale = 0.9)

### Perform DE/Pathway Enrichment
SS_ICU_HC_DE <- function(comparison = c("healthy_control", "suspected_sepsis"),  incl_cell_props = FALSE ) {
  
  met <- meta %>% 
    filter(condition %in% comparison) %>% 
    mutate(condition = factor(condition, levels = comparison))
  
  exp <- expr %>% 
    dplyr::select(one_of(met$sample_identifier))
  
  # How many observations are there 
  no_pats <- met %>% group_by(condition) %>% summarize(n= n(), .groups = "drop")
  cat(paste0("There are ",  no_pats$n[1], " ", no_pats$condition[1], 
             " patients and ", 
             no_pats$n[2], " ", no_pats$condition[2], " patients.\n\n" ))
  
  des = "~condition + age + gender"
  
  # Include cell proportions  
  if (incl_cell_props) {
  deconv_pca <- deconv_res %>% 
    filter(sample_identifier %in% met$sample_identifier) %>% 
    remove_zero_var_cols() %>% 
    column_to_rownames(var = "sample_identifier") %>% 
    perform_pca() 
  
  met <- met %>% left_join(deconv_pca$x, by = "sample_identifier")
  
  pcs_to_incl <- deconv_pca$pov %>% explain_95(perc = 90)
  pcs_to_incl <- paste0("PC", 1:pcs_to_incl)
  
  des <- paste0(des, " + ",paste(pcs_to_incl, collapse = " + "))
  
  }
  
  if (length(unique(met$sequencing_month_year)) > 1 ) {
    paste0(des, " + sequencing_month_year")
  }
  
  # Check if things are in the right order
  right_order <- all(colnames(exp) == met$sample_identifier)
  if (!right_order) {
    cat("Wrong order.")
    stop()
  }

  ### DE 
  de_res <- de(exp, met, des = des, main_covar = "condition", filt = FALSE)
  
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

res <- list(
  ss_vs_hc = SS_ICU_HC_DE(c("healthy_control", "suspected_sepsis"), des = "~condition + sequencing_month_year + age + gender", incl_cell_props = T),
  icu_vs_hc = SS_ICU_HC_DE(c("healthy_control", "icu"),  des = "~condition + sequencing_month_year + age + gender", incl_cell_props = T),
  ss_vs_icu = SS_ICU_HC_DE(c("suspected_sepsis", "icu"),  des = "~condition + age + gender", incl_cell_props = T)
)

# Get DE Gene Numbers
res %>% 
  map(~.x$de) %>% 
  map(~map(.x, ~de_gene_numbers(.x)))

### Plot pathways 
pthwy_htmaps <- function(de_list, pthwy_group ){
  
  # Set up DF
  pthwy_plot <- de_list %>% 
    map(~.x$pthwy) %>% 
    map(~map(.x, ~bind_rows(.x, .id = "direction" ))) %>% 
    map(~bind_rows(.x, .id ="comparison"  )) %>% 
    bind_rows() %>% 
    separate(BgRatio,into = c("M", "N")) %>%
    mutate(Ratio = Count/as.numeric(M) ) %>% 
    left_join(pathway_hier,  by = c("ID" = "enr_pathway")) 
  
  comparison_names = c("SS vs\n HC", "ICU vs\n HC", "ICU vs\n SS")
  names(comparison_names) = unique(pthwy_plot$comparison)
  
  # Plot
  plt <- pthwy_plot %>% 
    filter(top_level_pathway_descrip == pthwy_group) %>% 
    mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down"))) %>% 
    #mutate(comparison = factor(comparison, levels =  names(comparison_names), labels =  comparison_names)) %>% 
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

plts <- list()
keep = c("Immune System", "Programmed Cell Death", "Metabolism", "Hemostasis", "Cellular responses to external stimuli",
         "Signal Transduction", "Transport of small molecules")
for (i in keep) {
  plts[[i]] <- res %>% 
    pthwy_htmaps(pthwy_group = i)
}

plts$`Immune System`
plts$Hemostasis
# ggsave("./de/figures/ss_icu_hc_pthwy_enr_Immune_System.png", plts$`Immune System`, width = 9, height = 5.5 )
# ggsave("./de/figures/ss_icu_hc_pthwy_enr_Programmed_Cell_Death.png", plts$`Programmed Cell Death`, width = 6, height = 2.5)
# ggsave("./de/figures/ss_icu_hc_pthwy_enr_Metabolism.png", plts$Metabolism, width = 6.5, height = 6)
# ggsave("./de/figures/ss_icu_hc_pthwy_enr_Hemostasis.png", plts$Hemostasis, width = 6.5, height = 3.5)
# ggsave("./de/figures/ss_icu_hc_pthwy_enr_Cellular_Response.png", plts$`Cellular responses to external stimuli`, width = 7, height = 5)
# ggsave("./de/figures/ss_icu_hc_pthwy_enr_Signal_Transduction.png", plts$`Signal Transduction`, width = 7, height = 12)
