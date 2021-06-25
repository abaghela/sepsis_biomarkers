#### Script to measure enrichment of ICU signatures in SS Train and Test patients.
rm(list = ls())
source("./misc/helper_functions.R")

### Read in signatures 
all_sigs <- read_rds("./train_test_mods/results/DE_ICA_sigs.RDS")
all_sigs <- all_sigs[str_detect(names(all_sigs),"ICU")] # Only interested in the ICU signatures

SS_ROAST_Enrich <- function(SS_cohort , outcomes, sigs) {
  
  ### Set up expression and meta data depending on SS TR or TE selected. 
  if (SS_cohort %in% c("er_tr", "er_te")) {
    expr <- tr_te_dat[[SS_cohort]]$expr
    meta <- tr_te_dat[[SS_cohort]]$meta 
  } else if (SS_cohort == "all") {
    
    expr <- tr_te_dat[["er_tr"]]$expr %>% 
      full_join(tr_te_dat[["er_te"]]$expr, by = "ensembl_gene_id")
    meta <- tr_te_dat[["er_tr"]]$meta %>% 
      bind_rows(tr_te_dat[["er_te"]]$meta)
    
  } else {
    return("Must select either 'er_tr' or 'er_te' or 'all'. ")
  }
  
  res <- list()
  for (i in outcomes) {
    met <- meta
    
    if (is.null(met[[i]])){
      cat(paste0(i, " does not exist in meta, moving on... \n"))
      next
    }
    
    # Print the outcome
    cat(paste0("Outcome: ", i, "\n"))
    
    # Set up meta and expression data
    met <- met %>% 
      dplyr::rename(comparison = all_of(i)) %>% 
      filter(!is.na(comparison))
    
    exp <- expr %>% 
      dplyr::select(one_of("ensembl_gene_id", met$sample_identifier)) %>% 
      column_to_rownames(var = "ensembl_gene_id") %>% 
      as.matrix() %>% 
      varianceStabilizingTransformation()
    
    # Check if samples are in the right order
    right_order <-  all(colnames(exp) == met$sample_identifier)
    
    if (!right_order) {
      cat("Things are not in the right order")
      next
    }
    
    # How many observations are there?
    no_pats <- met %>% group_by(comparison) %>% summarize(n= n(), .groups = "drop")
    cat(paste0("There are ",  no_pats$n[1], " ", no_pats$comparison[1], 
               " patients and ", 
               no_pats$n[2], " ", no_pats$comparison[2], " patients.\n\n" ))
    
    # ROAST
    #sigs_for_outcome <- sigs[str_detect(names(sigs), i)]
    sigs_for_outcome <- sigs
    
    mod.mat <- model.matrix(~comparison, data = met)
    res[[i]] <- limma::mroast(exp, sigs_for_outcome, design = mod.mat, nrot = 10000)
    
  }
  return(res)
}

### Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

# Outcomes of interest.
outcomes <- c("endotype2class", "survive","icu_adm", "culture", "qsofa_High_Low", "High_Low", "HighInt_Low", "High_IntLow")

# cohort_sets <- c("er_tr", "er_te")
# res <- cohort_sets %>% 
#   map(~SS_ROAST_Enrich(.x, outcomes, all_sigs)) %>% 
#   set_names(cohort_sets)
res <- SS_ROAST_Enrich("all", outcomes, all_sigs)

res %>% write_rds("./train_test_mods/results/SS_ROAST_Enrich.rds")

#### PLOT
SS_ROAST_PLOT <- function(roast_res, gene_sets) {
  
  outcomes <- c("survive", "icu_adm", "culture", "High_Low", "HighInt_Low")
  
  # Get results into a data frame
  roast_enr_df <- roast_res %>% 
    keep(names(.) %in% outcomes) %>% 
    map(~rownames_to_column(.x, var = "signature")) %>% 
    bind_rows(.id = "comparison") 
  
  # PICK GENE SETS
  # roast_enr_df %>% 
  #   group_by(comparison) %>% 
  #   top_n(-10, PValue) %>% 
  #   View()
  
  roast_enr_plt <- roast_enr_df %>% 
    dplyr::select(one_of("comparison", "signature", "NGenes", "PropDown", "PropUp", "PValue")) %>% 
    filter(signature %in% gene_sets) %>% 
    pivot_longer(cols = c("PropUp", "PropDown"), "direction") %>% 
    mutate(comparison = factor(comparison, levels = outcomes)) %>% 
    ggplot(aes(x = comparison, y = value, fill = direction)) + 
    geom_bar(position="stack", stat="identity") + 
    facet_grid(cols = vars(signature)) +
    theme(axis.text.x = element_text(angle = 90)) + 
    ylab("Proportion of Gene Set Up/Down Regulated") + xlab("") +
    scale_fill_manual("",values = alpha(c("red", "darkgreen"), 0.7))
  
  return(roast_enr_plt)
  
}

res_plt <- SS_ROAST_PLOT(res, c("ICU_Mod_040", "ICU_Mod_051", "ICU_Mod_021", "ICU_survive_DE_UP"))
ggsave("./train_test_mods/figures/ICU_SIGS_ENRICHED_IN_ER.png", res_plt, scale = 0.9, height = 5, width = 9)
