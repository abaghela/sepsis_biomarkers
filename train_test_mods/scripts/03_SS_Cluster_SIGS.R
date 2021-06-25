# Build Elastic Nets for mortality/severity . 

rm(list = ls())
source("./misc/helper_functions.R")

### Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))
all_sigs <- read_rds("./train_test_mods/results/DE_ICA_sigs.RDS")

gsva_enrichment <- function(SS_cohort) {
  
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
  
  ### Read in signatures - currently interested in the reduced one
  #reduced_elasnet_sigs <- read_rds("./train_test_mods/results/ICU_Train_Mods.RDS")
  # SIG <- reduced_elasnet_sigs$survive$ICU_combine_DE_UP %>%
  #   extract_coefs() %>%
  #   pull(ensembl_gene_id)
  #SIG <- colnames(reduced_elasnet_sigs$survive$ICU_survive_DE_UP$trainingData)[-1]
  SIG <- all_sigs$ICU_survive_DE_UP 
  SIG2 <- all_sigs$ICU_Mod_051
  
  # Check if samples are in the right order
  right_order <-  all(colnames(expr)[-1] == meta$sample_identifier)
  
  if (!right_order) {
    cat("Things are not in the right order")
    next
  }
  
  ### Normalize Expression data 
  expr_vst <- expr %>% 
    dplyr::select(one_of("ensembl_gene_id", meta$sample_identifier)) %>% 
    column_to_rownames(var = "ensembl_gene_id") %>% 
    as.matrix() %>% 
    varianceStabilizingTransformation() %>% 
    sva::ComBat(batch = meta$sequencing_month_year) %>% 
    as.data.frame()
  # Batch Correct?
  
  # Filter to genes of interest
  expr_filt <- expr_vst %>% 
    rownames_to_column(var = "ensembl_gene_id") %>% 
    filter(ensembl_gene_id %in% SIG2) %>% 
    left_join(universe, by = "ensembl_gene_id") %>% 
    dplyr::select(!one_of("ensembl_gene_id", "entrezgene_id", "description")) %>% 
    column_to_rownames(var = "hgnc_symbol")
  
  # Cluster
  expr_cluster <- expr_filt %>% 
    t() %>% scale(center = TRUE, scale = TRUE) %>%
    dist() %>% hclust() %>% cutree(k = 4) %>% 
    as.data.frame() %>% 
    set_names("cluster") %>% 
    rownames_to_column(var = "sample_identifier") %>% 
    mutate(Group = case_when(cluster == 1 ~ "One", cluster == 2 ~ "Two", cluster == 3 ~ "Three", cluster == 4 ~ "Four",
                             TRUE ~ NA_character_))
  
  # Plot
  meta_sig <- expr_cluster %>% 
    left_join(meta) %>% 
    dplyr::select(one_of("sample_identifier", "Group", "icu_adm", "survive", "culture", "sofa_24"))
  
  column_ha <- HeatmapAnnotation(
    SOFA = anno_barplot(meta_sig$sofa_24),
    Culture = meta_sig$culture,
    ICU_ADM = meta_sig$icu_adm,
    Mortality = meta_sig$survive,
    SignatureExprGroup = meta_sig$Group
                                 )
  icu_sig_clustered <- expr_filt %>% 
    t() %>% scale(center = TRUE, scale = TRUE) %>% t() %>% 
    ComplexHeatmap::Heatmap(column_split = meta_sig$Group,
      show_column_names = FALSE, 
                            show_row_names = FALSE, 
                            name = "Z-score",
                            top_annotation = column_ha)
  
  png("./train_test_mods/figures/ICU_MORTALITY_SIG_IN_ER.png", width = 900, height = 600)
  print(icu_sig_clustered)
  dev.off()
  
  meta %>% 
    left_join(expr_cluster) %>% 
    group_by(Group, icu_adm) %>% 
    summarize(n())
  meta %>% 
    left_join(expr_cluster) %>% 
    group_by(Group, culture) %>% 
    summarize(n())
  meta %>% 
    left_join(expr_cluster) %>% 
    group_by(Group, survive) %>% 
    summarize(n())
  meta %>% 
    left_join(expr_cluster) %>% 
    group_by(Group) %>% 
    summarize(mean(sofa_24, na.rm = TRUE))
  meta %>% 
    left_join(expr_cluster) %>% 
    group_by(Group) %>% 
    summarize(mean(sofa_72, na.rm = TRUE))
  
  # GSVA
  
  # ROAST
  
  # Predict it. 
  expr_filt %>% first_five()
  

  
  
  
}