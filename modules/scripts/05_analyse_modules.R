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

### Read in modules, module enrichment
ICA_res  <- read_rds("./modules/results/ICA_res.rds")
ICA_mod_plot <- read_rds("./modules/results/ICA_mod_plot.rds")

######### Lets plot results for the corrected data
### How many modules are there?
ICA_res$expr_cor$ICA_mods_filt %>% length()

### Look at module overlap
ICA_jac_clust <- ICA_res %>% 
  map(~get_jac_mat(.x$ICA_mods_filt) %>% 
        t() %>% 
        vegan::vegdist(method = "jaccard", binary = TRUE, diag = TRUE) %>% 
        hclust)

png("./modules/figures/modules_expr_cor_jaccard_clust.png", width = 900, height = 400, units = "px")
plot(ICA_jac_clust$expr_cor)
dev.off()

### Look at module pathways
ggsave("./modules/figures/modules_expr_cor_reactome.png", ICA_mod_plot$expr_cor$reactome, width = 10, height = 6 )
ggsave("./modules/figures/modules_expr_cor_msigdb.png",ICA_mod_plot$expr_cor$msigdb, width = 10, height = 6 )

### Look at clinical associations 
clinical_assoc <-  read_rds("./modules/results/ICA_mod_eigen_clin_assoc.rds")
clinical_assoc_filt <-  read_rds("./modules/results/ICA_mod_eigen_clin_assoc_filt_plot.rds")

# Plot
ggsave("./modules/figures/modules_expr_cor_sofa24_eigengene.png", 
       clinical_assoc_filt$expr_cor$plot$first_at_ed_sofa, height = 2, width = 10, scale = 1.5)
ggsave("./modules/figures/modules_expr_cor_sofa72_eigengene.png", 
       clinical_assoc_filt$expr_cor$plot$worst_within_72_sofa, height = 2, width = 10, scale = 1.5 )

# Look at modules which are correlated to more outcomes
mods_to_keep <- clinical_assoc_filt$expr_vst$full_res %>% 
  map(~.x$comp) %>% 
  flatten_chr() %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 3) %>% 
  pull(".") %>% 
  as.character() %>% 
  str_replace("_", " ")

ICA_mod_plot$expr_cor$reactome$data %>% 
  filter(comp %in% mods_to_keep) %>% View()
ICA_mod_plot$expr_cor$msigdb$data %>% 
  filter(comp %in% mods_to_keep) %>% View()

### Plot heatmaps
meta_htmap <- meta %>% arrange(desc(outcome_icu_admission))

column_ha <- HeatmapAnnotation(SOFA24 = anno_barplot(meta_htmap$first_at_ed_sofa),
                               SOFA72 = anno_barplot(meta_htmap$worst_within_72_sofa),
                               Culture = meta_htmap$micro_blood_culture_pathogen,
                               ICU_Adm = meta_htmap$outcome_icu_admission
)

temp <- expr_mods$expr_vst %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  dplyr::select(one_of("ensembl_gene_id", meta_htmap$sample_identifier)) %>% 
  filter(ensembl_gene_id %in% ICA_res$expr_vst$ICA_mods_filt$Mod_34$ensembl_gene_id) %>% 
  left_join(universe) %>% 
  dplyr::select(!one_of("ensembl_gene_id", "description", "entrezgene_id")) %>% 
  filter(hgnc_symbol != "") %>% 
  column_to_rownames(var = "hgnc_symbol")
all(colnames(temp) == meta_htmap$sample_identifier)

png("./modules/figures/modules_expr_vst_mod_9.png", width = 900, height = 600)
print(ComplexHeatmap::Heatmap(as.matrix(temp),
                        name = "Expression",
                        show_column_names = FALSE,
                        show_row_names = FALSE,
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = TRUE,
                        cluster_column_slices = FALSE,
                        border = TRUE, 
                        top_annotation = column_ha
))
dev.off()

