# Get DE Genes correcting for cell proportions, age, gender, batch

rm(list = ls())

library(tidyverse)
library(magrittr)
library(functionjunction)
library(DESeq2)
library(viridis)
library(ggrepel)
library(UpSetR)
source("./helper_functions.R")

### Read in data
data <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_261120.RDS")
universe <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_261120.RDS")$universe
tr_te_dat <- read_rds("../paper_er_icu_sepsis/data_overview/tr_te_dat.RDS")
icu_dat <- read_rds("../paper_er_icu_sepsis/data_overview/icu_fir_sec_dat.RDS")

sev_model <- data.frame(
  model = c("severity", "severity", "severity", "severity"),
  group = rev(c("high", "intermediate", "low", "healthy_control")),
  color = rev(c("red", "darkorange1" ,"burlywood4", "black")),
  proper_name = rev(c("High", "Int", "Low", "HC"))
)

### Set up data
expr <- tr_te_dat$train$expr$raw %>% 
  left_join(icu_dat$first$expr$raw)

keep <- c("sample_identifier", "condition", "sequencing_month_year",  "outcome_icu_survival",  "outcome_mortality",
          "first_at_ed_sofa","worst_within_72_sofa", "outcome_sofa_second", "outcome_sofa_first", "outcome_sofa_third",
          "age", "gender", "micro_blood_culture_pathogen")
meta <- icu_dat$first$meta %>% 
  dplyr::select(one_of(keep)) %>% 
  dplyr::rename(sofa = "outcome_sofa_second", survive = "outcome_icu_survival") %>% 
  mutate(sofa_progress = ifelse(outcome_sofa_third - outcome_sofa_first > 0, "worse",
                                ifelse(outcome_sofa_third == outcome_sofa_first, "same", "improve"))) %>% 
  mutate(survive = ifelse(survive == 1, "survived", "dead"))
meta <- tr_te_dat$train$meta %>% 
  dplyr::select(one_of(keep)) %>% 
  dplyr::rename(sofa = "first_at_ed_sofa", survive = "outcome_mortality") %>% 
  mutate(sofa_progress = ifelse(worst_within_72_sofa - sofa > 0, "worse", 
                                ifelse(worst_within_72_sofa == sofa, "same", "improve"))) %>% 
  mutate(survive = ifelse(survive == 0, "survived", "dead")) %>% 
  bind_rows(meta) 

meta %<>% mutate(gender = ifelse(gender == 0, "M", "F")) 

### Remove controls
meta_no_hc <- meta %>% filter(!condition == "healthy_control")
expr_no_hc <- expr %>% dplyr::select(one_of("ensembl_gene_id", meta_no_hc$sample_identifier))
all(meta_no_hc$sample_identifier == colnames(expr_no_hc)[-1])

### Get severity groups
meta_no_hc %<>% 
  mutate(sofa_24_only_sev_group = 
           ifelse(condition == "healthy_control", "healthy_control", 
                  ifelse( sofa >= 5, "high", 
                          ifelse(sofa >= 2, "intermediate","low")))) %>% 
  mutate(sofa_24_only_sev_group = ifelse(is.na(sofa_24_only_sev_group), NA, .$sofa_24_only_sev_group)) %>% 
  mutate(sofa_24_only_sev_group = factor(sofa_24_only_sev_group, levels = c("low", "intermediate", "high")))

### Now make outcomes to model for DE. HighInt vs Low, High vs IntLow, Died vs Survive.
meta_no_hc %<>% 
  mutate(HighInt_Low = ifelse(is.na(sofa_24_only_sev_group) , NA, 
                              ifelse(sofa_24_only_sev_group %in% c("high", "intermediate"), "HighInt", "Low"))) %>% 
  mutate(HighInt_Low = factor(HighInt_Low, levels = c("Low", "HighInt"))) %>% 
  mutate(High_IntLow = ifelse(is.na(sofa_24_only_sev_group) , NA, 
                              ifelse(sofa_24_only_sev_group %in% c("high"), "High", "IntLow"))) %>% 
  mutate(High_IntLow = factor(High_IntLow, levels = c("IntLow", "High"))) %>%
  mutate(High_Low_BC = ifelse(is.na(sofa_24_only_sev_group), NA,
                              ifelse(sofa_24_only_sev_group == "high" & micro_blood_culture_pathogen == 1, "High_BC",
                                     ifelse(sofa_24_only_sev_group == "low" & micro_blood_culture_pathogen == 0, "Low_BC",
                                            ifelse(sofa_24_only_sev_group == "intermediate", NA, NA))))) %>% 
  mutate(High_Low_BC = factor(High_Low_BC, levels = c("Low_BC", "High_BC"))) %>% 
  mutate(survive = factor(survive, levels = c("survived", "dead"))) %>% 
  mutate(sofa_progress = factor(sofa_progress, levels = c("improve", "same", "worse")))

all(colnames(expr_no_hc)[-1] == meta_no_hc$sample_identifier)

### PCA
pca_meta <- meta_no_hc %>% filter(!is.na(sofa_24_only_sev_group))
pca_expr <- expr_no_hc %>% dplyr::select(one_of("ensembl_gene_id", pca_meta$sample_identifier))
pca <- pca_expr %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  as.matrix() %>% 
  DESeq2::varianceStabilizingTransformation()  %>% 
  sva::ComBat(batch = pca_meta$sequencing_month_year, mod  = model.matrix(~as.character(sofa_24_only_sev_group), data=pca_meta)) %>%
  t() %>% 
  perform_pca()

pca$x %>%
  left_join(meta_no_hc) %>% 
  ggplot(aes(x = PC1, y = PC2, color = sofa_24_only_sev_group, fill = sofa_24_only_sev_group)) +
  geom_point(size = 4, colour="black", pch=21 ) + 
  xlab(paste0("PC1 (", pca$pov[1]*100, "%)" )) + 
  ylab(paste0("PC2 (", pca$pov[2]*100, "%)" )) + 
  scale_fill_manual("", labels = as.character(sev_model$proper_name)[-1], values = as.character(sev_model$color)[-1] ) +
  scale_color_manual("", labels = as.character(sev_model$proper_name)[-1], values = as.character(sev_model$color)[-1]) +
  theme_minimal() +
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13))

### Add in cell proportion information
deconv_res <- read_rds("../paper_er_icu_sepsis/deconvolution/deconv_res.rds")
deconv_res_icu <- read_rds("../paper_er_icu_sepsis/deconvolution/icu_deconv_res.rds")

deconv_res_tr_icu_df <- deconv_res$train$cibersort %>% 
  left_join(deconv_res_icu$first$cibersort) %>% 
  dplyr::select(one_of("cell_type", meta_no_hc$sample_identifier)) %>% 
  mutate(cell_type = janitor::make_clean_names(cell_type)) %>% 
  column_to_rownames(var = "cell_type") %>% 
  t()
deconv_res_pca <- deconv_res_tr_icu_df %>% perform_pca()
deconv_res_pca$pov %>% explain_95()
deconv_res_pca_x <- deconv_res_pca$x %>% dplyr::select(one_of("sample_identifier", paste0("PC", 1:15)))

# Add into cell proportions to meta
meta_no_hc %<>% left_join(deconv_res_pca_x)

### DE
cols_to_keep <- c("sequencing_month_year", "age", "gender", paste0("PC", 1:15))
#outcomes <- c("sofa_progress","survive", "sofa_24_only_sev_group", "HighInt_Low", "High_IntLow")
outcomes <- c("High_Low_BC","sofa_24_only_sev_group", "HighInt_Low", "High_IntLow")


for (i in outcomes){
  numbers <- meta_no_hc %>% group_by_at(i) %>% summarize(n=n(), .groups = 'drop') 
  print(numbers)
  rm(numbers)
}

res_de <- list() ; res_pthwy <- list()
for (i in outcomes) {
  cat(paste0("Outcome: ", i, "\n"))
  
  met <- meta_no_hc %>% 
    dplyr::select(one_of("sample_identifier", i, cols_to_keep)) %>% 
    dplyr::rename(comparison = i) %>% 
    filter(!is.na(comparison))
  
  exp <- expr_no_hc %>% 
    dplyr::select(one_of("ensembl_gene_id", met$sample_identifier)) %>% 
    column_to_rownames(var = "ensembl_gene_id")
  
  right_order <-  all(colnames(exp) == met$sample_identifier)
  
  if (!right_order) {
    cat("Things are not in the right order")
    stop()
  }
  
  de_res <- de(counts = exp, meta = met, FC = 1, 
               des = paste0("comparison + ", paste(cols_to_keep, collapse = " + ")), 
               main_covar = "comparison", filt = FALSE)
  
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
  res_de[[i]] <- de_res
  res_pthwy[[i]] <- list(
    all = pthwy,
    up = pthwy_up,
    down = pthwy_down)
}

# Save results
res_de %>% write_rds(paste0("./results/severity_combine_de_tr_icu_res.rds"))
res_pthwy %>% write_rds(paste0("./results/severity_combine_pthwy_enr_tr_icu_res.rds"))

# Get numbers
res_de %>% 
  map(~map(.x, ~filter(.x, de == "de"))) %>% 
  map(~map(.x, ~nrow(.x)))
res_de %>% 
  map(~map(.x, ~filter(.x, de == "de" & direction == "up"))) %>% 
  map(~map(.x, ~nrow(.x)))


# Save results
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
