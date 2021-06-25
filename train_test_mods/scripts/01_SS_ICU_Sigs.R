# Script to combine all signatures such that they can be used for enrichment statistics or prediction in the SS or ICU validation cohorts.
rm(list = ls())

#### Read in signatures
# ICU
ICU_ICA <-  read_rds("./modules/results/ICU_ICA_mod_sigs.rds")$expr_uncorr
ICU_ICA_sigs <- ICU_ICA %>% 
  dplyr::select(one_of("sig_name_no_clin", "mod_genes")) %>% 
  distinct()
ICU_ICA_sigs <- ICU_ICA_sigs$mod_genes %>% 
  set_names(ICU_ICA_sigs$sig_name_no_clin)

ICU_DE <- read_rds("./de/results/ICU_DE_res.rds") %>% 
  map(~map(.x$de, ~filter(.x, de == "de")))

ICU_DE_ALL <- ICU_DE %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  map(~.x$ensembl_gene_id) %>% 
  set_names(paste0("ICU_",names(.), "_DE_ALL"))

ICU_DE_UP <- ICU_DE %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  map(~filter(.x, direction == "up")) %>% 
  map(~.x$ensembl_gene_id) %>% 
  set_names(paste0("ICU_",names(.), "_DE_UP"))

ICU_DE_DN <- ICU_DE %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  map(~filter(.x, direction == "down")) %>% 
  map(~.x$ensembl_gene_id) %>% 
  set_names(paste0("ICU_",names(.), "_DE_DN"))

# SS
SS_ICA <- read_rds("./modules/results/SS_ICA_mod_sigs.rds")$expr_uncorr
SS_ICA_sigs <- SS_ICA %>% 
  dplyr::select(one_of("sig_name_no_clin", "mod_genes")) %>% 
  distinct()
SS_ICA_sigs <- SS_ICA_sigs$mod_genes %>% 
  set_names(SS_ICA_sigs$sig_name_no_clin)

SS_DE <- read_rds("./de/results/SS_DE_res.rds") %>% 
  map(~map(.x$de, ~filter(.x, de == "de")))

SS_DE_ALL <- SS_DE %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  map(~.x$ensembl_gene_id) %>% 
  set_names(paste0("SS_",names(.), "_DE_ALL"))

SS_DE_UP <- SS_DE %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  map(~filter(.x, direction == "up")) %>% 
  map(~.x$ensembl_gene_id) %>% 
  set_names(paste0("SS_",names(.), "_DE_UP"))

SS_DE_DN <- SS_DE %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  map(~filter(.x, direction == "down")) %>% 
  map(~.x$ensembl_gene_id) %>% 
  set_names(paste0("SS_",names(.), "_DE_DN"))

all_sigs <- ICU_ICA_sigs %>% 
  append(SS_ICA_sigs) %>% 
  append(ICU_DE_ALL) %>% 
  append(SS_DE_ALL) %>% 
  append(ICU_DE_UP) %>% 
  append(SS_DE_UP) %>% 
  append(ICU_DE_DN) %>% 
  append(SS_DE_DN) 

all_sigs %>% write_rds("./train_test_mods/results/DE_ICA_sigs.RDS")
