rm(list = ls())

library(tidyverse)
library(magrittr)
library(functionjunction)
universe <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_261120.RDS")$universe

res_de <-  read_rds(paste0("./results/severity_combine_de_tr_icu_res.rds"))

top_pos <- res_de %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  bind_rows(.id = "outcome") %>% 
  filter(de == "de" & direction == "up") %>% 
  group_by(comparison) %>% 
  top_n(300, fc) %>% 
  mutate(comparison = paste0(comparison, "_POS")) %>% 
  split(.$comparison) %>% 
  keep(names(.) %in% c("high_low_POS", "High_IntLow_POS", "High_BC_Low_BC_POS") )

top_pos_neg <-  res_de %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  bind_rows(.id = "outcome") %>% 
  filter(de == "de" ) %>% 
  group_by(comparison) %>% 
  top_n(300, abs(fc)) %>% 
  split(.$comparison) %>% 
  keep(names(.) %in% c("high_low", "High_IntLow", "High_BC_Low_BC") )

genes_to_use <- top_pos %>% 
  append(top_pos_neg)
genes_to_use %<>% map(~.x$ensembl_gene_id)
genes_to_use %>% write_rds("./results/severity_signatures.rds")
