rm(list = ls())
source("./helper_functions.R")

### Read in differential expression/pathway enrichment results
res <- read_rds("./de/results/de_tr_res.rds")

### Combine results and Plot
res_pthwy_df <- res %>% 
  map(~.x$pthwy) %>% 
  map(~map(.x, ~bind_rows(.x, .id = "comparison"))) %>% 
  map(~bind_rows(.x, .id = "direction")) %>% 
  bind_rows(.id = "outcome")

res_msigdb_df <- res %>% 
  map(~.x$msigdb) %>% 
  map(~map(.x, ~map(.x, ~bind_rows(.x, .id = "database")))) %>% 
  map(~map(.x, ~bind_rows(.x, .id = "comparison"))) %>% 
  map(~bind_rows(.x, .id = "direction")) %>% 
  bind_rows(.id = "outcome")

### Plot 
# Outcomes: Reactome
keep = list("icu_vs_non_icu" = "ICU vs\n non-ICU", 
            "dead_vs_survive" = "Dead vs\n Survive" )
res_pthwy_df %>% 
  filter(comparison %in% names(keep)) %>% 
  mutate(comparison = factor(comparison, levels = names(keep),  labels = flatten_chr(keep))) %>% 
  plot_reactome_pathway()
# Outcomes: MSigDB
res_msigdb_df %>% 
  filter(database == "MSigDB_Hallmark_2020") %>% 
  filter(comparison %in% names(keep)) %>% 
  mutate(comparison = factor(comparison, levels = names(keep),  labels = flatten_chr(keep))) %>% 
  plot_msigdb_pathways()

# Severity: Reactome
keep = list("high_int_vs_low" = "High + Int vs\n Low", 
            #"pos_vs_neg" = "BC Pos vs\n Neg",
            "high_int_pos_vs_low_neg" = "High + Int (BC Pos) vs\n Low (BC Neg)")
res_pthwy_df %>% 
  filter(comparison %in% names(keep)) %>% 
  mutate(comparison = factor(comparison, levels = names(keep),  labels = flatten_chr(keep))) %>% 
  plot_reactome_pathway()
# Severity: MSigDB
res_msigdb_df %>% 
  filter(database == "MSigDB_Hallmark_2020") %>% 
  filter(comparison %in% names(keep)) %>% 
  mutate(comparison = factor(comparison, levels = names(keep),  labels = flatten_chr(keep))) %>% 
  plot_msigdb_pathways()

# Severity Extreme Phenotype: Reactome
keep = list("high_vs_low" = "High vs\n Low", 
            #"pos_vs_neg" = "BC Pos vs\n Neg",
            "high_pos_vs_low_neg" = "High (BC Pos) vs\n Low (BC Neg)")
res_pthwy_df %>% 
  filter(comparison %in% names(keep)) %>% 
  mutate(comparison = factor(comparison, levels = names(keep),  labels = flatten_chr(keep))) %>% 
  plot_reactome_pathway()
# Severity: MSigDB
res_msigdb_df %>% 
  filter(database == "MSigDB_Hallmark_2020") %>% 
  filter(comparison %in% names(keep)) %>% 
  mutate(comparison = factor(comparison, levels = names(keep),  labels = flatten_chr(keep))) %>% 
  plot_msigdb_pathways()

"high_pos_vs_low_neg" = "High (BC Pos) vs\n Low (BC Neg)"


labels = c("Pos vs\n Neg","ICU vs\n non-ICU", "Dead vs\n Survive")
### Volcano plot 
# res_de$sofa_24_only_sev_group$high_low %>% vp(colour = "red", labs_to_incl = 10, fc_cutoff = 1.09, title = "High vs Low")
# ggsave("./de/figures/de_high_vs_low.tiff", scale = 1.2)
# ggsave("./de/figures/de_high_vs_low.png", scale = 1.2)
