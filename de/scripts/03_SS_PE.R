rm(list = ls())
source("./misc/helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

expr <- tr_te_dat$er_tr$expr
meta <- tr_te_dat$er_tr$meta 
all(colnames(expr)[-1] == meta$sample_identifier)

### Add in Endotype classification
endotype_res <- read_rds("../paper_er_icu_sepsis/endotype/consensus_clust_res/kmed_endotype_all_res_clusters.rds")
endotype_res <- endotype_res$top_5_percent$k_5 %>% 
  dplyr::select(one_of("sample_identifier", "cluster")) %>% 
  mutate(endotype2class = ifelse(cluster %in% c("cluster_2", "cluster_3"), "NPS_INF", "IHD_IFN_ADA")) %>% 
  mutate(endotype2class = factor(endotype2class, levels = c("IHD_IFN_ADA", "NPS_INF")))
meta %<>% left_join(endotype_res, by = "sample_identifier")

### Read in DE Results
res <- read_rds("./de/results/SS_DE_res.rds")

### Plot pathways 
<<<<<<< HEAD
outcomes <- c("culture", "High_Low", "icu_adm", "survive")
comparison_names <- c("Pos vs\n Neg", "High vs\n Neg", "ICU vs\n Non-ICU", "Dead vs\n Survived")
names(comparison_names) <-  c("pos_vs_neg", "high_vs_low",  "icu_vs_non_icu", "dead_vs_survive")

=======
outcomes <- c("endotype2class", "survive","icu_adm", "culture", "qsofa_High_Low", "High_Low")
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff

pthwy_plot <- res %>% 
  map(~.x$pthwy) %>% 
  map(~map(.x, ~bind_rows(.x, .id = "direction"))) %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
<<<<<<< HEAD
  bind_rows(.id = "outcome") %>% 
=======
  bind_rows() %>% 
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  separate(BgRatio,into = c("M", "N")) %>%
  mutate(Ratio = Count/as.numeric(M) ) %>% 
  left_join(pathway_hier,  by = c("ID" = "enr_pathway")) 

<<<<<<< HEAD
pthwy_plot$outcome %>% unique()
pthwy_plot %>% 
  filter(outcome %in% outcomes) %>% 
  group_by(comparison, direction, one_lower_level_pathway_descrip_clean_name ) %>% 
  top_n(3, Ratio) %>% 
  mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down"))) %>% 
  mutate(comparison = factor(comparison, levels =  names(comparison_names), labels =  comparison_names)) %>% 
=======
pthwy_plot %>% 
  #filter(top_level_pathway_descrip == pthwy_group) %>% 
  mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down"))) %>% 
  # mutate(comparison = factor(comparison, levels =  names(comparison_names), labels =  comparison_names)) %>% 
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  mutate(one_lower_level_pathway_descrip_clean_name = str_wrap(one_lower_level_pathway_descrip_clean_name, width = 20)) %>% 
  mutate(Description = str_wrap(Description, width = 60)) %>% 
  ggplot(aes(x = direction, y = Description, fill = Ratio)) + 
  geom_tile() + 
<<<<<<< HEAD
  facet_grid(cols = vars(comparison), rows = vars(top_level_pathway_descrip, one_lower_level_pathway_descrip_clean_name), 
             scale = "free_y", space = "free_y") +
=======
  facet_grid(cols = vars(comparison), rows = vars(top_level_pathway_descrip, one_lower_level_pathway_descrip_clean_name), scale = "free_y", space = "free_y") +
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA), plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(angle = 0)) + 
  scale_fill_viridis(direction = -1) + 
  ylab("") + xlab("")
<<<<<<< HEAD
ggsave("./de/figures/SS_pthwy_enr.png", scale = 0.9, width = 13, height = 8)
=======
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
