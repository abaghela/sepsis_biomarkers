rm(list = ls())
source("./misc/helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

expr <- tr_te_dat$icu$expr
meta <- tr_te_dat$icu$meta 
all(colnames(expr)[-1] == meta$sample_identifier)

### Read in DE Results
res <- read_rds("./de/results/ICU_DE_res.rds")

### Plot pathways 
# outcomes <- c( "High_Low", "survive")
# comparison_names <- c( "High vs\n Neg",  "Dead vs\n Survived")
# names(comparison_names) <-  c( "high_vs_low",   "dead_vs_survive")

pthwy_plot <- res %>% 
  map(~.x$pthwy) %>% 
  map(~map(.x, ~bind_rows(.x, .id = "direction"))) %>% 
  map(~bind_rows(.x, .id = "comparison")) %>% 
  bind_rows(.id = "outcome") %>% 
  separate(BgRatio,into = c("M", "N")) %>%
  mutate(Ratio = Count/as.numeric(M) ) %>% 
  left_join(pathway_hier,  by = c("ID" = "enr_pathway")) 

pthwy_plot %>% 
  #filter(outcome %in% outcomes) %>% 
  # group_by(comparison, direction, one_lower_level_pathway_descrip_clean_name ) %>% 
  # top_n(3, Ratio) %>% 
  mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down"))) %>% 
  #mutate(comparison = factor(comparison, levels =  names(comparison_names), labels =  comparison_names)) %>% 
  mutate(one_lower_level_pathway_descrip_clean_name = str_wrap(one_lower_level_pathway_descrip_clean_name, width = 20)) %>% 
  mutate(Description = str_wrap(Description, width = 60)) %>% 
  ggplot(aes(x = direction, y = Description, fill = Ratio)) + 
  geom_tile() + 
  facet_grid(cols = vars(outcome), rows = vars(top_level_pathway_descrip, one_lower_level_pathway_descrip_clean_name), scale = "free_y", space = "free_y") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA), plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(angle = 0)) + 
  scale_fill_viridis(direction = -1) + 
  ylab("") + xlab("")
#ggsave("./de/figures/ICU_pthwy_enr.png", scale = 0.9, width = 11, height = 5)


