rm(list = ls())
source("./misc/helper_functions.R")

### Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

expr <- tr_te_dat$er_tr$expr
meta <- tr_te_dat$er_tr$meta 
all(colnames(expr)[-1] == meta$sample_identifier)

### Read in Modules
SS_ICA_mod_plot <- read_rds("./modules/results/SS_ICA_mod_plot.rds")
ICU_ICA_mod_plot <- read_rds("./modules/results/ICU_ICA_mod_plot.rds")

### Which signatures are worth of exploration
SS_validated_sigs <- read_rds("./validation/results/SS_validation_enrichment.rds")
SS_validated_sigs %>% 
  map(~top_n(.x, -5, PValue)) %>% View()


SS_ICA_mod_plot$expr_uncorr$reactome
SS_ICA_mod_plot$expr_uncorr$msigdb

ICU_ICA_mod_plot$expr_uncorr$reactome
ICU_ICA_mod_plot$expr_uncorr$msigdb

hey <- SS_validated_sigs %>% 
  map(~filter(.x, PValue <= 0.05))
hey$culture %>% 
  rownames_to_column(var = "signature") %>% 
  ggplot(y = signature, x )
hey$High_Low
