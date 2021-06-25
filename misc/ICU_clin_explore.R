rm(list = ls())
source("./misc/helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

expr <- tr_te_dat$icu$expr
meta <- tr_te_dat$icu$meta 
all(colnames(expr)[-1] == meta$sample_identifier)

### Severity - Outcomes Progression 
alluv_plt <- meta %>% 
  dplyr::select(one_of("sample_identifier","sofa_sev_24", "survive")) %>% 
  group_by(sofa_sev_24, survive) %>% 
  summarize(Freq = n()) %>% 
  ungroup() %>% 
  na.omit()

png("./misc/figures/ICU_Sev_Outcome_Alluvial.png", width = 400, height = 300, units = "px")
print(alluvial::alluvial(alluv_plt[,c(1,2)], freq = alluv_plt$Freq))
dev.off()

### Severity - Outcomes Association 
broom::tidy(lm(sofa_24 ~ survive, data = meta))
meta %>% 
  dplyr::select(one_of("sample_identifier", "sofa_24", "survive")) %>% 
  na.omit() %>% 
  ggplot(aes(x = survive, y = sofa_24)) + 
  geom_boxplot()
ggsave("./misc/figures/ICU_SOFA_outcomes.png", scale = 0.9)
