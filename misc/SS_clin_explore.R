rm(list = ls())
source("./misc/helper_functions.R")

### Read in data - remove healthy controls
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
deconv_res <- read_rds("./create_tr_te/deconv_res.rds") %>% 
  map(~map(.x, ~as.data.frame(.x)))

expr <- tr_te_dat$er_tr$expr
meta <- tr_te_dat$er_tr$meta 
all(colnames(expr)[-1] == meta$sample_identifier)

### Severity - Outcomes Association 
sev_outcome_assoc <- meta %>% 
  dplyr::select(one_of("sample_identifier", "sofa_24", "at_ed_qsofa", "icu_adm", "survive" ))

broom::tidy(lm(at_ed_qsofa ~ survive, data = sev_outcome_assoc))
broom::tidy(lm(at_ed_qsofa ~ icu_adm, data = sev_outcome_assoc))

broom::tidy(lm(sofa_24 ~ survive, data = sev_outcome_assoc))
broom::tidy(lm(sofa_24 ~ icu_adm, data = sev_outcome_assoc))

sev_outcome_assoc %>% 
  pivot_longer(cols = c("survive", "icu_adm"), "outcome") %>% 
  filter(!is.na(sofa_24)) %>% 
  filter(!is.na(value)) %>% 
<<<<<<< HEAD
  ggplot(aes(x = value, y = sofa_24)) + geom_boxplot() + 
  facet_grid(cols = vars(outcome), scales = "free") + 
  ylab("SOFA") + xlab("")
ggsave("./misc/figures/SS_SOFA_outcomes.png", scale = 0.9)

### Severity - Outcomes Progression 
alluv_plt <- meta %>% 
  dplyr::select(one_of("sample_identifier", "culture", "sofa_sev_24", "icu_adm", "survive")) %>% 
  group_by( culture, sofa_sev_24, icu_adm, survive) %>% 
=======
  ggplot(aes(x = value, y = sofa_24)) + geom_boxplot() + facet_grid(cols = vars(outcome), scales = "free") 

### Severity - Outcomes Progression 
alluv_plt <- meta %>% 
  dplyr::select(one_of("sample_identifier","qsofa_sev_adm", "culture", "sofa_sev_24", "icu_adm", "survive")) %>% 
  group_by(qsofa_sev_adm, culture, sofa_sev_24, icu_adm, survive) %>% 
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff
  summarize(Freq = n()) %>% 
  ungroup() %>% 
  na.omit()

<<<<<<< HEAD
png("./misc/figures/SS_Sev_Outcome_Alluvial.png", width = 700, height = 300, units = "px")
print(alluvial::alluvial(alluv_plt[,c(1:4)], freq = alluv_plt$Freq))
dev.off()
=======
alluvial::alluvial(alluv_plt[,c(1:5)], freq = alluv_plt$Freq)
>>>>>>> 36a7168835f86e6a246bd8b4e62b095b5923a0ff

# DE
meta_de <- meta %>% 
  mutate(low_sev_outcome_poor = ifelse(sofa_sev_24 %in% c("low") & (icu_adm == "icu" | survive == "dead" ), "low_sev_poor_outcome",
                                ifelse(sofa_sev_24 %in% c("low") & (icu_adm == "non_icu" & survive == "survive" ), "low_sev_fair_outcome", NA))) %>% 
  dplyr::select(one_of("sample_identifier", "sofa_sev_24", "icu_adm","survive", "low_sev_outcome_poor", "age", "gender", "sequencing_month_year")) %>% 
  filter(!is.na(low_sev_outcome_poor))

deconv_pca <- deconv_res$er_tr$cibersort %>% 
    rownames_to_column(var = "sample_identifier") %>% 
    filter(sample_identifier %in% meta_de$sample_identifier) %>% 
    remove_zero_var_cols() %>% 
    column_to_rownames(var = "sample_identifier") %>% 
    perform_pca() 
  
meta_de <- meta_de %>% left_join(deconv_pca$x,  by = "sample_identifier")
  
pcs_to_incl <- deconv_pca$pov %>% explain_95(perc = 90)
cat(paste0("Cell Proportion PCs used: ", pcs_to_incl, "\n"))
pcs_to_incl <- paste0("PC", 1:pcs_to_incl)
  
des <- paste0( "~ low_sev_outcome_poor + age + gender + sequencing_month_year", " + ",paste(pcs_to_incl, collapse = " + "))
  
expr_de <- expr %>% 
  dplyr::select(one_of("ensembl_gene_id", meta_de$sample_identifier)) %>% 
  column_to_rownames(var = "ensembl_gene_id")

all(colnames(expr_de) == meta_de$sample_identifier)

de_res <- de(counts = expr_de, meta_de, des = des, main_covar = "low_sev_outcome_poor")

de_res$low_sev_poor_outcome_vs_low_sev_fair_outcome %>% filter(de == "de") %>% View()
pthwy <- de_res %>% 
  map(~filter(.x, de == "de")) %>% 
  map(~split(.x, .x$direction)) %>% 
  map(~map(.x, ~pathway_enrichment(.x$entrezgene_id, p_val = 0.05, universe_list = universe$entrezgene_id)))
msigdb <-  de_res %>% 
  map(~filter(.x, de == "de")) %>% 
  map(~split(.x, .x$direction)) %>% 
  map(~map(.x, ~go_enrichment(.x$hgnc_symbol, p_val = 0.05))) 


pthwy$low_sev_poor_outcome_vs_low_sev_fair_outcome$down
pthwy$low_sev_poor_outcome_vs_low_sev_fair_outcome$up

msigdb$low_sev_poor_outcome_vs_low_sev_fair_outcome$down$MSigDB_Hallmark_2020
msigdb$low_sev_poor_outcome_vs_low_sev_fair_outcome$up$MSigDB_Hallmark_2020
