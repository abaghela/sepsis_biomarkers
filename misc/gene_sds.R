rm(list = ls())
source("./helper_functions.R")

### Read in data
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")

read_proc_sd <- function(cohort){
  
  # Read in 
  expr <- tr_te_dat[[cohort]]$expr 
  meta <- tr_te_dat[[cohort]]$meta 
  all(colnames(expr)[-1] == meta$sample_identifier)
  
  # VST
  process <- expr %>% 
    column_to_rownames(var = "ensembl_gene_id") %>% 
    as.matrix() %>% 
    DESeq2::varianceStabilizingTransformation() 
  
  # Combat
  batches <- length(unique(meta$sample_location))
  
  if (batches > 1) {
    process <- process %>% 
      sva::ComBat(batch = meta$sequencing_month_year) 
  }

  # Get SD 
  sd <- process %>% 
    t() %>% 
    as.data.frame() %>% 
    map_df(~sd(.x))
  
  return(sd)

}

res <- c("er_tr", "icu", "hc") %>% 
  map(~read_proc_sd(cohort = .x)) 
names(res) <- c("er_tr", "icu", "hc")

# PLot
res %>% 
  map(~t(.x) %>% as.data.frame() %>% rownames_to_column(var = "ensembl_gene_id")) %>% 
  bind_rows(.id = "group") %>% 
  mutate(group = factor(group, levels = c("hc", "er_tr", "icu"), labels = c("HC", "SS", "ICU"))) %>% 
  ggplot(aes(x = V1, fill = group)) + geom_histogram(color="#e9ecef", alpha=0.8, position = 'identity') +
  xlab("Standard Deviation")
ggsave("./de/figures/ss_icu_hc_gene_sd.png", width = 8, height = 3)


