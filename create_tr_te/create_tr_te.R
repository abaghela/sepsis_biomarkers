source("./helper_functions.R")

### Read in data 
data <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_261120.RDS")
universe <- data$universe

expr <- data$expr$raw_unfilt$samples_merged
meta <- data$meta$samples_merged

# Remove samples of interest
meta %<>% 
  filter(!sample_identifier %in% c("sepnet1324b", "sepnet1214b", "sepnet1179b", "sepnet155b") & time_point != "second") 

# Put everything in order
expr %<>% dplyr::select(one_of("ensembl_gene_id", meta$sample_identifier))
all(colnames(expr)[-1] == meta$sample_identifier)

# Filter Samples & Genes
remove <- "novel transcript|novel protein|Novel transcript|Novel|novel|immunoglobulin lambda constant|immunoglobulin lambda variable|immunoglobulin lambda joining|immunoglobulin heavy constant|immunoglobulin heavy variable|immunoglobulin kappa constant|immunoglobulin kappa variable|immunoglobulin kappa joining |non-protein|pseudogene|non-functional|inactive|TEC|TEC gene|artifact|uncharacterized"
universe_filt <- universe %>% filter(!grepl(remove, description)) %>% filter(!is.na(entrezgene_id))

expr %<>% 
  filter(ensembl_gene_id %in% universe_filt$ensembl_gene_id) %>% 
  column_to_rownames(var = "ensembl_gene_id") %>% 
  remove_low_count_genes() %>% 
  sample_lib_size_filter() %>% 
  rownames_to_column(var = "ensembl_gene_id")
  
meta %>% 
  filter(sample_identifier %in% colnames(expr)) %>% 
  group_by(condition, sample_location, sequencing_month_year, time_point) %>% summarize(n = n())

all(colnames(expr)[-1] == meta$sample_identifier)

#### ER Train 
er_train_meta <- meta %>% 
  filter(sample_location %in% c("australia_2", "colombia", "houston", "netherlands", "vancouver", "hancock_lab")) 
er_train_expr <- expr %>% 
  dplyr::select(one_of("ensembl_gene_id", er_train_meta$sample_identifier))
er_train_meta %>%  group_by(condition, sample_location, sequencing_month_year, time_point) %>% summarize(n = n())
all(colnames(er_train_expr)[-1] == er_train_meta$sample_identifier)

#### ER Test
er_test_meta <- meta %>% 
  filter(sample_location %in% c("australia")) 
er_test_expr <- expr %>% 
  dplyr::select(one_of("ensembl_gene_id", er_test_meta$sample_identifier))
er_test_meta %>%  group_by(condition, sample_location, sequencing_month_year, time_point) %>% summarize(n = n())
all(colnames(er_test_expr)[-1] == er_test_meta$sample_identifier)

#### ICU
icu_meta_full <- data$meta_icu_covid
icu_meta <- meta %>% 
  dplyr::select(one_of("sample_identifier", "condition", "sequencing_month_year", "sample_location")) %>% 
  filter(sample_location %in% c("toronto")) %>% 
  left_join(icu_meta_full)
icu_expr <- expr %>% 
  dplyr::select(one_of("ensembl_gene_id", icu_meta$sample_identifier ))
icu_meta %>%  group_by(condition, sample_location, sequencing_month_year, time_point) %>% summarize(n = n())
all(colnames(icu_expr)[-1] == icu_meta$sample_identifier)

tr_te_dat <- list(
  er_tr = list(expr = er_train_expr, meta = er_train_meta),
  er_te = list(expr = er_test_expr, meta = er_test_meta),
  icu = list(expr = icu_expr, meta = icu_meta)
  )
tr_te_dat %>% map(~all(colnames(.x$expr)[-1] == .x$meta$sample_identifier))
rm(er_train_expr, er_train_meta, er_test_expr, er_test_meta, icu_meta_full, icu_meta, icu_expr)

tr_te_dat %>% write_rds("./create_tr_te/tr_te_dat.RDS")
