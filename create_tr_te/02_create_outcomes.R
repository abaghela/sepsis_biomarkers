rm(list = ls())
source("./helper_functions.R")

### Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")

# The ER cohorts have the same metadata/clinical profiles 
all(colnames(tr_te_dat$er_tr$meta) == colnames(tr_te_dat$er_te$meta))

# The ICU cohort has slightly different column names for the same data 
all(colnames(tr_te_dat$er_tr$meta) == colnames(tr_te_dat$icu$meta))

##### Create outcomes and clean up covariates of interest
## Create column names which are common between ER and ICU patients 
# ER
tr_te_dat$er_tr$meta <- tr_te_dat$er_tr$meta %>% 
  dplyr::rename(sofa_24 = "first_at_ed_sofa", sofa_72 = "worst_within_72_sofa", 
                survive = "outcome_mortality", icu_adm = "outcome_icu_admission",
                culture = "micro_blood_culture_pathogen")
tr_te_dat$er_te$meta <- tr_te_dat$er_te$meta %>% 
  dplyr::rename(sofa_24 = "first_at_ed_sofa", sofa_72 = "worst_within_72_sofa", 
                survive = "outcome_mortality", icu_adm = "outcome_icu_admission",
                culture = "micro_blood_culture_pathogen")
# ICU
tr_te_dat$icu$meta <- tr_te_dat$icu$meta %>% 
  dplyr::rename(sofa_24 = "outcome_sofa_second", sofa_72 = "outcome_sofa_third", 
                survive = "outcome_icu_mortality", culture = "micro_blood_culture_pathogen") %>% 
  dplyr::mutate(icu_adm = ifelse(patient_location == "icu", 1, 0)) 

## Write function to create new column names 

create_outcomes <- function(meta) {
  
  # Gender
  meta %<>% 
    mutate(gender = ifelse(gender == 0, "M", "F"))
  
  # Outcomes: ICU admission, blood culture, survival
  meta %<>% 
    # Survival
    mutate(survive = ifelse(survive == 1, "dead", "survive")) %>% 
    mutate(surive = factor(survive, levels = c("survive", "dead"))) %>% 
    # ICU admission
    mutate(icu_adm = ifelse(icu_adm == 1, "icu", "non_icu")) %>% 
    mutate(icu_adm = factor(icu_adm, levels = c("non_icu", "icu"))) %>% 
    # Blood Culture
    mutate(culture = ifelse(culture == 1, "pos", "neg")) %>% 
    mutate(culture = factor(culture, levels = c("neg", "pos")))
  
  
  # Severity: SOFA scores
  meta %<>% 
    # Severity @ 24H
    mutate(sofa_sev_24 = 
             ifelse(is.na(sofa_24), NA, 
                    ifelse(sofa_24 >= 5, "high", ifelse(sofa_24 >= 2 & sofa_24 < 5, "int", "low")))) %>% 
    mutate(sofa_sev_24 = factor(sofa_sev_24, levels = c("low", "int", "high"))) %>% 
    # Severity Progression from 24H to 72H
    mutate(sofa_sev_prog = 
             ifelse(is.na(sofa_72), NA,
                    ifelse(sofa_72 - sofa_24 > 0, "worse", 
                                  ifelse(sofa_72 == sofa_24, "same", "improve")))) %>% 
    mutate(sofa_sev_prog = factor(sofa_sev_prog, levels = c("improve", "same", "worse")))
  
  # Composite scores
  meta %<>% 
    # High + Int vs Low
    mutate(HighInt_Low = ifelse(is.na(sofa_sev_24) , NA, 
                                ifelse(sofa_sev_24 %in% c("high", "int"), "high_int", "low"))) %>% 
    mutate(HighInt_Low = factor(HighInt_Low, levels = c("low", "high_int"))) %>% 
    # High vs Int + Low
    mutate(High_IntLow = ifelse(is.na(sofa_sev_24) , NA, 
                                ifelse(sofa_sev_24 %in% c("high"), "high", "int_low"))) %>% 
    mutate(High_IntLow = factor(High_IntLow, levels = c("int_low", "high"))) %>%
    # High vs Low
    mutate(High_Low = ifelse(is.na(sofa_sev_24) , NA, 
                                ifelse(sofa_sev_24 %in% c("high"), "high", 
                                       ifelse(sofa_sev_24 %in% c("low"), "low", NA)))) %>% 
    mutate(High_Low = factor(High_Low, levels = c("low", "high"))) %>% 
    # High + Int (Culture Pos) vs Low (Culture Neg)
    mutate(HighInt_Low_BC = ifelse(is.na(sofa_sev_24), NA,
                                ifelse(sofa_sev_24 %in% c("high", "int") & culture == "pos", "high_int_pos",
                                       ifelse(sofa_sev_24 == "low" & culture == "neg", "low_neg", NA)))) %>% 
    mutate(HighInt_Low_BC = factor(HighInt_Low_BC, levels = c("low_neg", "high_int_pos")))
  
  
  return(meta)
}


tr_te_dat$er_tr$meta <- create_outcomes(tr_te_dat$er_tr$meta)
tr_te_dat$er_te$meta <- create_outcomes(tr_te_dat$er_te$meta)
tr_te_dat$icu$meta <- create_outcomes(tr_te_dat$icu$meta)
  
tr_te_dat %>% write_rds("./create_tr_te/tr_te_dat.RDS")
  
  # dplyr::select(one_of("sample_identifier", "condition", "survive", "icu_adm", "culture", "sofa_sev_24", "sofa_sev_prog",
  #                      "HighInt_Low", "High_IntLow", "High_Low", "High_Low_BC")) %>% View()
