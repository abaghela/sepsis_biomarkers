rm(list = ls())
source("./helper_functions.R")

# Read in data 
tr_te_dat <- read_rds("./create_tr_te/tr_te_dat.RDS")
meta <- tr_te_dat$er_tr$meta %>% filter(condition != "healthy_control")

# Which clinical variables are we interested in 
meta_vars_of_int <- meta %>% 
  mutate(outcome_icu_admission = ifelse(outcome_icu_admission == 1, "1", "0")) %>% 
  mutate(outcome_icu_admission = factor(outcome_icu_admission, levels = c("0", "1"))) %>% 
  mutate(micro_blood_culture_pathogen = ifelse(micro_blood_culture_pathogen == 1, "1", "0")) %>% 
  mutate(micro_blood_culture_pathogen = factor(micro_blood_culture_pathogen, levels = c("0", "1"))) %>% 
  dplyr::select(one_of("sample_identifier", "first_at_ed_sofa",  "worst_within_72_sofa", "outcome_icu_admission", 
                       "micro_blood_culture_pathogen")) 

# Read in Module eigenvectors
ICA_eigenvec <- read_rds("./modules/results/ICA_mod_eigen.rds")

# Perform LMs
lm_clin_eig <- function(ICA_eig) {
  
  # Check things are in the right order 
   right_order <- ICA_eig %>% 
    map(~all(.x$sample_identifier == meta_vars_of_int$sample_identifier)) %>% 
    unlist() %>% 
    all()
  
  if (isFALSE(right_order)) {
    return(c("Things are not in the right order"))
    
  }
   
   # What variables do we have?
   lm_vars <- colnames(meta_vars_of_int[,-1])
   
   # Initiate loop - loop through variables, and then module
   res <- list()
   for (var in lm_vars) {
    for (mod in names(ICA_eig)) {
      
      # Create a new df with the var and eigenvector the module
      lm_df <- ICA_eig[[mod]] %>% 
        left_join(meta_vars_of_int, by = "sample_identifier") %>% 
        dplyr::select(one_of("sample_identifier","PC1", var)) %>% 
        na.omit()
      
      # What is the class of the variable?
      var_class <- class(meta_vars_of_int[[var]])
    
      if(var_class %in% c("numeric")) {
        
        # GLM - var ~ PC1
        res[[var]][[mod]] <- glm(formula = formula(paste0(var, " ~ PC1")), family = gaussian, data = lm_df ) %>% 
          broom::tidy()
        } else {
        res[[var]][[mod]] <- glm(formula = formula(paste0(var, " ~ PC1")), family = binomial, data = lm_df ) %>% 
          broom::tidy()
        }
      }
     }
   return(res)
   } 


clinical_assoc <- ICA_eigenvec %>% map(~lm_clin_eig(.x))

plot_signif_vars <- function(clin_assoc, ICA_eig){
  mod_clin_pval <- clin_assoc %>%
    map(~bind_rows(.x, .id = "comp")) %>%
    map(~filter(.x, term == "PC1")) %>%
    map(~mutate(.x, adj.p.value = p.adjust(p.value, "BH"))) %>%
    map(~filter(.x, adj.p.value <= 0.05)) %>%
    map(~arrange(.x, adj.p.value)) %>% 
    discard(~nrow(.x) == 0)
  
  mod_clin_pval_top <- mod_clin_pval %>% 
    map(~top_n(.x, -5, adj.p.value))
  
  # Plot numeric factors 
  
  res <- list()
  res[["full_res"]] <- mod_clin_pval
  for (var in names(mod_clin_pval_top)) {
    
    # What is the class of the variable?
    var_class <- class(meta_vars_of_int[[var]])
    
    if(var_class %in% c("numeric")) {
      
      res[["plot"]][[var]] <- ICA_eig %>%
        bind_rows(.id = "comp") %>%
        filter(comp %in% mod_clin_pval_top[[var]]$comp) %>%
        left_join(meta_vars_of_int) %>% 
        dplyr::select(one_of("comp", "sample_identifier", "PC1", var)) %>% 
        pivot_longer(cols = -c("comp", "sample_identifier", "PC1"), names_to = "variable") %>%
        ggplot(aes(PC1, value, color = value)) + geom_point(position = "jitter") + 
        facet_grid(rows = vars(variable), cols = vars(comp), scales = "free_y") + 
        geom_smooth(method = "lm") +
        scale_color_viridis(option = "magma", direction = -1) +
        theme_minimal()
    
    } else if (var_class %in% c("factor")) {
      res[["plot"]][[var]] <- NULL
      # df_plot <- ICA_eigenvec$expr_vst %>%
      #   bind_rows(.id = "comp") %>%
      #   filter(comp %in% mod_clin_pval[[var]]$comp) %>%
      #   left_join(meta_vars_of_int) %>% 
      #   dplyr::select(one_of("comp", "sample_identifier", "PC1", var)) %>% 
      #   pivot_longer(cols = -c("comp", "sample_identifier", "PC1"), names_to = "variable") %>%
      #   filter(!is.na(value)) %>% 
      #   mutate_if(is.factor, ~ifelse(.x =="1", 1, 0))
      # 
      # h = df_plot %>% group_by(value) %>% 
      #   mutate(breaks = cut(PC1, breaks = seq(-30, 30, 1), labels = seq(-29.75, 30, 1), include.lowest = TRUE ),
      #          breaks = as.numeric(as.character(breaks))) %>% 
      #   group_by(value, breaks) %>% 
      #   summarize(n = n()) %>% 
      #   mutate(pct = ifelse(value == 0, n/sum(n), 1 - n/sum(n)))
      # 
      # 
      # ggplot() + 
      #   geom_segment(data=h, size=4, show.legend=FALSE,
      #                aes(x=breaks, xend=breaks, y=value, yend=pct, colour=factor(value))) +
      #   geom_segment(dat=df_plot[df_plot$value==0,], 
      #                aes(x=PC1, xend=PC1, y=0, yend=-0.02), size=0.2, colour="grey30") +
      #   geom_segment(dat=df_plot[df_plot$value==1,], 
      #                aes(x=PC1, xend=PC1, y=1, yend=1.02), size=0.2, colour="grey30") + 
      #   geom_line(data=data.frame(x=seq(-30, 30, 1), 
      #                             y=predict(glm(value ~ PC1, family="binomial", data=df_plot), 
      #                                       newdata=data.frame(PC1=seq(-30, 30, 1)),
      #                                       type="response")), 
      #             aes(x,y), colour="grey50", lwd=1) +
      #   scale_y_continuous(limits=c(-0.02,1.02)) +
      #   scale_x_continuous(limits=c(-30,30)) +
      #   facet_grid()
      #   theme_bw(base_size=12)
      # 
      #   
      # 
      # temp_df_plot %>% 
      #   ggplot(aes(PC1, value)) + geom_point() + 
      #   facet_grid(rows = vars(variable), cols = vars(comp), scales = "free_y") + 
      #   geom_smooth(aes(PC1, value), method="glm", method.args=list(family="binomial"), se=FALSE) +
      #   #geom_histogram(data=temp_df_plot[temp_df_plot$value==0,], aes(x=PC1)) +
      #   #geom_histogram(data=temp_df_plot[temp_df_plot$value==1,], aes(x=PC1))
      #   #scale_color_viridis(option = "magma", direction = -1) +
      #   theme_minimal()
      
    }
  }
  return(res)
}

clinical_assoc_filt <- clinical_assoc %>% map2(ICA_eigenvec, ~plot_signif_vars(.x,.y ))

clinical_assoc %>% write_rds("./modules/results/ICA_mod_eigen_clin_assoc.rds")
clinical_assoc_filt %>% write_rds("./modules/results/ICA_mod_eigen_clin_assoc_filt_plot.rds")
