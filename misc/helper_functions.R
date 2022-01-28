### Packages used
library(functionjunction)
library(tidyverse)
library(magrittr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(UpSetR)
library(enrichR)
library(ReactomePA)
library(caret)
library(doParallel)

pathway_hier <- read_csv("/mnt/analysis1/ImportantFiles/Human/enr_pathway_highest_level_clean_names.csv")
pathway_enrichment <- function(gene_list, universe_list, p_val = 0.05, ID = "entrez") {
  
  if(length(gene_list) == 0){
    return(NULL)
  }
  
  gene_list = na.omit(gene_list)
  universe_list = na.omit(universe_list)
  
  if(ID == "ensembl") {
    gene_list <- data.frame(ensembl_gene_id = gene_list) %>% left_join(universe) %>% pull("entrezgene_id")
    
  }
  
  # Set res to NULL
  res <- NULL
  
  # Get results
  tryCatch({
    res <- ReactomePA::enrichPathway(gene = as.character(gene_list), universe = as.character(universe_list))@result
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")} )
  
  if (is.null(res)){
    return(NULL)
  } else {
    res <- res %>% dplyr::mutate(geneID = stringr::str_split(geneID, pattern = "/")) %>% dplyr::filter(p.adjust <= p_val)
    return(res)
  }
}

go_enrichment <- function(gene_list, p_val = 0.05, ID = "hgnc_symbol"){
  
  # Check if the gene list is empty
  if(length(gene_list) == 0){
    return(NULL)
  }
  
  gene_list = na.omit(gene_list)
  
  # Set res to NULL
  res <- NULL
  
  # Get results
  tryCatch({
    res <- enrichR::enrichr(genes = as.character(gene_list), 
                            databases = c("MSigDB_Hallmark_2020", 
                                          "GO_Molecular_Function_2018", 
                                          "GO_Cellular_Component_2018", 
                                          "GO_Biological_Process_2018") 
    )
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")} )
  
  if (is.null(res)){
    return(NULL)
  } else {
    res <- res %>% map(~dplyr::filter(.x, Adjusted.P.value <= p_val))
    return(res)
  }
}

perform_pca <- function(expr){
  df <- expr %>% 
    as.data.frame() %>% 
    #rownames_to_column(var = "sample_identifier") %>% 
    #dplyr::select(one_of("sample_identifier", features)) %>% 
    #column_to_rownames(var = "sample_identifier") %>% 
    prcomp(center = TRUE, scale. = TRUE)
  
  PoV <- round(df$sdev^2/sum(df$sdev^2), digits = 3)
  x <- df$x %>% as.data.frame() %>% rownames_to_column(var = "sample_identifier") 
  rot <- df$rotation %>% as.data.frame() %>% rownames_to_column(var = "ensembl_gene_id")
  
  res <- list(x = x, rot = rot,  pov = PoV)
  cat(paste0("PCA ---- DONE.\n"))
  return(res)
}

explain_95 <- function(ev, perc){
  for (ind in 1:length(ev)){
    var_expl <- sum(ev[1:ind])
    if(var_expl > perc/100){
      return(ind)
    }
  }
}

de <- function(counts, meta,  PADJ= 0.05, FC= 1.5, des, main_covar, versus_hc = FALSE, hc_name = "healthy_control", filt = TRUE) {
  
  # Load required package 
  #require("DESeq2")
  
  # Function to clean up DE results
  clean_de_res <- function(res){
    res_clean <- res %>% 
      as.data.frame() %>% 
      rownames_to_column(var= "ensembl_gene_id") %>% 
      left_join(universe, by = "ensembl_gene_id") %>% 
      mutate(Abs_LFC = abs(log2FoldChange)) %>% 
      mutate(fc= sign(log2FoldChange)* 2^(Abs_LFC)) %>% 
      mutate(padj = ifelse(is.na(padj), 1, padj )) %>% 
      mutate(de = ifelse(padj <= PADJ & (fc >= FC | fc <= -FC), "de", "non_de")) %>% 
      mutate(direction = ifelse( fc >= 0 & de == "de", "up" ,ifelse(fc <= 0 & de == "de", "down", "non_de")))
    
    if(filt == TRUE) {
      res_clean_filt <- res_clean %>%  dplyr::filter(de == "de")
      return(res_clean_filt)
    }
    return(res_clean)
  }
  
  # Set up function for categorical main covariate
  de_factor <- function() {
    
    comparisons <- combn(levels(colData(dds)[[main_covar]]), 2, simplify = FALSE)
    
    if (isTRUE(versus_hc)) {
      comparisons %<>% keep(~any(c(hc_name) %in% .x))
    }
    
    FinalResults <- list()
    for (i in 1:length(comparisons)) {
      res <- DESeq2::results(dds, contrast = c(main_covar, comparisons[[i]][2], comparisons[[i]][1]))
      res %<>% clean_de_res()
      
      name <- paste0(comparisons[[i]][2],"_vs_", comparisons[[i]][1])
      FinalResults[[name]] <- res
    }
    return(FinalResults)
    
  }
  
  # Set up function for numerical main covariate
  de_numeric <- function() {
    FinalResults <- list()
    res <- DESeq2::results(dds, contrast = list(main_covar))
    res %<>% clean_de_res()
    FinalResults[[main_covar]] <- res
    return(FinalResults)
  }
  
  # Set up dds
  cat("Creating dds object...\n")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta, design = as.formula(des))
  
  # What is the class of the main covariate? 
  main_covar_class <- class(SummarizedExperiment::colData(dds)[[main_covar]])
  
  # Perform DE based on main covariate class 
  if (main_covar_class == "numeric") {
    cat("\n The main covariate is a numeric. \n")
    # Perform LMs 
    dds <- DESeq2::DESeq(dds, parallel = TRUE)
    return(de_numeric())
  } else if (main_covar_class == "factor") {
    # Perform LMs 
    dds <- DESeq2::DESeq(dds, parallel = TRUE)
    cat("\n Main covariate is a factor. \n")
    return(de_factor())
  } else {
    cat("\n Please make variable a factor or numeric! \n")
    return(NULL)
  }
}


de_gene_numbers <- function(de_res, comparison_name = NULL){
  return(data.frame(
    all = de_res %>% filter(de == "de") %>% nrow(),
    up = de_res %>% filter(de == "de" & direction == "up") %>% nrow(),
    down = de_res %>% filter(de == "de" & direction == "down") %>% nrow(),
    row.names = comparison_name)
  )
  }

plot_reactome_pathway <- function(df){
  pathway_hier <- read_csv("/mnt/analysis1/ImportantFiles/Human/enr_pathway_highest_level_clean_names.csv")
  keep_class <- c("Adaptive", "Innate", "Cytokine Signaling", "Signaling by GPCR" , 
                  "Platelet activation, signaling and aggregation") 
  remove_class <- c("Cell Cycle, Mitotic", NA, "Cell Cycle Checkpoints", "O2/CO2 exchange in erythrocytes")
  
  df <- df %>% 
    separate(BgRatio,into = c("M", "N")) %>%
    mutate(Ratio = Count/as.numeric(M) ) %>% 
    left_join(pathway_hier, by = c("ID" = "enr_pathway")) %>% 
    dplyr::select(one_of("comparison","direction", "p.adjust", 
                         "enr_pathway_descrip_clean_name","Ratio", "one_lower_level_pathway_descrip_clean_name",
                         "geneID")) %>% 
    filter(direction %in% c("down", "up")) %>% 
    mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down"))) %>% 
    #filter(one_lower_level_pathway_descrip_clean_name %in% keep_class) %>% 
    filter(!one_lower_level_pathway_descrip_clean_name %in% remove_class) %>% 
    mutate(one_lower_level_pathway_descrip_clean_name = str_wrap(one_lower_level_pathway_descrip_clean_name, width = 20)) %>% 
    mutate(enr_pathway_descrip_clean_name = str_wrap(enr_pathway_descrip_clean_name, width = 40))
  
  res <- ggplot(df, aes(y = enr_pathway_descrip_clean_name, x = direction, fill = Ratio, size = -log10(p.adjust))) +
    facet_grid(rows = vars(one_lower_level_pathway_descrip_clean_name), 
               cols = vars(comparison), scales = "free_y", space = "free_y" ) +
    geom_point(color = "black", pch = 22) +
    theme_light() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(color = "black", size = 12),
          legend.position = "right",
          strip.text.y = element_text(angle =0, color = "black" ),
          strip.text = element_text(color = "black", size = 12)
    ) +
    scale_fill_viridis(direction = -1) +
    scale_size(name = "Adj p val \n (-log10)") +
    ylab("") + xlab("")
  return(res)
}

plot_msigdb_pathways <- function(df){
  df <- df %>% 
    separate(Overlap,into = c("M", "N")) %>%
    mutate(Ratio = as.numeric(M)/as.numeric(N) ) %>% 
    filter(direction %in% c("down", "up")) %>% 
    mutate(direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down")))
  
  res <- ggplot(df, aes(x = direction, y = Term, fill = Ratio, size = -log10(Adjusted.P.value))) + 
    geom_point() + 
    facet_grid(cols = vars(comparison), scales = "free_y", space = "free_y") +
    geom_point(color = "black", pch = 22) +
    theme_light() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(color = "black", size = 12),
          legend.position = "right",
          strip.text.y = element_text(angle = 0, color = "black" ),
          strip.text = element_text(color = "black", size = 12)
    ) +
    scale_fill_viridis(direction = -1) +
    scale_size(name = "Adj p val \n (-log10)") +
    ylab("") + xlab("")
  return(res)
}

vp <- function(dat, title, colour, fc_cutoff = 1, pval = 0.05, labs_to_incl = 5){
  
  top_pos <- dat %>% filter(de == "de") %>% top_n(labs_to_incl, fc)
  top_neg <- dat %>% filter(de == "de") %>% top_n(labs_to_incl, fc*-1)
  
  x_axis_neg = ceiling(min(top_neg$log2FoldChange)) ; x_axis_pos = ceiling(max(top_pos$log2FoldChange))
  y_axis = ceiling(-log10(min(c(top_pos$padj, top_neg$padj)))) 
  df <- dat %>% dplyr::rename(Direction = "direction") %>% 
    mutate(Direction = ifelse(Direction == "down", "Down", ifelse(Direction == "up", "Up", "Non DE"))) %>% 
    mutate(label_to_incl = ifelse(hgnc_symbol %in% c(top_pos$hgnc_symbol, top_neg$hgnc_symbol), hgnc_symbol, NA))
  
  fig <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj) )) + 
    geom_point(aes(color = Direction), alpha = 0.3) + 
    geom_text_repel(aes(label = label_to_incl), fontface = "bold") +
    scale_color_manual(values = c(colour, "grey" ,colour)) + 
    geom_vline(xintercept = c(-log2(fc_cutoff) , log2(fc_cutoff) ), linetype = "dotted") +
    geom_hline(yintercept = c( -log10(pval) ), linetype = "dotted") +
    theme_minimal(base_size = 15) %+replace%
    theme(legend.position = "none", axis.text.x = element_text(color = "black", size = 13), 
          axis.text.y = element_text(color = "black", size = 13)) +
    ylab("") + xlab(title) +
    #geom_text(x = 2, y = y_axis -1, label = "Up", size = 5) +
    #geom_text(x = -2.5, y = y_axis - 1, label = "Down", size = 5 ) +
    xlim((x_axis_pos)*-1 - 1, x_axis_pos + 1) + ylim(0, y_axis) 
  return(fig)
}

extract_coefs <- function(mod) {
  
  response_levs = length(mod$levels)
  
  # Multi-class 
  if (response_levs > 2) {
    coefs <-  coef(mod$finalModel, mod$bestTune$lambda) %>% 
      map(~as.matrix(.x) %>% as.data.frame() %>% rownames_to_column(var = "ensembl_gene_id")) %>% 
      bind_rows(.id = "classes") %>% 
      dplyr::rename(coef = "1") %>% 
      filter(!ensembl_gene_id == "(Intercept)") %>% 
      filter(coef != 0) %>% 
      left_join(universe, by = "ensembl_gene_id")
    return(coefs)
  } else{
    coefs <-  coef(mod$finalModel, mod$bestTune$lambda) %>% 
      as.matrix() %>% as.data.frame() %>% rownames_to_column(var = "ensembl_gene_id") %>%
      dplyr::rename(coef = "1") %>%
      filter(!ensembl_gene_id == "(Intercept)") %>%
      filter(coef != 0) %>%
      left_join(universe, by = "ensembl_gene_id")
  }
  return(coefs)    
}

#### Misc
universe <- read_rds("../../sepsis_rnaseq_all/final/counts_meta_10_read_filt_261120.RDS")$universe

# Get Jaccard Index Matrix
jaccard <- function(M, user1, user2) {
  sums = rowSums(M[,c(user1, user2)])
  
  similarity = length(sums[sums==2])
  total = length(sums[sums==1]) + similarity
  
  similarity/total
}

get_jac_mat <- function(list, backgroud){
  mat <- list %>% 
    map(~dplyr::select(.x, one_of("ensembl_gene_id"))) %>% 
    bind_rows(.id = "comp") %>% 
    mutate(present = 1) %>% 
    spread(comp, present) %>% 
    replace(is.na(.),0) %>% 
    column_to_rownames(var = "ensembl_gene_id")
  return(mat)
}
