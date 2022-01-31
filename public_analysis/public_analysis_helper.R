# PCA 
perform_pca <- function(expr_mat){
    # Perform PCA with prcomp
  df <- expr_mat %>% 
    as.data.frame() %>% 
    prcomp(center = TRUE, scale. = TRUE)
  
  PoV <- round(df$sdev^2/sum(df$sdev^2), digits = 3)
  x <- df$x %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "sample_identifier") 
  rot <- df$rotation %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "ensembl_gene_id")
  
  res <- list(x = x, rot = rot,  pov = PoV)
  cat(paste0("PCA ---- DONE.\n"))
  return(res)
}

# Quick helper function to extract the number of PCs explaining X percent of variation in the data 
explain_95 <- function(ev, perc){
  for (ind in 1:length(ev)){
    var_expl <- sum(ev[1:ind])
    if(var_expl > perc/100){
      return(ind)
    }
  }
}

quiet <- function(..., messages=FALSE, cat=FALSE){
    if(!cat){
        sink(tempfile())
        on.exit(sink())
    }
    out <- if(messages) eval(...) else suppressMessages(eval(...))
    out
}