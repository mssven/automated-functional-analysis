


check_cluster_correlation <- function(indices, M, threshold_mean, threshold_min) {
  
  if (length(indices) < 2) return(TRUE) # clusters of size 1 are fine
    subM <- M[indices, indices]
    avg_cor <- mean(subM[upper.tri(subM)], na.rm = TRUE)
    min_cor <- min(subM[upper.tri(subM)], na.rm = TRUE)

  return(avg_cor >= threshold_mean & min_cor>=threshold_min)

  }


get_overlap_matrix <- function( pws, database ){
  

  overlap_matrix <- matrix(0, nrow=length(pws), ncol=length(pws), 
                           dimnames=list(names(pws), names(pws)))
  
  for( i in 1:length(pws)){
    
 
      list_of_genes_1 <- database[[pws[i]]]

    go_term_1 <- pws[i]

    
    for( j in 1:length(pws)){
      

        list_of_genes_2 <- database[[pws[j]]]
      
        
        go_term_2 <- pws[j]
      overlap_matrix[i, j] <- length(intersect(list_of_genes_1, list_of_genes_2))/(min(c(length(list_of_genes_1),  length(list_of_genes_2))))
      
    }
  }
  
  
  M <- as.data.frame(overlap_matrix)
  
  rownames(M) <- pws
  colnames(M) <- pws
  
  
  return(M)
}


get_sub_clusters_ORA <- function( M, df_fa_in, threshold_mean, threshold_min ){
  
  
  genes_in_set <- df_fa_in$gene_count
  set_size <- df_fa_in$set_size
  set_size_in_universe <- df_fa_in$set_size_in_universe
  fdr <- df_fa_in$pv_bh
  pws <- df_fa_in$functional_category
  hc <- hclust(d = as.dist(1-M), method = "ward.D2")
  dend <- as.dendrogram(hc)
  leaves <- labels(dend)
  idx <- match(leaves, rownames(M))
  
  count <- 1
  clusters <- list()
  
  dend <- as.dendrogram(hc)
  not_yet_clusters <- list(dend)
  
  
  
  while (length(not_yet_clusters) > 0) {
    
    # Take the first element from the queue
    
    nyc <- not_yet_clusters[[1]]
    not_yet_clusters <- not_yet_clusters[-1]
    
    # Get indices of leaves under this branch
    leaves <- labels(nyc)
    idx <- match(leaves, rownames(M))
    
    
    # Check correlation or stop if leaf
    if (check_cluster_correlation(idx, M, threshold_mean , threshold_min ) || is.leaf(nyc)) {
      clusters[[count]] <- idx   # store indices (or nyc if you want sub-dendrograms)
      count <- count + 1
    } else {
      # Add children back to queue for further checking
      not_yet_clusters <- c(not_yet_clusters, list(nyc[[1]], nyc[[2]]))
    }
  }
  
  cluster_labels <- rep(NA, nrow(M))
  
  for (i in seq_along(clusters)) {
    cluster_labels[clusters[[i]]] <- i
  }
  
  
  
  ord <- hc$order  # keep row order consistent with heatmap
  
  
  
  term_clusters <- data.frame(
    pw        = pws[ord],
    cluster   = cluster_labels[ord],
    fdr       = fdr[ord],
    set_size  = set_size[ord], 
    set_size_in_universe  = set_size_in_universe[ord], 
    genes_in_set = genes_in_set[ord]
  )

  
  return(term_clusters)
}



get_sub_GSEA <- function( M, fgseaRes_in, threshold_mean, threshold_min ){
  
  genes_in_set <- unlist(lapply(fgseaRes_in$leadingEdge, length))
  NES <- fgseaRes_in$NES
  set_size <- fgseaRes_in$size
  fdr <- fgseaRes_in$padj
  pws <- fgseaRes_in$pathway
  
  ph <- pheatmap::pheatmap(
    M,
    clustering_method = "ward.D2", 
    silent = TRUE
  )
  
  # Extract row order
  idx <- ph$tree_row$order
  
  count <- 1
  clusters <- list()
  
  dend <- as.dendrogram(ph$tree_row)
  not_yet_clusters <- list(dend)
  
  par(mar = c(2, 2, 2, 30))  # bottom margin bigger
  plot(dend, horiz = TRUE, cex = 0.2)   
  
  
  while (length(not_yet_clusters) > 0) {
    
    # Take the first element from the queue
    
    nyc <- not_yet_clusters[[1]]
    not_yet_clusters <- not_yet_clusters[-1]
    
    # Get indices of leaves under this branch
    leaves <- labels(nyc)
    idx <- match(leaves, rownames(M))
    
    
    # Check correlation or stop if leaf
    if (check_cluster_correlation(idx, M, threshold_mean , threshold_min ) || is.leaf(nyc)) {
      clusters[[count]] <- idx   # store indices (or nyc if you want sub-dendrograms)
      count <- count + 1
    } else {
      # Add children back to queue for further checking
      not_yet_clusters <- c(not_yet_clusters, list(nyc[[1]], nyc[[2]]))
    }
  }
  
  cluster_labels <- rep(NA, nrow(M))
  
  for (i in seq_along(clusters)) {
    cluster_labels[clusters[[i]]] <- i
  }
  
  
  
  ord <- ph$tree_row$order  # keep row order consistent with heatmap
  
  
  
  term_clusters <- data.frame(
    pw        = pws[ord],
    cluster   = cluster_labels[ord],
    fdr       = fdr[ord],
    set_size  = set_size[ord],
    NES  = NES[ord],
    genes_in_set = genes_in_set[ord]
  )
  
  
  
  return(term_clusters)
}

