

##https://www.gsea-msigdb.org/gsea/msigdb/mouse/genesets.jsp

#library(Rmpfr)
#library(fgsea)


get_database <- function( database, folder ){
  
  if(database == 'gobp'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/msigdb_v2024.1.Mm_GMTs/m5.go.bp.v2024.1.Mm.symbols.gmt"))
  }
  
  if(database == 'gocc'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/msigdb_v2024.1.Mm_GMTs/m5.go.cc.v2024.1.Mm.symbols.gmt"))
  }
  
  if(database == 'gomf'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/msigdb_v2024.1.Mm_GMTs/m5.go.mf.v2024.1.Mm.symbols.gmt"))
  }
  
  if(database == 'reactome'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/msigdb_v2024.1.Mm_GMTs/m2.cp.reactome.v2024.1.Mm.symbols.gmt"))
  }
  
  if(database == 'reactome_human'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
  }
  
  if(database == 'gobp_human'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/c5.go.bp.v2025.1.Hs.symbols.gmt"))
  }
  
  if(database == 'gocc_human'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/c5.go.cc.v2025.1.Hs.symbols.gmt"))
  }
  
  return(pathways)
}


genes_ORA <- function(genes_drawn, universe, database, min_size, max_size){
  
  genes_drawn <- c(na.omit(genes_drawn))
  universe <- c(na.omit(universe))
  
  set_size <- c()
  set_size_in_universe <- c()
  
  pv <- c()
  pv_calc <- c()
  pw <- c()
  gene_ratio <- c()
  geneset <- c()
  geneset_drawn <- c()
  gene_count <- c()
  geneset_drawn_in_universe <- c()
  
  geneset_tmp <- c()
  for(pw_i in names(database)){
    
    
    
    
    

      genes_i <- unname(database[pw_i])[[1]]
    
    
    geneset_i <- paste(intersect(genes_i, universe), collapse="/")
    
    geneset_tmp  <- c(geneset_tmp , geneset_i)
    
  }
  
  unique_indices <- which(!duplicated(geneset_tmp ))
  

  for(pw_i in names(database[unique_indices])){
    
    

      genes_i <- unname(database[pw_i])[[1]]
    
    
    N <- length(universe)
    n <- length(intersect(genes_drawn, universe))
    K <- length(intersect(genes_i, universe))
    k <- length(intersect(intersect(genes_i, genes_drawn), universe))
    
    if((K>=min_size) & (K<=max_size)){
      
      geneset <- c(geneset, paste(intersect(genes_i, universe), collapse="/"))
      geneset_drawn <- c(geneset_drawn, paste(intersect(genes_i, genes_drawn), collapse="/"))
      geneset_drawn_in_universe <- c(geneset_drawn_in_universe, paste(intersect(intersect(genes_i, genes_drawn), universe), collapse="/"))
      pw <- c(pw, pw_i)
      pv <- c(pv, Rmpfr::asNumeric(1-exp(Rmpfr::mpfr(phyper(k-1, K, N-K, n, lower.tail = TRUE, log.p = TRUE), 2000))))
      gene_ratio <- c(gene_ratio, as.double(k)/as.double(K))
      gene_count <- c(gene_count, k)
      set_size <- c(set_size, length(genes_i))
      set_size_in_universe <- c(set_size_in_universe, K)
      
    }  
  }

  df_out <- data.frame(functional_category = pw, pv, gene_ratio, geneset_drawn,geneset_drawn_in_universe, geneset, gene_count, set_size, set_size_in_universe)

  df_out$pv_bh <- p.adjust(df_out$pv, method="BH")
  
  return(df_out[order(log(df_out$pv)), ])
}





