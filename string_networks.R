



identify_interactions <- function( my_geneset, string_database_location, species, threshold, mapped_proteins ){
  
  string_db <- STRINGdb::STRINGdb$new(version="12", species=species, score_threshold=0, input_directory=string_database_location, link_data='detailed')
  mapped_proteins <- string_db$map(data.frame(protein = my_geneset), "protein", removeUnmappedRows = TRUE)
  all_interactions_in <- string_db$get_interactions(mapped_proteins$STRING_id)
  interactions <- all_interactions_in %>% distinct(from, to, .keep_all = TRUE)
  
  interactions$score <- 1 - (1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$database/1000)
  
  interactions$alternative_score <- 1 -(1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$neighborhood/1000)*
    (1-interactions$fusion/1000)*(1-interactions$cooccurence/1000)*(1-interactions$textmining/1000)
  
  
  
  return(list(interactions=interactions, string_db = string_db))
}


construct_network <- function( interactions, min_cluster_size, string_db, protein_highlight, score_threshold, mapped_proteins){
  
  interactions <- interactions[interactions$score>=score_threshold, ]
  
  edges <- interactions %>% dplyr::select(from = from, to = to, score = score)
  
  interaction_graph <- tidygraph::tbl_graph(edges = edges, directed = FALSE)
  
  comp <- igraph::components(interaction_graph)
  
  big_comps <- which(comp$csize >= min_cluster_size)
  
  interaction_graph_filtered <- igraph::induced_subgraph(
    interaction_graph,
    vids = igraph::V(interaction_graph)[comp$membership %in% big_comps]
  )
  

  
  
  communities <- igraph::cluster_louvain(interaction_graph_filtered)
  igraph::V(interaction_graph_filtered)$community <- communities$membership
  
  
  
  if(length(protein_highlight)>0){
    highlight <- c()
    
    for(ens_i in igraph::V(interaction_graph_filtered)$name){
      
      prot_i <- mapped_proteins[mapped_proteins$STRING_id == ens_i, ]$protein
      highlight <- c(highlight, prot_i %in% protein_highlight )
      
      
      
    }
    igraph::V(interaction_graph_filtered)$highlight <- highlight
    }
  
  
  return(list(
    graph = interaction_graph_filtered,
    communities = communities, 
    highlights = highlight
  ))
}



prepare_labels <- function(cat_sub, df_ora, communities){
  
  label_out <- list()
  db_out <- c()
  label_genes <- c()
  label_genes_out <- list()
  
  
  n_communities <- length(cat_sub$communities$original_genes)
  
  
  for(i in 1:n_communities){
    label <- c()
    ### Get the label with gene names
    gene_text <- stringr::str_wrap(paste(cat_sub$communities$original_genes[[i]], collapse = " "), width = 40)
    gene_single_text <- stringr::str_wrap(paste(cat_sub$communities$singleton_genes[[i]], collapse = " "), width = 40)
    
    if(length(gene_single_text)>0){
      label_genes <- c(label_genes, paste0("\n\n", gene_text, "\n SINGLE GENES: ", gene_single_text))
    }else{
      label_genes <- c(label_genes, paste0("\n\n", gene_text))
    }
    
    label_genes_out[[i]] <- label_genes
    
    for( db_j in unique(df_ora$database) ){
      
      
      for( com_i in unique(df_ora$category_cluster)){
        
        ind <- (df_ora$category_cluster == com_i) & (df_ora$selected) & (df_ora$database == db_j) & (df_ora$community == i) 
        df_sub <- df_ora[ind, ]
        if(sum(ind)>0){
        
          if(dim(df_sub)[1]>0){
            label_i <- strsplit(df_sub$functional_category, "_")[[1]]
            label_i <- label_i[2:length(label_i)]
          }else{
            label_i <- NA
          }
          
          
          label_i <- tolower(label_i)
          label_i <- stringr::str_wrap(label_i, width = 40)
          
          
          
          db_out <- c(db_out, db_j)
          label <- c(label, paste0(paste0(label_i, collapse =" "), " (", db_j, ")", "\n"))
        }
        
      }
      
    }    
    
    #labels_k <- tolower(labels_k)
    
    label_out[[i]] <- c(paste0(label, collapse=""), "\n\n")  
  } #n_communities 
  return(label_out)
}



get_community_representatives <- function( cat_sub, universe, selection, databases_tested, threshold_mean, threshold_min, ora_min, ora_max ){
  
  n_databases <- length(databases_tested)
  
  for(j in 1:n_databases){
    
    
    
    
    n_communities <- length(cat_sub$communities$all_genes)
    
    for(i in 1:n_communities){
      
      
      df_ora_all <- genes_ORA(cat_sub$communities$all_genes[[i]], universe, db, ora_min, ora_max)
      df_ora <- df_ora_all[df_ora_all$pv_bh<=0.05, ]
      
      if(dim(df_ora)[1]>1){
        overlap_matrix <- get_overlap_matrix( df_ora$functional_category, db )
        df_cluster_ora <- get_sub_clusters_ORA( overlap_matrix, df_ora, threshold_mean, threshold_min )
        df_ora$category_cluster <- df_cluster_ora$cluster
        
        if( selection == "fdr"){
          df_ora <- df_ora %>%
            group_by(category_cluster) %>%
            mutate(
              selected = pv_bh == min(pv_bh, na.rm = TRUE)
            ) %>%
            ungroup()
        }
        
        if( selection == "set_size"){
          df_ora <- df_ora %>%
            group_by(category_cluster) %>%
            mutate(
              selected = pv_bh == min(pv_bh, na.rm = TRUE)
            ) %>%
            ungroup()
        }
        
        df_ora <- as.data.frame(df_ora)
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        if(i == 1){ df_save <- df_ora }else{ df_save <- rbind(df_save, df_ora)}
      }else{
        df_ora$category_cluster <- rep(1, dim(df_ora)[1])
        df_ora$selected <- rep(TRUE, dim(df_ora)[1])
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        
        if(i == 1){ df_save <- df_ora }else{ df_save <- rbind(df_save, df_ora)}
        
      }
      
      
    }
    
    
    
    df_save$database <- rep(strsplit(names(db)[1], '_')[[1]][1], dim(df_save)[1])
    
    if(j == 1){df_ora_combined <- df_save } else{ df_ora_combined <- rbind(df_ora_combined, df_save)}
  } # j 
  
  return(df_ora_combined)
}


identify_category_subclusters <- function(communities, df_singleton, string_db,interaction_graph_filtered, singleton_threshold){
  
  genes_not_in_cluster <- c()
  #genes_in_cluster_singleton <- c()
  #genes_in_cluster_OG <- c()
  #all_genes_in_cluster <- c()
  all_genes <- list()
  original_genes <- list()
  singleton_genes <- list()
  
  for (community_id in unique(communities$membership)) {
    
    #print(paste0("community", community_id))
    # Get the protein list for the current community
    community_proteins <- igraph::V(interaction_graph_filtered)$name[igraph::V(interaction_graph_filtered)$community == community_id]
    
    
    
    df <- string_db$add_proteins_description(data.frame(STRING_id = community_proteins))
    #print(community_id)
    #print(community_id)

      df$preferred_name <- df$preferred_name
    
    
    
 
    
    
    ### add back singletons to the network
    ind_cluster_single <- (df_singleton$cluster_out_max == community_id) & (df_singleton$score_out_max>singleton_threshold)
    ind_not_in_cluster_single <- (df_singleton$cluster_out_max == community_id) & (df_singleton$score_out_max<=singleton_threshold)
    
    genes_not_in_community <- df_singleton[ind_not_in_cluster_single, ]$gene_out 
    genes_single <- df_singleton[ind_cluster_single, ]$gene_out 
    genes_in_community <- df$preferred_name
    
    genes_not_in_cluster <- c(genes_not_in_cluster, genes_not_in_community )
    
    all_genes[[community_id]] <- c(df$preferred_name, genes_single)
    original_genes[[community_id]] <-  df$preferred_name
    
    if(length(genes_single>0)){
      singleton_genes[[community_id]] <-  genes_single
    }else{
      singleton_genes[[community_id]] <-  list(NA)
    }
    
  }
  
  list_out <- list(
    not_community = list(
      genes_not_in_cluster = genes_not_in_cluster
    ),
    communities = list(
      all_genes       = all_genes, 
      original_genes  = original_genes, 
      singleton_genes = singleton_genes
    )
  )
  
  return(list_out)
}



assign_subnetwork_membership <- function(my_geneset, interactions, communities, string_db, mapped_proteins){
  
  
  
  
  all_stringID <- c(interactions$from, interactions$to) 
  all_stringID_unique <- unique(all_stringID)
  
  genes_in_network <- c()
  network_membership <- c()
  
  
  for( i in 1:length(all_stringID_unique)){
    
    si <- communities$names[i]
    
    membership <- communities$membership[i]
    
    genes_in_network <- c(genes_in_network, mapped_proteins[mapped_proteins$STRING_id==si, ]$protein)
    network_membership <- c(network_membership, membership)
  }
  genes_not_in_network <- setdiff(my_geneset, genes_in_network)
  
  names <- c()
  cluster <- c()
  string_name <- c()
  
  ens2membership <- setNames( communities$membership, communities$names )
  
  #### make dataframe output specifying subcluster membership
  for(gene in my_geneset){
    
    ens <- mapped_proteins[mapped_proteins$protein==gene, ]$STRING_id
    
    if(length(ens)==0){
      cluster <- c(cluster, "not_in_cluster")
      string_name <- c(string_name, paste(gene, "(not_in_string)"))
      names <- c(names, gene)
    }else{
      for(ens_i in ens){
        
        if(ens_i %in% communities$names){
          
          
          membership <- ens2membership[[ens_i]]
        
          
          
          
          cluster <- c(cluster, membership)
          names <- c(names, gene)
          string_name <- c(string_name, ens_i)
          
          #print(gene)
        }else if(gene %in% genes_not_in_network){
          
          cluster <- c(cluster, "not_in_cluster")
          names <- c(names, gene)
          string_name <- c(string_name, ens_i)
          
        }else{
          #cluster <- c(cluster, "not_located")
          #print(gene)
        }
        
      }
    } 
    
  }
  
  df_out <- data.frame(cluster, gene = names, string_name )
  df_out <- df_out[order(df_out$cluster), ]
  
  return(list(df_out = df_out, mapped_proteins = mapped_proteins))
  
}


place_singletons <- function( df_out, mapped_proteins, interactions ){
  
  ens_out <- c()
  cluster_out <- c()
  score_out <- c()
  gene_out <- c()
  
  cluster_out_max <- c()
  score_out_max <- c()
  
  cluster_out_mean <- c()
  score_out_mean <- c()
  
  for(ens_i in df_out[df_out$cluster=='not_in_cluster', ]$string_name){
    
    
    ind_from <- which(interactions$from == ens_i)
    ind_to <- which(interactions$to == ens_i)
    
    from_connections <- interactions$to[ind_from]
    from_scores <- interactions$alternative_score[ind_from]
    
    to_connections <- interactions$from[ind_to]
    to_scores <- interactions$alternative_score[ind_to]
    
    
    #print(interactions[which(interactions$from == ens_i | interactions$to == ens_i), ])
    
    clusters <- unique(df_out$cluster)
    n_clusters <- length(clusters)
    
    cluster_score <- setNames(rep(0,n_clusters), clusters)
    cluster_score_max <- setNames(rep(0,n_clusters), clusters)
    cluster_score_mean <- setNames(rep(0,n_clusters), clusters)
    
    cluster_score_from <- setNames(rep(0,n_clusters), clusters) 
    
    #if(F){
    if(length(from_connections)>0){
      for( j in 1:length(from_connections)){
        
        ens_j <- from_connections[j]
        score_j <- from_scores[j]
        
        ind <- ens_j==df_out$string_name
        if(sum(ind)>0){
          
          cluster_i <- df_out$cluster[ind]
          
          cluster_score_from[[cluster_i]]<- cluster_score_from[[cluster_i]] +  score_j
          cluster_score_max[[cluster_i]]<- max(c(cluster_score_max[[cluster_i]], score_j))
          cluster_score_mean[[cluster_i]]<- cluster_score_mean[[cluster_i]] +  score_j/(length(from_connections) + length(to_connections))
          
          #print(paste(ens_j, score_j, df_out$string_name[ind], df_out$cluster[ind]))
          
        }
      }
    }
    #}
    if(length(to_connections)>0){
      cluster_score_to <- setNames(rep(0,n_clusters), clusters) 
      for( j in 1:length(to_connections)){
        ens_j <- to_connections[j]
        score_j <- to_scores[j]
        
        ind <- ens_j==df_out$string_name
        if(sum(ind)>0){
          
          cluster_i <- df_out$cluster[ind]
          n_cluster <- sum(df_out$cluster==cluster_i, na.rm=T)
          cluster_score_to[[cluster_i]]<- cluster_score_to[[cluster_i]] +  score_j
          cluster_score[[cluster_i]]<- cluster_score[[cluster_i]] +  score_j
          cluster_score_max[[cluster_i]]<- max(c(cluster_score_max[[cluster_i]], score_j))
          cluster_score_mean[[cluster_i]]<- cluster_score_mean[[cluster_i]] +  score_j/(length(from_connections) + length(to_connections))
          #print(paste(ens_i, ens_j, score_j, df_out$string_name[ind], df_out$cluster[ind]))
          
        }
      }
    }
    
    
    exclude_name <- "not_in_cluster"  
    filtered_list <- cluster_score[setdiff(names(cluster_score), exclude_name)]
    filtered_list_max <- cluster_score_max[setdiff(names(cluster_score_max), exclude_name)]  
    filtered_list_mean <- cluster_score_mean[setdiff(names(cluster_score_mean), exclude_name)]  
    
    max_name <- names(filtered_list)[which.max(unlist(filtered_list))]
    max_value <- max(unlist(filtered_list))
    
    maxmax_name <- names(filtered_list_max)[which.max(unlist(filtered_list_max))]
    maxmax_value <- max(unlist(filtered_list_max))
    
    meanmax_name <- names(filtered_list_mean)[which.max(unlist(filtered_list_mean))]
    meanmax_value <- max(unlist(filtered_list_mean))
    
    ens_out <- c(ens_out, ens_i)
    cluster_out <- c(cluster_out, max_name)
    score_out <- c(score_out, max_value)
    
    cluster_out_max <- c(cluster_out_max, maxmax_name)
    score_out_max <- c(score_out_max, maxmax_value)
    
    cluster_out_mean <- c(cluster_out_mean, meanmax_name)
    score_out_mean <- c(score_out_mean, meanmax_value)
    
    gn <- mapped_proteins[mapped_proteins$STRING_id==ens_i, ]$protein
    
    gene_out<- c(gene_out, ifelse(length(gn)>0, gn, ens_i))
    
  }
  
  df_singleton <- data.frame(gene_out, cluster_out, score_out, cluster_out_max, score_out_max, cluster_out_mean, score_out_mean)
  
  return(df_singleton)
  
}

save_network_structure <- function(cat_sub, directory, reference){
  
  gene_out <- c()
  com <- c()
  category <- c()
  
  n_communities <- length(cat_sub$communities$original_genes)
  
  original_genes <- cat_sub$communities$original_genes
  singleton_genes <- cat_sub$communities$singleton_genes
  outside_genes <- cat_sub$not_community$genes_not_in_cluster
  
  for(i in 1:n_communities){
    
    community_org_i <- original_genes[[i]]
    community_single_i <- singleton_genes[[i]]
    
    for( gene_i in community_org_i){
      
      gene_out <- c(gene_out, gene_i)
      com <- c(com, i)
      category <- c(category, "original")
      
    } # original genes
    
    for( gene_i in community_single_i){
      
      gene_out <- c(gene_out, gene_i)
      com <- c(com, i)
      category <- c(category, "singleton")
      
    } #singleton genes
  }
  
  
  
  for( gene_i in outside_genes){
    
    gene_out <- c(gene_out, gene_i)
    com <- c(com, NA)
    category <- c(category, "Not_in_network")
    
  } # singleton genes
  
  df_out <- data.frame(gene = gene_out, community = com, network_membership = category)
  
  df_out2 <- df_out[!is.na(df_out$gene), ]
  
  write.xlsx(df_out2, file = paste0(directory, "network_membership_",reference,".xlsx"))
}

