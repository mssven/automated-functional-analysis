


library(dplyr)



automated_functional_analysis <- function( my_geneset, 
                               universe, 
                               threshold_mean, 
                               threshold_min, 
                               databases_tested, 
                               min_cluster, 
                               species, 
                               ora_min, 
                               ora_max, 
                               directory, 
                               reference, 
                               string_score_threshold,
                               singleton_threshold,
                               network_height, 
                               network_width,
                               protein_highlight = c()  ){
  
  source("D:/Function_for_functional_analysis/R/functional_analysis_functions.R")
  source("D:/Function_for_functional_analysis/R/functions_term_reduction.R")
  source("D:/Function_for_functional_analysis/R/string_networks.R")
 
  
  string_database_location <- "D:/Function_for_functional_analysis/String"
  

  out_ii <- identify_interactions( my_geneset, string_database_location, species)  
  

  interactions <- out_ii$interactions
  string_db <- out_ii$string_db
  
  mapped_proteins <- string_db$map(data.frame(protein = my_geneset), "protein", removeUnmappedRows = TRUE)

  out_cn <- construct_network( interactions, min_cluster, string_db, protein_highlight, string_score_threshold, mapped_proteins )

  communities <- out_cn$communities
  interaction_graph_filtered <- out_cn$graph
  
  #mapped_proteins <- string_db$map(data.frame(protein = my_geneset), "protein", removeUnmappedRows = TRUE)
  
  # Find singletons 
  all_stringID <- c(interactions$from, interactions$to) 
  all_stringID_unique <- unique(all_stringID)

  genes_in_network <- c()
  network_membership <- c()
  
  
  out_asm <- assign_subnetwork_membership(my_geneset, interactions, communities, string_db, mapped_proteins)
  df_out <- out_asm$df_out
  
 
  df_singleton <- place_singletons( df_out, mapped_proteins, interactions )

  layout_fr  <- ggraph::create_layout(interaction_graph_filtered, layout = 'fr')
  
  selection <- "fdr"


  cat_sub <- identify_category_subclusters(communities, df_singleton, string_db, interaction_graph_filtered,  singleton_threshold)
  df_ora <- get_community_representatives( cat_sub, universe, selection, databases_tested, threshold_mean, threshold_min, ora_min, ora_max )
  labels <- prepare_labels(cat_sub, df_ora, communities)
  
  write.csv(df_ora, file = paste0(directory, "/overrepresentation_analysis_", reference, ".csv"))
  
  save_network_structure(cat_sub, directory, reference)
  
  gene_name <- c()
  for(str_i in layout_fr$name){
    gene_name <- c(gene_name, mapped_proteins$protein[mapped_proteins$STRING_id == str_i])
  }
  
  layout_fr$gene_name <- gene_name
  
  
  if(length(protein_highlight)>0){
    layout_fr$highlight <- igraph::V(interaction_graph_filtered)$highlight
  }
  
  
  
  g02 <- ggraph(layout_fr)  + 
    geom_edge_link(color = "grey40", width = 0.5, alpha = 0.7) +       # Draw the edges
    geom_node_point(
      aes(fill = as.factor(community)), 
      size = 5, 
      shape = 21, 
      stroke = 0.7, 
      color = "black"
    )+
    geom_node_text(
      aes(label = gene_name),
      repel = TRUE,        # avoids text overlapping
      size = 5
    ) + 
    theme_bw() + 
    theme(
      legend.text = ggtext::element_markdown(),  # <- this enables HTML/Markdown rendering
      axis.text.x = ggtext::element_markdown(), 
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(), 
      legend.position = "none"
    ) 
  
  if(length(protein_highlight)>0){
    g_figure <- ggraph::ggraph(layout_fr) +
      ggraph::geom_edge_link(edge_colour = "black") +
      
      # Outer ring: LMB1
      
      ggraph::geom_node_point(
        aes(fill  = as.factor(community)),
        color = "black",
        size  = 4,
        shape = 21,
        stroke = ifelse(layout_fr$highlight, 3, 0),
        #show.legend = FALSE
      ) + 
      ggplot2::theme_bw()
    
    g_figure <- g_figure + 
      ggplot2::scale_fill_discrete(
        name   = "",
        labels = labels
      )
    
    
  }else{
  g_figure <- ggraph::ggraph(layout_fr)  + 
    ggraph::geom_edge_link(edge_colour = 'black') +
    ggraph::geom_node_point(ggplot2::aes(color = as.factor(community)), size = 5) +
    ggplot2::theme_bw()
  
  g_figure <- g_figure + 
    ggplot2::scale_color_discrete(
      name   = "",
      labels = labels
    )
  }
  


  
  
  g_figure2 <- g_figure + 
    geom_node_text(
      aes(label = gene_name),
      repel = TRUE,        # avoids text overlapping
      size = 5
    )# + 
    #theme_bw() + 
    #theme(
    #  legend.text = ggtext::element_markdown(),  # <- this enables HTML/Markdown rendering
    #  axis.text.x = ggtext::element_markdown(), 
    #  panel.background = element_blank(),
    #  panel.grid = element_blank(),
    #  axis.ticks = element_blank(),
    #  axis.text = element_blank(),
    #  axis.title = element_blank()#, 
      #legend.position = "none"
    #) 
  
  ggsave( plot = g_figure, filename = paste0(directory, "/ppi_network_", reference, ".png"), width = network_width, height = network_height )
  ggsave( plot = g_figure2, filename = paste0(directory, "/ppi_network_with_labels_", reference, ".png"), width = network_width, height = network_height )
  
  return(list(figure = g_figure, results = df_ora, network_structure = cat_sub, network_layout = layout_fr, labels = labels, mapped_proteins))
  
}

  
  
  
  
replot_network <- function(layout_fr, directory, network_height, network_width, label_size, legend_size, labels){
  
highlight <- layout_fr$highlight
community <- layout_fr$community

if(length(highlight)>0){
  
  highlight <- layout_fr$highlight
  g_figure <- ggraph::ggraph(layout_fr) +
    ggraph::geom_edge_link(edge_colour = "black") +
    
    # Outer ring: LMB1
    
    ggraph::geom_node_point(
      aes(fill  = as.factor(community)),
      color = "black",
      size  = 4,
      shape = 21,
      stroke = ifelse(layout_fr$highlight, 3, 0),
      #show.legend = FALSE
    ) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      legend.text  = ggplot2::element_text(size = legend_size),
      legend.title = ggplot2::element_text(size = legend_size)
    )
  
  g_figure <- g_figure + 
    ggplot2::scale_fill_discrete(
      name   = "",
      labels = labels
    )
  
  
}else{
  g_figure <- ggraph::ggraph(layout_fr)  + 
    ggraph::geom_edge_link(edge_colour = 'black') +
    ggraph::geom_node_point(ggplot2::aes(color = as.factor(community)), size = 5) +
    ggplot2::theme_bw()
  
  g_figure <- g_figure + 
    ggplot2::scale_color_discrete(
      name   = "",
      labels = labels
    )
}

g_figure2 <- g_figure + 
  geom_node_text(
    aes(label = gene_name),
    repel = TRUE,        # avoids text overlapping
    size = label_size
  )


ggsave( plot = g_figure, filename = paste0(directory, "/ppi_network_replot_", reference, ".png"), width = network_width, height = network_height )
ggsave( plot = g_figure2, filename = paste0(directory, "/ppi_network_with_labels_replot_", reference, ".png"), width = network_width, height = network_height )


}
  




  