library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the cnetplot function
cnetplot <- function(df, showCategory = 5, foldChange = NULL, layout = "kk", colorEdge = FALSE, circular = FALSE, node_label = "all", cex_category = 1, cex_gene = 1, cex_label_category = 1, cex_label_gene = 1) {
  # Ensure the dataframe has the required columns
  required_columns <- c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
  if (!all(required_columns %in% colnames(df))) {
    stop("Dataframe must contain the following columns: ", paste(required_columns, collapse = ", "))
  }
  # Clean the Term column by removing everything after "R-HSA"
  df$Term <- sub("R-HSA.*", "", df$Term)
  # Clean the Term column by removing everything after "R-HSA"
  df$Term <-sub(" \\(GO:.*", "", df$Term)
  # Process the dataframe to extract gene sets
  #df <- df %>% arrange(P.value) %>% head(showCategory)
  geneSets <- setNames(strsplit(as.character(df$Genes), ";"), df$Term)  # Split genes by ";"
  
  # Create a graph from the gene sets
  edges <- do.call(rbind, lapply(names(geneSets), function(term) data.frame(from = term, to = geneSets[[term]], stringsAsFactors = FALSE)))
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Add foldChange information if provided
  if (!is.null(foldChange)) {
    V(g)$foldChange <- foldChange[V(g)$name]
  }
  
  # Set node sizes and apply the size increase factor of 10
  V(g)$size <- 1
  V(g)$size[V(g)$name %in% names(geneSets)] <- sapply(geneSets, length) * 10
  
  # Highlight genes that are common in multiple pathways
  gene_count <- table(edges$to)
  common_genes <- names(gene_count[gene_count > 1])
  V(g)$color <- ifelse(V(g)$name %in% common_genes, "red", "orchid")
  
  # Determine edge color
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- geom_edge_link(aes(color = category), alpha = 0.8)
  } else {
    edge_layer <- geom_edge_link(color = 'darkgrey', alpha = 0.8)
  }
  
  # Create the ggraph plot
  p <- ggraph(g, layout = layout, circular = circular) +
    edge_layer +
    geom_node_point(aes(size = size, color = color)) +
    scale_size(range = c(3, 8) * cex_category) +
    theme_void() +  # Use theme_void() to remove the background and axes
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      axis.line = element_blank(),  # Remove axis lines
      axis.text = element_blank(),  # Remove axis text
      axis.ticks = element_blank(), # Remove axis ticks
      panel.grid = element_blank()  # Remove grid lines
    )
  
  # Add node labels
  if (node_label != "none") {
    label_data <- p$data
    if (node_label == "category") {
      label_data <- label_data[V(g)$name %in% names(geneSets), ]
    } else if (node_label == "gene") {
      label_data <- label_data[!V(g)$name %in% names(geneSets), ]
    }
    p <- p + geom_node_text(aes(label = name), repel = TRUE, size = 5 * cex_label_category, data = label_data)
  }
  
  return(p)
}

# Example usage:
# Assuming df is your dataframe and foldChange is a named vector of fold changes for genes
# df <- read.csv("your_data.csv")
# foldChange <- c("gene1" = 1.5, "gene2" = -0.8, "gene3" = 2.0, ...)
# p <- cnetplot(df, foldChange = foldChange)
# print(p)
# p <- cnetplot(enrichres, foldChange = enrichres$Odds.Ratio)
# print(p)
