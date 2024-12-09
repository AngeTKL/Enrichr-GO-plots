


# Define the function
clean_term_names <- function(term_names) {
  # Load required library
  library(stringr)
  
  # Remove "(GO:" and everything after it
  term_names <- str_replace(term_names, "\\s*\\(GO:.*$", "")
  
  # Remove "P00" and everything after it
  term_names <- str_replace(term_names, "P00.*$", "")
  
  # Trim any leading or trailing whitespace
  term_names <- str_trim(term_names)
  
  return(term_names)
}

# Function to compute similarity matrix
get_similarity_matrix <- function(terms) {
  term_matrix <- sapply(terms, function(term) {
    str_split(term, " ")[[1]]
  }, simplify = FALSE)
  
  # Create a binary matrix for presence/absence of unique words
  unique_words <- unique(unlist(term_matrix))
  binary_matrix <- sapply(terms, function(term) {
    presence_vector <- as.integer(unique_words %in% str_split(term, " ")[[1]])
    return(presence_vector)
  })
  
  
  
  # Convert the binary matrix to a distance matrix
  similarity_matrix <- 1 - as.matrix(dist(t(binary_matrix), method = "binary"))  # Transpose for proper orientation
  
  # Ensure that the similarity matrix is square and assign row/column names
  if (nrow(similarity_matrix) == length(terms) && ncol(similarity_matrix) == length(terms)) {
    rownames(similarity_matrix) <- terms
    colnames(similarity_matrix) <- terms
  } else {
    stop("The dimensions of the similarity matrix do not match the length of the terms.")
  }
 
  
  return(similarity_matrix)
}


# Function to determine the optimal number of clusters
determine_optimal_clusters <- function(similarity_matrix, max_clusters = 10) {
  
  # Load required libraries
  library(cluster)  # For silhouette analysis
  library(factoextra)  # For fviz_nbclust
  
  # Convert the similarity matrix to a distance matrix
  distance_matrix <- as.dist(1 - similarity_matrix)
  
  # Get the number of terms
  n_terms <- nrow(similarity_matrix)
  
  # Set a range of potential cluster numbers
  max_clusters <- min(max_clusters, n_terms-1)  # Limit max clusters to number of terms
  cluster_numbers <- 2:max_clusters  # At least 2 clusters
  
  # Elbow method to determine the optimal number of clusters
  wss <- sapply(cluster_numbers, function(k) {
    if (k > n_terms || k < 1) {
      return(NA)  # Return NA for invalid cluster number
    }
    kmeans(distance_matrix, centers = k, nstart = 25)$tot.withinss
  })
  
  # Plot the elbow method results
 p= plot(cluster_numbers, wss, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of Clusters", ylab = "Total Within-Cluster Sum of Squares",
       main = "Elbow Method for Optimal Clusters")
  
  # Silhouette method to determine the optimal number of clusters
  silhouette_scores <- sapply(cluster_numbers, function(k) {
    if (k > n_terms || k <= 1) {
      return(NA)  # Return NA for invalid cluster number
    }
    km <- kmeans(distance_matrix, centers = k, nstart = 25)
    silhouette_score <- mean(silhouette(km$cluster, distance_matrix)[, 3])
    return(silhouette_score)
  })
  
  # Plot silhouette method results
  plot(cluster_numbers, silhouette_scores, type = "b", pch = 19, col = "blue",
       frame = FALSE, xlab = "Number of Clusters", ylab = "Average Silhouette Score",
       main = "Silhouette Method for Optimal Clusters")
  
  # Determine the optimal number of clusters based on maximum silhouette score
  optimal_clusters <- cluster_numbers[which.max(silhouette_scores)]
  
  return(optimal_clusters)
}

shorten_label <- function(x, label_format = 30) {
  # Initialize an empty vector to store shortened labels
  shortened_labels <- character(length(x))
  
  for (i in seq_along(x)) {
    # Capitalize each word in the term
    capitalized <- tools::toTitleCase(x[i])
    
    # Check if the capitalized label is longer than the specified length
    if (nchar(capitalized) <= label_format) {
      shortened_labels[i] <- capitalized
    } else {
      # Cut the label to the specified length
      shortened_labels[i] <- substr(capitalized, 1, label_format)
      
      # Find the last space within the specified length to ensure it cuts at the last whole word
      last_space <- max(gregexpr(" ", shortened_labels[i])[[1]])
      
      # If there's a space, truncate at the last space
      if (last_space > 0) {
        shortened_labels[i] <- substr(shortened_labels[i], 1, last_space - 1)
      }
    }
  }
  
  return(shortened_labels)
}


# Main function to generate tree plot from enrichment results
treeplot.enrichResult <- function(enrich_result, 
                                  showCategory = 30,
                                  color = "Adjusted.P.value",
                                  fontsize = 8,
                                  label_format = 30,
                                  cluster_params = list(),
                                  hilight_params = list(hilight = TRUE, align = "both"),
                                  offset_params = list(tiplab = 0.1),
                                  ...) {
  
  
  library(ggtree)
  # Filter and prepare data
  enrich_result <- enrich_result[order(enrich_result[[color]]), ]  # Sort by color variable
  if (nrow(enrich_result) > showCategory) {
    enrich_result <- enrich_result[1:showCategory, ]  # Limit to top showCategory terms
  }
  # Clean the names of the terms
  enrich_result$Term <- clean_term_names(enrich_result$Term)
  
  # Calculate the number of genes associated with each term
  enrich_result$Count <- sapply(strsplit(enrich_result$Genes, ";"), length)
  
  
  p1 <- ggplot(enrich_result, aes(x = Count, y = reorder(Term, Count))) + 
    geom_bar(stat = "identity", fill = "skyblue") + 
    theme_minimal() +
    
    labs(x = "Number of Genes", y = "Enriched Terms", title = "Enrichment Analysis results :  Number of genes per enriched terms ") +
    theme(axis.text.y = element_text(size = 15),  # Increase size of y-axis text (terms)
          axis.text.x = element_text(size = 10))  # Increase size of x-axis text (counts)
  
  # Handle the case where there is only one term
  if (nrow(enrich_result) <=2) {
    message("Less than 2 terms present. Tree plot is not possible.")
    # Return a simple plot showing the single term
    p <- ggplot(enrich_result, aes(x = 1, y = 1)) +
      geom_text(aes(label = Term), size = 5) +
      theme_void() +
      ggtitle("Single Term: Tree plot  Not Applicable")
    return(p1)
  }
  
  
  
  
  # Compute the similarity matrix based on the cleaned terms
  termsim <- get_similarity_matrix(as.vector(enrich_result$Term))
  
  # Check dimensions
  dim(termsim)
  # Determine the optimal number of clusters
  optimal_clusters <- determine_optimal_clusters(termsim)

  # Check if optimal_clusters is NULL or less than 1
  # Check if optimal_clusters is empty, NA, NULL, or invalid
  if (length(optimal_clusters) == 0 || is.null(optimal_clusters) || is.na(optimal_clusters) || optimal_clusters < 1 || optimal_clusters >= nrow(enrich_result)) {
    message("Clustering not possible with the given data. Generating empty image  instead.")
    p <- ggplot() + theme_void() + ggtitle("No valid clusters could be generated.")
    return(p)
  }
  
  # Perform hierarchical clustering using Ward's method
  hc <- hclust(as.dist(1 - termsim), method = "ward.D2")
  
  # Cut the tree into the optimal number of clusters
  clus <- cutree(hc, k = optimal_clusters)
  
  # Prepare data for plotting
  plot_data <- data.frame(
    label = sapply(enrich_result$Term, function(x) shorten_label(x, label_format= 30)),  # Shorten labels
    color = enrich_result[[color]],  # Extract color variable
    count = enrich_result$Count,  # Number of genes per term
    cluster = factor(clus[hc$order])  # Order clusters according to the hierarchical clustering
  )
  
  # Plot the hierarchical tree with ggtree
  p <- ggtree(hc, layout = "rectangular") %<+% plot_data +  # Merge data with plot
    geom_tiplab(aes(label = label), size = fontsize, align = TRUE, offset = offset_params$tiplab) +  # Add labels
    geom_tippoint(aes(color = color, size = count), size = 4) +  # Add points representing gene count
    scale_color_continuous(name = color, low = "blue", high = "red") +  # Color scale based on the color variable
    theme_tree2() +  # Apply tree-specific theme
    theme(legend.position = "none" ,  #legend.position = "top", legend.text = element_text(size = 4),  # Reduced legend text size
          #legend.title = element_text(size = 8),  # Reduced legend title size
          #legend.key.size = unit(0.5, 'cm'),
          axis.text.x = element_blank(), axis.title.x = element_blank())  # Remove x-axis labels and title
  
  # Adjust the zoom and label size
  max_length <- 1.0  # Increase max_length to zoom out
  p <- p + coord_cartesian(xlim = c(-max_length, max_length * 1.5)) +  # Adjust plot limits
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))  # Set plot margins
  
  # Highlight clusters with similar colors
  if (hilight_params$hilight) {
    unique_clusters <- unique(plot_data$cluster)  # Identify unique clusters
    group_color <- rainbow(length(unique_clusters))  # Generate colors for each cluster
    names(group_color) <- unique_clusters  # Assign colors to clusters
    
    p <- p + geom_hilight(node = which(hc$order %in% which(clus == unique_clusters[1])), 
                          fill = group_color[1], alpha = 0.3) +  # Highlight first cluster
      lapply(unique_clusters[-1], function(cluster) {  # Loop through other clusters
        geom_hilight(node = which(hc$order %in% which(clus == cluster)), 
                     fill = group_color[cluster], alpha = 0.3)  # Highlight each cluster
      })
  }
  
  # # Add legend without cluster labels
  # p <- p + labs(size = "Number of Genes"#, color = "Adjusted P-value"
  #               ) +  # Add labels to the legend
  #   theme(legend.position = "top")  # Position the legend
  
  # Return the plot and data
  return( plot=p) 
}




# Main function to generate tree plot from enrichment results
barplot.enrichResult <- function(enrich_result, 
                                  showCategory = 30,
                                  color = "Adjusted.P.value",
                                  fontsize = 8,
                                  label_format = 30,
                                  cluster_params = list(),
                                  hilight_params = list(hilight = TRUE, align = "both"),
                                  offset_params = list(tiplab = 0.1),
                                  ...) {
  
  
  library(ggtree)
  # Filter and prepare data
  enrich_result <- enrich_result[order(enrich_result[[color]]), ]  # Sort by color variable
  if (nrow(enrich_result) > showCategory) {
    enrich_result <- enrich_result[1:showCategory, ]  # Limit to top showCategory terms
  }
  # Clean the names of the terms
  enrich_result$Term <- clean_term_names(enrich_result$Term)
  
  # Calculate the number of genes associated with each term
  enrich_result$Count <- sapply(strsplit(enrich_result$Genes, ";"), length)
  
  
  p1 <- ggplot(enrich_result, aes(x = Count, y = reorder(Term, Count))) + 
    geom_bar(stat = "identity", fill = "skyblue") + 
    # Add text labels showing the counts on top of the bars
    geom_text(aes(label = Count), hjust = -0.2, size = 5) + 
    theme_minimal() +
    labs(x = "Number of Genes", y = "Enriched Terms", title = "Enrichment Analysis results :  Number of genes per enriched terms ") +
    theme(axis.text.y = element_text(size = 15),  # Increase size of y-axis text (terms)
          axis.text.x = element_text(size = 10)) + # Increase size of x-axis text (counts)
  
     # Increase size of x-axis text (counts)
    # Ensure x-axis only shows whole numbers
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = c(0, 0.1)))
  
  return(p1)
}
