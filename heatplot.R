##'Explanation:
##'Function Definition: heatplot.gseaResult is defined to create a heat plot directly from gseaResult object.
##'Data Processing: The function extracts genes, creates a matrix indicating gene presence across categories, and converts this to a long format suitable for ggplot.
#Extracts top categories/terms based on showCategory.
#Splits genes and creates a matrix (gene_matrix) indicating gene presence across categories.
#Converts gene_matrix to a dataframe (gene_df) suitable for plotting with ggplot2.
#Reshapes gene_df from wide to long format (gene_df_long) using tidyr::gather.
#Plotting:
#  Creates a heatmap using ggplot2.
#Cells are colored based on gene presence (Presence) or optionally based on foldChange values if provided.
#Customization:
#  Adjusts plot appearance including labels, theme (theme_minimal), and axis text orientation.
##'Usage:
##'To use this function, assuming enrichres is your dataframe:
# Required libraries
# Required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)  # For color palette
# Required libraries


# Define the heatplot.gseaResult function
heatplot <- function(x, showCategory = 30, foldChange = NULL) {
  
  # Extract the top categories
  top_categories <- x[1:showCategory, ]
  
  # Split genes and create a data frame with gene presence
  gene_df <- top_categories %>%
    mutate(Genes = strsplit(as.character(Genes), ";")) %>%
    unnest(Genes) %>%
    distinct(Term, Genes) %>%
    mutate(Presence = 1)
  # Remove " R-HSA-" and everything after it from Term
  gene_df$Term <- sub(" R-HSA-.*", "", gene_df$Term)
  # Clean the Term column by removing everything after "R-HSA"
  gene_df$Term <- sub(" \\(GO:.*", "", gene_df$Term)
  # Remove any rows where Term is NA (if any)
  gene_df <- gene_df[!is.na(gene_df$Term), ]
  
  # Create gene presence heatmap
  p = ggplot(gene_df, aes(x = Genes, y = Term, fill = Term)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = viridis_pal(option = "plasma")(length(unique(gene_df$Genes))),
                      name = "Genes",
                      breaks = unique(gene_df$Term),guide = "none") +
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=15),
          axis.text.y = element_text(size = 30),
          plot.margin = margin(10, 40, 20, 20, "pt"),  # Adjust plot margins
          plot.background = element_rect(fill = "white")) +  # Set plot background color
    labs(x = NULL, y = NULL)
    
  
  # Optionally add foldChange gradient
  if (!is.null(foldChange)) {
    fold_df <- gene_df %>%
      filter(Presence == 1) %>%
      mutate(foldChange = ifelse(Presence == 1, foldChange[Genes], NA_real_))
    
    p <- ggplot(fold_df, aes(x = Genes, y = Term, fill = foldChange)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = median(na.omit(foldChange)),
                           name = "Fold Change", guide = "none") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size=15),
            axis.title.y = element_text(size = 15),  # Optional: Increase y-axis title size
            plot.margin = margin(10, 40, 20, 20, "pt"),  # Adjust plot margins
            axis.text.y = element_text(size = 12),
            plot.background = element_rect(fill = "white")) +  # Set plot background color
      labs(x = NULL, y = NULL)
  } 
 
  
  # Save the plot with adjusted size
  ggsave("heatplot.png", plot = p, width = 10, height = 6, dpi=150,  units = "cm", limitsize = FALSE)
  
  return(p)
    return(p)
  } 






#p <- heatplot.gseaResult(enrichres , foldChange = enrichres$Odds.Ratio)
#print(p)
