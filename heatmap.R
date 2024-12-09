##############- HEATMAP ###################
# Sample data creation (replace this with your actual data)
#library(data.table)
# set working directory 
#setwd('K:/UKB Projects/Hopewell Cardiovascular Therapies/Projects/AF_Diabetes_MR/Analysis/EUR') 
#dt =fread("K:/UKB Projects/Hopewell Cardiovascular Therapies/Projects/AF_Diabetes_MR/Analysis/EUR/MR/Output/Results/TableAllclustersDM_AF_basic.xlsx")
# Assuming your dataset is named 'data' and you want to rename the rows where the 'method' column has the value "Inverse variance weighted (fixed effects)"
#data[data$method == "Inverse variance weighted (fixed effects)", "method"] <- "IVW fixed effects"
#data[data$method == "Inverse variance weighted (multiplicative random effects)", "method"] <- "IVW random effects"



#NEED y=Method, x= exposure, fill = log(pval))) 
heatmap = function (data){
# Load necessary library
library(ggplot2)

  # Create a new column for annotations based on p-values
  data$annotation <- ifelse(data$pval < 0.05, "*", 
                            ifelse(data$pval < 0.01, "**", 
                                   ifelse(data$pval < 0.001, "***", "")))

# Convert Method and Cluster to factors to ensure proper ordering in the heatmap
data$method   <- factor(data$method  , levels = unique(data$method))
data$exposure <- factor(data$exposure, levels = unique(data$exposure))

# Create the heatmap
heatmap_plot <- ggplot(data, aes(x = exposure, y = method  , fill = log(pval))) +
  geom_tile(color = "white") +
 # geom_text(aes(label = annotation), color = "black", size = 4) +  # Add annotations 
  scale_fill_gradient(low = "green", high = "pink", name = "P-Value") +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1, , lineheight = 0.7),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 14)
       ) +
  labs(title = "Heatmap of Methods and Clusters Based on P-Values",
       x = "Cluster",
       y = "Method") +
  geom_text(aes(label = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", ""))), 
                #color = "black"
                ), size = 3, hjust = 0.5, vjust = 0.5)

 
# Print the heatmap
print(heatmap_plot)

}


#heatmap(dt)


#plot()

