#enrichr 
#install_github("guokai8/EnrichR")
#remotes::install_github("guokai8/Enrichr")

library(devtools)
library(enrichR)
listEnrichrSites()
#> OxEnrichr ... Connection is available!
setEnrichrSite("Enrichr") # Human genes
#CALL PLOT FUNCTION 
#for heatplot 
source("K:/UKB Projects/Hopewell Cardiovascular Therapies/Projects/AF_Diabetes_MR/Code/R_programmes/MR_analysis/Basic functions/heatplot.R")
#for 
source("K:/UKB Projects/Hopewell Cardiovascular Therapies/Projects/AF_Diabetes_MR/Code/R_programmes/MR_analysis/Basic functions/cnet_plot.R")
#map plot 
source("K:/UKB Projects/Hopewell Cardiovascular Therapies/Projects/AF_Diabetes_MR/Code/R_programmes/MR_analysis/Basic functions/mapplot.R")
#Then find the list of all available databases from Enrichr.

dbs <- listEnrichrDbs()
head(dbs)

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022")
enriched <- enrichr(candidate_list, dbs)
# i can view the results here 

#Plot Enrichr GO-BP output. (Plotting function contributed by I-Hsuan Lin)
#bar plot 
plotEnrich(enriched[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")

enrichres= subset(enriched[[4]], enriched[[4]]$Adjusted.P.value<0.05)

plotEnrich(enrichres, showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
p <- cnetplot(enrichres, foldChange = enrichres$Odds.Ratio)
print(p)
#heat plot 
p <- heatplot(enrichres , foldChange = enrichres$Odds.Ratio)
print(p)

# Enrichment map 
emapplot(enrichres)
