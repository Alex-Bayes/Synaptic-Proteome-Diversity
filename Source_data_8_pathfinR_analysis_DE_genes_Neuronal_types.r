#File to identified enriched terms in DE genes between neuronal types using the pathfinR algorythm. 
library(pathfindR)
library(openxlsx)
library(biomaRt)

library('org.Hs.eg.db')
library(writexl)
library(plyr)
library(tidyverse)

library(ComplexHeatmap)
library(circlize)
library(STRINGdb)
library(igraph)

#regions refer to neuronal classes
og_regions <- c("CA1", "CA1_ProS", "CA2", "CA3", "DG", "CTSUB", "NPSUB", "SUB_ProS")

for (t in 1:8){
protein_info_og <- read.csv(paste0(getwd(), "../expression_info/neuron_df_PSD_", og_regions[t], ".csv"))
hippocampus_sig <- read.xlsx(paste0(getwd(), "../expression_info/Taula_Resum_PSD_Tipus.xlsx"), colNames = TRUE, sheet= t)

regions<- unique(hippocampus_sig$Neuronal.type)


for (z in 1:length(regions)){

region<- regions[z]
protein_info <- protein_info_og
rhet_significant <- hippocampus_sig
rhet_significant$P.value[rhet_significant$Neuronal.type != regions[z]] <- 1
rhet_significant <- rhet_significant[,-3]

protein_info[,1] <- toupper(protein_info[,1])

rownames(protein_info) <- protein_info[,1]
rhet_significant$Gene.Name <- toupper(rhet_significant$Gene.Name)

protein_info <- protein_info[,c(-1,-2)]

protein_info <- as.matrix(protein_info)
zscaled_prot_info<- log2(protein_info + 1)

dbs <- c("C2", "C5")

for (y in 1:length(dbs)){
mmu_genes <- readRDS(paste0("mmu_MSigDB_", dbs[y], "_genes.RDS"))
mmu_descriptions <- readRDS(paste0("mmu_MSigDB_", dbs[y], "_descriptions.RDS"))
pathtosif<- paste0(getwd(),"/musmusculus_pin_combined.sif")

input_df<- data.frame(Gene.symbol= rhet_significant$Gene.Name,  adj.P.Val= rhet_significant$P.value)
input_df$Gene.symbol <- toupper(input_df$Gene.symbol)
input_df$adj.P.Val <- as.numeric(input_df$adj.P.Val)


if (dim(input_df[input_df$adj.P.Val<=0.05,])[1] >= 3){
output_df <- run_pathfindR(input = input_df,
                           convert2alias = FALSE, iterations = 50,
                           gene_sets = "Custom", plot_enrichment_chart = FALSE,
                           custom_genes = mmu_genes,
                           custom_descriptions =  mmu_descriptions, list_active_snw_genes = FALSE, pin_name_path = pathtosif, 
                           output_dir = paste0("Allen_TypesUnique", regions[z], "_",dbs[y] ))

write.csv(output_df, paste0("../Allen_TypesUnique", regions[z], "_",dbs[y],"/RAW_Allen_TypesUnique", regions[z], "_",dbs[y],".csv"), row.names = F)

if (dim(output_df)[1] != 0){
output_df$Up_regulated <- toupper(output_df$Up_regulated)
output_df$Down_regulated <- toupper(output_df$Down_regulated)
output_df <- output_df[lapply(str_split(output_df$Up_regulated, ","), length) >= 3,] #filter by length greater or equal than 3
output_df <- output_df[output_df$occurrence >= 13,]
}
if (length(output_df$ID) > 5){
output_df <- cluster_enriched_terms(
  output_df,
  method = "hierarchical",
  plot_clusters_graph = FALSE,
  use_description = FALSE,
  use_active_snw_genes = FALSE
)
}
write.csv(output_df, paste0("../Allen_TypesUnique", regions[z], "_",dbs[y],"/Allen_TypesUnique", regions[z], "_",dbs[y],".csv"), row.names = F)


}
}
}
}

