#File to identified enriched terms in DE genes between neuronal classes using the pathfinR algorythm. 
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


protein_info_og <- read.csv(paste0(getwd(), "../expression_info/geneexpression_byregion_ALLEN.csv"))
hippocampus_sig <- read.xlsx(paste0(getwd(), "../expression_info/PSD_Regions_highest.xlsx"), colNames = TRUE)

#regions refer to neuronal classes
regions<- unique(hippocampus_sig$Highest_In)
regions <- regions[! regions == 'NS']

for (z in (length(regions)):(length(regions))){
  region<- regions[z]
  protein_info <- protein_info_og
  rhet_significant <- hippocampus_sig
  rhet_significant$qval[rhet_significant$Highest_In != regions[z]] <- 1
  rhet_significant <- rhet_significant[,-3]

  protein_info[,1] <- toupper(protein_info[,1])
  protein_info<- protein_info[!duplicated(protein_info[ , 1]),]
  
  rownames(protein_info) <- protein_info[,1]
  rhet_significant$Description <- toupper(rhet_significant$Description)
  
  protein_info <- protein_info[,-1]
  protein_info <- as.matrix(protein_info)
  
  zscaled_prot_info <-  protein_info
  
  dbs <- c("C2", "C5", "HPO")
  
  for (y in 1:length(dbs)){

    db = dbs[y]
    mmu_genes <- readRDS(paste0("mmu_MSigDB_", dbs[y], "_genes.RDS"))
    mmu_descriptions <- readRDS(paste0("mmu_MSigDB_", dbs[y], "_descriptions.RDS"))
    pathtosif<- "../musmusculus_pin_combined.sif"

    
    input_df<- data.frame(Gene.symbol= rhet_significant$Description,  adj.P.Val= rhet_significant$qval)
    input_df$Gene.symbol <- toupper(input_df$Gene.symbol)
    input_df$adj.P.Val <- as.numeric(input_df$adj.P.Val)
    input_df <- input_df[!is.na(input_df$adj.P.Val),]
    
    output_df <- run_pathfindR(input = input_df,
                               convert2alias = FALSE, iterations = 50,
                               gene_sets = "Custom",plot_enrichment_chart = FALSE,
                               custom_genes = mmu_genes,visualize_enriched_terms=FALSE,
                               custom_descriptions =  mmu_descriptions, list_active_snw_genes = TRUE, pin_name_path = pathtosif, output_dir = paste0("../Allen_Classes_Highest" , region, "_", db))
    
    
    write.csv2(output_df,paste0("../Allen_Classes_Highest" , region, "_", db, "/RAW_Allen_Classes_Highest", region, "_", db, ".csv"), row.names = FALSE)
    
    output_df$Up_regulated <- toupper(output_df$Up_regulated)
    output_df$Down_regulated <- toupper(output_df$Down_regulated)
    output_df <- output_df[lapply(str_split(output_df$Up_regulated, ","), length) >= 3,] #filter by length greater or equal than 3
    output_df <- output_df[output_df$occurrence >= 13,]
    
    if (length(output_df$ID) > 5){
      output_df <- cluster_enriched_terms(
        output_df,
        method = "hierarchical",
        plot_clusters_graph = FALSE,
        use_description = FALSE,
        use_active_snw_genes = FALSE
      )
    }
    
    write.csv2(output_df, paste0("../Allen_Classes_Highest" , region, "_", db, "/Allen_Classes_Highest", region, "_", db, ".csv"), row.names = FALSE)    
    
    dir.create(paste0("../Allen_Classes_Highest", regions[z], "_",dbs[y], "/heatmaps"))
    output_df <- read.csv2(paste0("../Allen_Classes_Highest" , region, "_", db, "/Allen_Classes_Highest", region, "_", db, ".csv"))
    
    enrichment_df<- output_df
    nterms_sig<- length(enrichment_df[,1])
    if (nterms_sig > 0){
    sig_genes <- input_df[input_df$adj.P.Val<=0.05,]
    
    for (i in 1:nterms_sig){
      
      term1 <- enrichment_df$ID[i]
      
      term1_desc <- enrichment_df$Term_Description[i]
      genest1 <- mmu_genes[[term1]]
      
      if (sum(toupper(genest1) %in%  row.names(zscaled_prot_info))>1){
        
        heatdf <- zscaled_prot_info[row.names(zscaled_prot_info) %in% toupper(genest1),]
        
        
        
        genes_sig_heat<- row.names(heatdf[row.names(heatdf) %in% sig_genes$Gene.symbol,,drop=F])
        genes_nosig_heat <- row.names(heatdf[!(row.names(heatdf) %in% sig_genes$Gene.symbol),, drop=F])
        
        sig_on_heat<- heatdf[row.names(heatdf) %in% sig_genes$Gene.symbol,,drop=F]
        nosig_on_heat<- heatdf[!(row.names(heatdf) %in% sig_genes$Gene.symbol),,drop=F]
        
        if (nrow(sig_on_heat) >1 ){
          sig_on_heat <- sig_on_heat[order(rowMeans(sig_on_heat[,1:3]),decreasing = T),]
        }
        
        if (nrow(nosig_on_heat)> 1 ){
          nosig_on_heat <- nosig_on_heat[order(rowMeans(nosig_on_heat[,1:3]),decreasing = T),]
        }
        heatdf <- rbind(sig_on_heat,nosig_on_heat)
        
        labs <- row.names(heatdf) %in% sig_genes$Gene.symbol
        labs <- factor(ifelse(labs == T, "Significant", "Not significant"), levels= c("Significant", "Not significant"))
        
        
        
        if (dim(heatdf)[1] > 200)  {
          h= 20
          l=4
          
        }else if (dim(heatdf)[1] > 100){  
          h=16
          l=4
        }else{
          h=10
          l=7
        }
        
        
        h1= Heatmap(heatdf, name="Z-Score", column_title = term1_desc, cluster_rows = FALSE, cluster_row_slices=FALSE, cluster_column_slices=FALSE, show_column_dend = TRUE, row_split = labs, row_gap=unit(2.5,"mm"), column_gap = unit(2.5,"mm"), column_split= factor(regions,levels= regions), row_names_gp = grid::gpar(fontsize = l))
        pdf(file=sprintf(paste(getwd(), "../Allen_Classes_Highest", regions[z], "_",dbs[y], "/heatmaps/%s_Allen_Highest_%s_%s.pdf", sep=""),  term1_desc, regions[z], dbs[y]  ), height = h)
        draw(h1)
        
        
        dev.off()
        
      }
    }
    
  }
  }
}
