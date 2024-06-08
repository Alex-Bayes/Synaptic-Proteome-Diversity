#File to identify differentially expressed genes in a type of excitatory neurons.
library(Seurat)
library(stringr)
library(dplyr)
library(limma)
library(doParallel)
library(foreach)

samps = readRDS("samps_CA3.RDS") # Sample names.
prot_df = readRDS("prot_PSD_df_CA3.RDS") # Dataframe of protein detected per sample.
annot = readRDS("annot_CA3.RDS") # Sample annotation.

# Iteration process.
random_iters = function(samps,annot,prot_df,j,rand_samps,results_list = list()) {
  sub_annot = annot[annot$sample_name %in% rand_samps[[j]],]
  # Neuron name.
  sub_samps = as.data.frame(samps[samps$sample_name %in% rand_samps[[j]],])
  colnames(sub_samps) = "sample_name"
  # Neuron expression.
  sub_prot_df = prot_df[,colnames(prot_df) %in% rand_samps[[j]]]

  sub_annot = as.data.frame(sub_annot[order(sub_annot$sample_name),])
  sub_prot_df = sub_prot_df[,order(colnames(sub_prot_df))]
  row.names(sub_annot) = sub_annot$sample_name

  # Data analysis with Seurat package.
  pbmc = CreateSeuratObject(sub_prot_df,min.cells = 3,min.genes = 200,project = "RNA-SeqHip",
                            meta.data = sub_annot) # Seurat object.
  pbmc = NormalizeData(pbmc) # Data normalization.

  # Differential expression analysis.
  labs = unique(pbmc@meta.data[["cluster_label"]]) # Labeling neuronal types.

  # DE analysis comparing each neuron type with the others. 
  for (k in 1:length((labs))){
    deg_res = FindMarkers(pbmc,group.by = "cluster_label",ident.1 = labs[k],
                          logfc.threshold = 0,min.pct = 0)

    results_list[[length(results_list) + 1]] = deg_res
    type = as.character(labs[k])
    iteration = as.character(j)
    name = str_c(type,iteration,sep = "-")
    names(results_list)[length(results_list)] = name
  }
  return(results_list)
}

# Iterations number.
n = 150
rand_samps = lapply(1:n, function(n) annot %>%
                      group_by(cluster_label) %>%
                      sample_n(size = 25,replace = F) %>%
                      pull(sample_name))


# Calling the function.
cl <- makeCluster(10)
registerDoParallel(cl)
results_list <- foreach(i=1:n,
                            .packages = c("Seurat","dplyr","limma","stringr" )) %dopar% random_iters(samps,annot,prot_df,i,rand_samps)
stopCluster(cl)

saveRDS(object = results_list,file = "Output_Types.RDS")









