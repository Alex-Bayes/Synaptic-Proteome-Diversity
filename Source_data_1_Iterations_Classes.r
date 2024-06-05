#File to identify differentially expressed genes in a class of excitatory neurons. 
library(Seurat)
library(stringr)
library(dplyr)
library(limma)
library(doParallel)
library(foreach)

samps = readRDS("samps_1.RDS") 
prot_df = readRDS("neuron_df_PSD_1.RDS") 
annot = readRDS("annot_1.RDS") 

samples_CA1 = grep("CA1$",annot$cluster_label)
annot$subclass_label[samples_CA1] = "CA1"

random_iters = function(samps,annot,prot_df,j,rand_samps,results_list = list()) {
  sub_annot = annot[annot$sample_name %in% rand_samps[[j]],]
  sub_samps = as.data.frame(samps[samps$sample_name %in% rand_samps[[j]],])
  colnames(sub_samps) = "sample_name"
  sub_prot_df = prot_df[,colnames(prot_df) %in% rand_samps[[j]]]

  sub_annot = as.data.frame(sub_annot[order(sub_annot$sample_name),])
  sub_prot_df = sub_prot_df[,order(colnames(sub_prot_df))]
  row.names(sub_annot) = sub_annot$sample_name

  pbmc = CreateSeuratObject(sub_prot_df,min.cells = 3,min.genes = 200,project = "RNA-SeqHip",
                            meta.data = sub_annot) 
  pbmc = NormalizeData(pbmc) 

  labs = unique(pbmc@meta.data[["subclass_label"]]) 

  for (k in 1:length((labs))){
    deg_res = FindMarkers(pbmc,group.by = "subclass_label",ident.1 = labs[k],
                          logfc.threshold = 0,min.pct = 0)

    results_list[[length(results_list) + 1]] = deg_res
    n_class = as.character(labs[k])
    iteration = as.character(j)
    name = str_c(n_class,iteration,sep = "-")
    names(results_list)[length(results_list)] = name
  }
  return(results_list)
}

n = 150
rand_samps = lapply(1:n, function(n) annot %>%
                      group_by(subclass_label) %>%
                      sample_n(size = 100,replace = F) %>%
                      pull(sample_name))


cl <- makeCluster(10)
registerDoParallel(cl)
results_list <- foreach(i=1:n,
                            .packages = c("Seurat","dplyr","limma","stringr" )) %dopar% random_iters(samps,annot,prot_df,i,rand_samps)
stopCluster(cl)

saveRDS(object = results_list,file = "Output_Classes.RDS")









