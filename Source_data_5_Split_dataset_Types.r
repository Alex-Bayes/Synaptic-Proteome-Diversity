#Code to split Allen Brain atlas transcriptomic data from excitatory neurons.
library(stringr)
library(dplyr)

samps = readRDS("samps.RDS")
prot_df = readRDS("prot_PSD_df.RDS")
annot = readRDS("annot.RDS") 

annot_CA1 = annot[grep("CA1$",annot$cluster_label),]
samples_CA1 = as.vector(annot_CA1$sample_name)
samps_CA1 = data.frame(annot_CA1$sample_name)
colnames(samps_CA1) = "sample_name"

annot_CA1PROS = annot[grep("CA1-ProS$",annot$cluster_label),]
samples_CA1PROS = as.vector(annot_CA1PROS$sample_name)
samps_CA1PROS = data.frame(annot_CA1PROS$sample_name)
colnames(samps_CA1PROS) = "sample_name"

annot_CA2 = annot[grep("CA2-IG-FC$",annot$cluster_label),]
samples_CA2 = as.vector(annot_CA2$sample_name)
samps_CA2 = data.frame(annot_CA2$sample_name)
colnames(samps_CA2) = "sample_name"

annot_CA3 = annot[grep("CA3$",annot$cluster_label),]
samples_CA3 = as.vector(annot_CA3$sample_name)
samps_CA3 = data.frame(annot_CA3$sample_name)
colnames(samps_CA3) = "sample_name"

annot_DG = annot[grep("DG$",annot$cluster_label),]
samples_DG = as.vector(annot_DG$sample_name)
samps_DG = data.frame(annot_DG$sample_name)
colnames(samps_DG) = "sample_name"

annot_CTSUB = annot[grep("CT SUB$",annot$cluster_label),]
samples_CTSUB = as.vector(annot_CTSUB$sample_name)
samps_CTSUB = data.frame(annot_CTSUB$sample_name)
colnames(samps_CTSUB) = "sample_name"

annot_NPSUB = annot[grep("NP SUB$",annot$cluster_label),]
samples_NPSUB = as.vector(annot_NPSUB$sample_name)
samps_NPSUB = data.frame(annot_NPSUB$sample_name)
colnames(samps_NPSUB) = "sample_name"

annot_SUBPROS = annot[grep("SUB-ProS$",annot$subclass_label),]
samples_SUBPROS = as.vector(annot_SUBPROS$sample_name)
samps_SUBPROS = data.frame(annot_SUBPROS$sample_name)
colnames(samps_SUBPROS) = "sample_name"

prot_df_CA1 = select(prot_df,one_of(samples_CA1))
prot_df_CA1PROS = select(prot_df,one_of(samples_CA1PROS))
prot_df_CA2 = select(prot_df,one_of(samples_CA2))
prot_df_CA3 = select(prot_df,one_of(samples_CA3))
prot_df_DG = select(prot_df,one_of(samples_DG))
prot_df_CTSUB = select(prot_df,one_of(samples_CTSUB))
prot_df_NPSUB = select(prot_df,one_of(samples_NPSUB))
prot_df_SUBPROS = select(prot_df,one_of(samples_SUBPROS))

saveRDS(object = prot_df_CA1,file = "prot_PSD_df_CA1.RDS")
saveRDS(object = annot_CA1,file = "annot_CA1.RDS")
saveRDS(object = samps_CA1,file = "samps_CA1.RDS")

saveRDS(object = prot_df_CA1PROS,file = "prot_PSD_df_CA1PROS.RDS")
saveRDS(object = annot_CA1PROS,file = "annot_CA1PROS.RDS")
saveRDS(object = samps_CA1PROS,file = "samps_CA1PROS.RDS")

saveRDS(object = prot_df_CA2,file = "prot_PSD_df_CA2.RDS")
saveRDS(object = annot_CA2,file = "annot_CA2.RDS")
saveRDS(object = samps_CA2,file = "samps_CA2.RDS")

saveRDS(object = prot_df_CA3,file = "prot_PSD_df_CA3.RDS")
saveRDS(object = annot_CA3,file = "annot_CA3.RDS")
saveRDS(object = samps_CA3,file = "samps_CA3.RDS")

saveRDS(object = prot_df_DG,file = "prot_PSD_df_DG.RDS")
saveRDS(object = annot_DG,file = "annot_DG.RDS")
saveRDS(object = samps_DG,file = "samps_DG.RDS")

saveRDS(object = prot_df_CTSUB,file = "prot_PSD_df_CTSUB.RDS")
saveRDS(object = annot_CTSUB,file = "annot_CTSUB.RDS")
saveRDS(object = samps_CTSUB,file = "samps_CTSUB.RDS")

saveRDS(object = prot_df_NPSUB,file = "prot_PSD_df_NPSUB.RDS")
saveRDS(object = annot_NPSUB,file = "annot_NPSUB.RDS")
saveRDS(object = samps_NPSUB,file = "samps_NPSUB.RDS")

saveRDS(object = prot_df_SUBPROS,file = "prot_PSD_df_SUBPROS.RDS")
saveRDS(object = annot_SUBPROS,file = "annot_SUBPROS.RDS")
saveRDS(object = samps_SUBPROS,file = "samps_SUBPROS.RDS")

