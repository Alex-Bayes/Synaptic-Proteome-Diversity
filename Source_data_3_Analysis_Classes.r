#File to analyze differentially expressed genes between neuronal classes
library(dplyr)
library(stringr)

results_list = readRDS("Output_Classes.RDS")
results_list = unlist(results_list,recursive = FALSE)

CA1_l = results_list[grep("^CA1-\\d",names(results_list))]
CA1PROS_l = results_list[grep("^CA1-ProS-",names(results_list))]
CA2_l = results_list[grep("^CA2-",names(results_list))]
CA3_l = results_list[grep("^CA3-",names(results_list))]
DG_l = results_list[grep("^DG-",names(results_list))]
CTSUB_l = results_list[grep("^CT SUB-",names(results_list))]
NPSUB_l = results_list[grep("^NP SUB-",names(results_list))]
SUBPROS_l = results_list[grep("^SUB-ProS-",names(results_list))]

p_ajust = function(x){
  p_val_raw = as.numeric(unlist(x$p_val))
  p_val_bh = p.adjust(p = p_val_raw,method = "BH")
  x = cbind(x,p_val_bh)
}

results_list = lapply(results_list,p_ajust)

res_raw = unlist(results_list,recursive = F)
fc_total = res_raw[grep("avg_log2FC",names(res_raw))]
fc_total = unlist(fc_total,recursive = F)
pv_total = res_raw[grep("p_val_bh",names(res_raw))]
pv_total = unlist(pv_total,recursive = F)

mean(fc_total)
mean(pv_total)

png(filename = "FC_Distribution.png",width = 800,height = 600)
distr_fc = hist(fc_total,main = "FC_Distribution",col = "lightsalmon",
                xlab = "Log2 FC",ylab = "Frequency")
dev.off()

png(filename = "PV_Distribution.png",width = 800,height = 600)
distr_pv = hist(pv_total,main = "PV_Distribution",col = "lightsalmon",
                xlab = "Adjusted p-value",ylab = "Frequency")
dev.off()

filtration = function(x,final_gene_set = list()){
  genes_filter = filter(x,(avg_log2FC > 0.6) | (avg_log2FC < -0.6) & (p_val_bh < 0.05))

  final_gene_set = rbind(final_gene_set,genes_filter)
  return(final_gene_set)
}

sig_df = lapply(results_list,filtration)

n = length(sig_df)
count_df = data.frame(matrix(ncol = 8,nrow = n/8))
classes = c("CA1","CA1-ProS","CA2","CA3","DG","CT SUB","NP SUB","SUB-ProS")
colnames(count_df) = classes

CA1_sig = sig_df[grep("^CA1-\\d",names(sig_df))]
CA1PROS_sig = sig_df[grep("^CA1-ProS-",names(sig_df))]
CA2_sig = sig_df[grep("^CA2-",names(sig_df))]
CA3_sig = sig_df[grep("^CA3-",names(sig_df))]
DG_sig = sig_df[grep("^DG-",names(sig_df))]
CTSUB_sig = sig_df[grep("^CT SUB-",names(sig_df))]
NPSUB_sig = sig_df[grep("^NP SUB-",names(sig_df))]
SUBPROS_sig = sig_df[grep("^SUB-ProS-",names(sig_df))]

sig_gens_CA1 =  unique(unlist(sapply(CA1_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA1PROS =  unique(unlist(sapply(CA1PROS_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA2 =  unique(unlist(sapply(CA2_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3 =  unique(unlist(sapply(CA3_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_DG =  unique(unlist(sapply(DG_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CTSUB =  unique(unlist(sapply(CTSUB_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_NPSUB =  unique(unlist(sapply(NPSUB_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_SUBPROS =  unique(unlist(sapply(SUBPROS_sig,function(x) row.names(x))),use.names=FALSE)

sig_gens_total = unique(c(sig_gens_CA1,sig_gens_CA1PROS,sig_gens_CA2,sig_gens_CA3,sig_gens_DG,
                          sig_gens_CTSUB,sig_gens_NPSUB,sig_gens_SUBPROS))

gens_totals = unique(unlist(sapply(results_list,function(x) row.names(x))),use.names=FALSE)
gens_never_sig = gens_totals[!(gens_totals%in%sig_gens_total)]

count_df$CA1 = sapply(CA1_sig,function(x) nrow(x))
count_df$`CA1-ProS` = sapply(CA1PROS_sig,function(x) nrow(x))
count_df$CA2 = sapply(CA2_sig,function(x) nrow(x))
count_df$CA3 = sapply(CA3_sig,function(x) nrow(x))
count_df$DG = sapply(DG_sig,function(x) nrow(x))
count_df$`CT SUB` = sapply(CTSUB_sig,function(x) nrow(x))
count_df$`NP SUB` = sapply(NPSUB_sig,function(x) nrow(x))
count_df$`SUB-ProS` = sapply(SUBPROS_sig,function(x) nrow(x))

genecount_df = data.frame(matrix(data = 0,ncol = 8,nrow = length(sig_gens_total)))
row.names(genecount_df) = sig_gens_total
colnames(genecount_df) = classes

for (i in 1:length(sig_gens_total)) {
  genecount_df$CA1[i]= sum(sapply(CA1_sig,function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`CA1-ProS`[i]= sum(sapply(CA1PROS_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$CA2[i]= sum(sapply(CA2_sig,function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$CA3[i]= sum(sapply(CA3_sig,function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$DG[i]= sum(sapply(DG_sig,function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`CT SUB`[i]= sum(sapply(CTSUB_sig,function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`NP SUB`[i]= sum(sapply(NPSUB_sig,function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`SUB-ProS`[i]= sum(sapply(SUBPROS_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
}

iterations = nrow(count_df) 
threshold = 0.9 * iterations 
gens_threshold_CA1 = rownames(genecount_df[which(genecount_df$CA1 >= threshold),])
gens_threshold_CA1PROS = rownames(genecount_df[which(genecount_df$`CA1-ProS` >= threshold),])
gens_threshold_CA2  = rownames(genecount_df[which(genecount_df$CA2 >= threshold),])
gens_threshold_CA3  = rownames(genecount_df[which(genecount_df$CA3 >= threshold),])
gens_threshold_DG  = rownames(genecount_df[which(genecount_df$DG >= threshold),])
gens_threshold_CTSUB  = rownames(genecount_df[which(genecount_df$`CT SUB` >= threshold),])
gens_threshold_NPSUB  = rownames(genecount_df[which(genecount_df$`NP SUB` >= threshold),])
gens_threshold_SUBPROS  = rownames(genecount_df[which(genecount_df$`SUB-ProS` >= threshold),])

# CA1.
class_gens_CA1 = function(y,gens_CA1 = data.frame(Gene = character(),
                                                  Express = character(),
                                                  FC = numeric(),
                                                  PV = numeric(),
                                                  stringsAsFactors = F)){
  ex = as.data.frame(sapply(CA1_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_CA1[nrow(gens_CA1) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA1[nrow(gens_CA1) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA1)
}

gens_class_CA1 = lapply(gens_threshold_CA1,class_gens_CA1)

gens_class_CA1 = data.frame(matrix(unlist(gens_class_CA1),nrow = length(gens_class_CA1),
                                   byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA1) = c("Gene","Express","Log2FC","P-value")

CA1_up = length(gens_class_CA1$Express[gens_class_CA1$Express == "UP"])
CA1_down = length(gens_class_CA1$Express[gens_class_CA1$Express == "DOWN"])

# CA1-ProS.
class_gens_CA1PROS = function(y,gens_CA1PROS = data.frame(Gen = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  ex = as.data.frame(sapply(CA1PROS_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_CA1PROS[nrow(gens_CA1PROS) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA1PROS[nrow(gens_CA1PROS) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA1PROS)
}

gens_class_CA1PROS = lapply(gens_threshold_CA1PROS,class_gens_CA1PROS)

gens_class_CA1PROS = data.frame(matrix(unlist(gens_class_CA1PROS),
                                       nrow = length(gens_class_CA1PROS),byrow = TRUE),
                                stringsAsFactors = F)
colnames(gens_class_CA1PROS) = c("Gene","Express","Log2FC","P-valor")

CA1PROS_up = length(gens_class_CA1PROS$Express[gens_class_CA1PROS$Express == "UP"])
CA1PROS_down = length(gens_class_CA1PROS$Express[gens_class_CA1PROS$Express == "DOWN"])

# CA2.
class_gens_CA2 = function(y,gens_CA2 = data.frame(Gen = character(),
                                                  Express = character(),
                                                  FC = numeric(),
                                                  PV = numeric(),
                                                  stringsAsFactors = F)){
  ex = as.data.frame(sapply(CA2_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_CA2[nrow(gens_CA2) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA2[nrow(gens_CA2) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA2)
}

gens_class_CA2 = lapply(gens_threshold_CA2,class_gens_CA2)

gens_class_CA2 = data.frame(matrix(unlist(gens_class_CA2),nrow = length(gens_class_CA2),
                                   byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA2) = c("Gene","Express","Log2FC","P-valor")

CA2_up = length(gens_class_CA2$Express[gens_class_CA2$Express == "UP"])
CA2_down = length(gens_class_CA2$Express[gens_class_CA2$Express == "DOWN"])

# CA3.
class_gens_CA3 = function(y,gens_CA3 = data.frame(Gen = character(),
                                                  Express = character(),
                                                  FC = numeric(),
                                                  PV = numeric(),
                                                  stringsAsFactors = F)){
  ex = as.data.frame(sapply(CA3_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_CA3[nrow(gens_CA3) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3[nrow(gens_CA3) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3)
}

gens_class_CA3 = lapply(gens_threshold_CA3,class_gens_CA3)

gens_class_CA3 = data.frame(matrix(unlist(gens_class_CA3),nrow = length(gens_class_CA3),
                                   byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3) = c("Gene","Express","Log2FC","P-valor")

CA3_up = length(gens_class_CA3$Express[gens_class_CA3$Express == "UP"])
CA3_down = length(gens_class_CA3$Express[gens_class_CA3$Express == "DOWN"])

# DG.
class_gens_DG = function(y,gens_DG = data.frame(Gen = character(),
                                                Express = character(),
                                                FC = numeric(),
                                                PV = numeric(),
                                                stringsAsFactors = F)){
  ex = as.data.frame(sapply(DG_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_DG[nrow(gens_DG) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_DG[nrow(gens_DG) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_DG)
}

gens_class_DG = lapply(gens_threshold_DG,class_gens_DG)

gens_class_DG = data.frame(matrix(unlist(gens_class_DG),nrow = length(gens_class_DG),
                                   byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_DG) = c("Gene","Express","Log2FC","P-valor")

DG_up = length(gens_class_DG$Express[gens_class_DG$Express == "UP"])
DG_down = length(gens_class_DG$Express[gens_class_DG$Express == "DOWN"])

# CT SUB.
class_gens_CTSUB = function(y,gens_CTSUB = data.frame(Gen = character(),
                                                      Express = character(),
                                                      FC = numeric(),
                                                      PV = numeric(),
                                                      stringsAsFactors = F)){
  ex = as.data.frame(sapply(CTSUB_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_CTSUB[nrow(gens_CTSUB) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CTSUB[nrow(gens_CTSUB) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CTSUB)
}

gens_class_CTSUB = lapply(gens_threshold_CTSUB,class_gens_CTSUB)

gens_class_CTSUB = data.frame(matrix(unlist(gens_class_CTSUB),nrow = length(gens_class_CTSUB),
                                     byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CTSUB) = c("Gene","Express","Log2FC","P-valor")

CTSUB_up = length(gens_class_CTSUB$Express[gens_class_CTSUB$Express == "UP"])
CTSUB_down = length(gens_class_CTSUB$Express[gens_class_CTSUB$Express == "DOWN"])

# NP SUB.
class_gens_NPSUB = function(y,gens_NPSUB = data.frame(Gen = character(),
                                                      Express = character(),
                                                      FC = numeric(),
                                                      PV = numeric(),
                                                      stringsAsFactors = F)){
  ex = as.data.frame(sapply(NPSUB_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_NPSUB[nrow(gens_NPSUB) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_NPSUB[nrow(gens_NPSUB) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_NPSUB)
}

gens_class_NPSUB = lapply(gens_threshold_NPSUB,class_gens_NPSUB)

gens_class_NPSUB = data.frame(matrix(unlist(gens_class_NPSUB),nrow = length(gens_class_NPSUB),
                                     byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_NPSUB) = c("Gene","Express","Log2FC","P-valor")

NPSUB_up = length(gens_class_NPSUB$Express[gens_class_NPSUB$Express == "UP"])
NPSUB_down = length(gens_class_NPSUB$Express[gens_class_NPSUB$Express == "DOWN"])

# SUB-ProS.
class_gens_SUBPROS = function(y,gens_SUBPROS = data.frame(Gen = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  ex = as.data.frame(sapply(SUBPROS_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  g_mean = mean(df)
  mitpv = mean(pv)
  if (g_mean > 0) {
    gens_SUBPROS[nrow(gens_SUBPROS) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_SUBPROS[nrow(gens_SUBPROS) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_SUBPROS)
}

gens_class_SUBPROS = lapply(gens_threshold_SUBPROS,class_gens_SUBPROS)

gens_class_SUBPROS = data.frame(matrix(unlist(gens_class_SUBPROS),
                                       nrow = length(gens_class_SUBPROS),byrow = TRUE),
                                stringsAsFactors = F)
colnames(gens_class_SUBPROS) = c("Gene","Express","Log2FC","P-valor")

SUBPROS_up = length(gens_class_SUBPROS$Express[gens_class_SUBPROS$Express == "UP"])
SUBPROS_down = length(gens_class_SUBPROS$Express[gens_class_SUBPROS$Express == "DOWN"])


rc_totals = c(length(sig_gens_CA1),length(sig_gens_CA1PROS),length(sig_gens_CA2),
              length(sig_gens_CA3),length(sig_gens_DG),length(sig_gens_CTSUB),
              length(sig_gens_NPSUB),length(sig_gens_SUBPROS))
rc_threshold = c(length(gens_threshold_CA1),length(gens_threshold_CA1PROS),length(gens_threshold_CA2),
               length(gens_threshold_CA3),length(gens_threshold_DG),length(gens_threshold_CTSUB),
               length(gens_threshold_NPSUB),length(gens_threshold_SUBPROS))
rc_up = c(CA1_up,CA1PROS_up,CA2_up,CA3_up,DG_up,CTSUB_up,NPSUB_up,SUBPROS_up)
rc_down = c(CA1_down,CA1PROS_down,CA2_down,CA3_down,DG_down,CTSUB_down,NPSUB_down,SUBPROS_down)

summary_table = data.frame(colnames(genecount_df),rc_totals,rc_threshold,rc_up,rc_down)
colnames(summary_table) = c("Class","Significant Genes","Threshold Genes","Threshold Genes UP",
                          "Threshold Genes DOWN")

thres_genes_up = unique(c(gens_class_CA1$Gen[gens_class_CA1$Express == "UP"],
                        gens_class_CA1PROS$Gen[gens_class_CA1PROS$Express == "UP"],
                        gens_class_CA2$Gen[gens_class_CA2$Express == "UP"],
                        gens_class_CA3$Gen[gens_class_CA3$Express == "UP"],
                        gens_class_DG$Gen[gens_class_DG$Express == "UP"],
                        gens_class_CTSUB$Gen[gens_class_CTSUB$Express == "UP"],
                        gens_class_NPSUB$Gen[gens_class_NPSUB$Express == "UP"],
                        gens_class_SUBPROS$Gen[gens_class_SUBPROS$Express == "UP"]))
thres_genes_down = unique(c(gens_class_CA1$Gen[gens_class_CA1$Express == "DOWN"],
                          gens_class_CA1PROS$Gen[gens_class_CA1PROS$Express == "DOWN"],
                          gens_class_CA2$Gen[gens_class_CA2$Express == "DOWN"],
                          gens_class_CA3$Gen[gens_class_CA3$Express == "DOWN"],
                          gens_class_DG$Gen[gens_class_DG$Express == "DOWN"],
                          gens_class_CTSUB$Gen[gens_class_CTSUB$Express == "DOWN"],
                          gens_class_NPSUB$Gen[gens_class_NPSUB$Express == "DOWN"],
                          gens_class_SUBPROS$Gen[gens_class_SUBPROS$Express == "DOWN"]))

png(filename = "CA1_Pattern.png",width = 800,height = 600)
pattern_CA1 = hist(genecount_df$CA1,breaks = 10,main = "Expression Pattern",ylim = c(0,1000),
                 xlab = colnames(genecount_df)[1],ylab = "Frequency",col = "lightblue")
dev.off()
pattern_t_CA1 = transform(table(cut(genecount_df$CA1,breaks = 10)))

png(filename = "CA1PROS_Pattern.png",width = 800,height = 600)
pattern_CA1PROS = hist(genecount_df$`CA1-ProS`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[2],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_CA1PROS = transform(table(cut(genecount_df$`CA1-ProS`,breaks = 10)))

png(filename = "CA2_Pattern.png",width = 800,height = 600)
pattern_CA2 = hist(genecount_df$CA2,breaks = 10,main = "Expression Pattern",ylim = c(0,1000),
                 xlab = colnames(genecount_df)[3],ylab = "Frequency",col = "lightblue")
dev.off()
pattern_t_CA2 = transform(table(cut(genecount_df$CA2,breaks = 10)))

png(filename = "CA3_Pattern.png",width = 800,height = 600)
pattern_CA3 = hist(genecount_df$CA3,breaks = 10,main = "Expression Pattern",ylim = c(0,1000),
                 xlab = colnames(genecount_df)[4],ylab = "Frequency",col = "lightblue")
dev.off()
pattern_t_CA3 = transform(table(cut(genecount_df$CA3,breaks = 10)))

png(filename = "DG_Pattern.png",width = 800,height = 600)
pattern_DG = hist(genecount_df$DG,breaks = 10,main = "Expression Pattern",ylim = c(0,1000),
                xlab = colnames(genecount_df)[5],ylab = "Frequency",col = "lightblue")
dev.off()
pattern_t_DG = transform(table(cut(genecount_df$DG,breaks = 10)))

png(filename = "CTSUB_Pattern.png",width = 800,height = 600)
pattern_CTSUB = hist(genecount_df$`CT SUB`,breaks = 10,main = "Expression Pattern",ylim = c(0,1000),
                   xlab = colnames(genecount_df)[6],ylab = "Frequency",col = "lightblue")
dev.off()
pattern_t_CTSUB = transform(table(cut(genecount_df$`CT SUB`,breaks = 10)))

png(filename = "NPSUB_Pattern.png",width = 800,height = 600)
pattern_NPSUB = hist(genecount_df$`NP SUB`,breaks = 10,main = "Expression Pattern",ylim = c(0,1000),
                   xlab = colnames(genecount_df)[7],ylab = "Frequency",col = "lightblue")
dev.off()
pattern_t_NPSUB = transform(table(cut(genecount_df$`NP SUB`,breaks = 10)))

png(filename = "SUBPROS_Pattern.png",width = 800,height = 600)
pattern_SUBPROS = hist(genecount_df$`SUB-ProS`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[8],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_SUBPROS = transform(table(cut(genecount_df$`SUB-ProS`,breaks = 10)))


UP_CA1 = gens_class_CA1$Gen[gens_class_CA1$Express == "UP"]
DOWN_CA1 = gens_class_CA1$Gen[gens_class_CA1$Express == "DOWN"]
UP_CA1PROS = gens_class_CA1PROS$Gen[gens_class_CA1PROS$Express == "UP"]
DOWN_CA1PROS = gens_class_CA1PROS$Gen[gens_class_CA1PROS$Express == "DOWN"]
UP_CA2 = gens_class_CA2$Gen[gens_class_CA2$Express == "UP"]
DOWN_CA2 = gens_class_CA2$Gen[gens_class_CA2$Express == "DOWN"]
UP_CA3 = gens_class_CA3$Gen[gens_class_CA3$Express == "UP"]
DOWN_CA3 = gens_class_CA3$Gen[gens_class_CA3$Express == "DOWN"]
UP_DG = gens_class_DG$Gen[gens_class_DG$Express == "UP"]
DOWN_DG = gens_class_DG$Gen[gens_class_DG$Express == "DOWN"]
UP_CTSUB = gens_class_CTSUB$Gen[gens_class_CTSUB$Express == "UP"]
DOWN_CTSUB = gens_class_CTSUB$Gen[gens_class_CTSUB$Express == "DOWN"]
UP_NPSUB = gens_class_NPSUB$Gen[gens_class_NPSUB$Express == "UP"]
DOWN_NPSUB = gens_class_NPSUB$Gen[gens_class_NPSUB$Express == "DOWN"]
UP_SUBPROS = gens_class_SUBPROS$Gen[gens_class_SUBPROS$Express == "UP"]
DOWN_SUBPROS = gens_class_SUBPROS$Gen[gens_class_SUBPROS$Express == "DOWN"]

gcount_sig = data.frame(genecount_df >= threshold)
gcount_sig = gcount_sig[rowSums(gcount_sig) != 0,]
colnames(gcount_sig) = c("CA1","CA1-ProS","CA2","CA3","DG","CT SUB","NP SUB","SUB-ProS")

gcount_sig$CA1[rownames(gcount_sig) %in% UP_CA1] = "UP"
gcount_sig$CA1[rownames(gcount_sig) %in% DOWN_CA1] = "DOWN"
gcount_sig$`CA1-ProS`[rownames(gcount_sig) %in% UP_CA1PROS] = "UP"
gcount_sig$`CA1-ProS`[rownames(gcount_sig) %in% DOWN_CA1PROS] = "DOWN"
gcount_sig$CA2[rownames(gcount_sig) %in% UP_CA2] = "UP"
gcount_sig$CA2[rownames(gcount_sig) %in% DOWN_CA2] = "DOWN"
gcount_sig$CA3[rownames(gcount_sig) %in% UP_CA3] = "UP"
gcount_sig$CA3[rownames(gcount_sig) %in% DOWN_CA3] = "DOWN"
gcount_sig$DG[rownames(gcount_sig) %in% UP_DG] = "UP"
gcount_sig$DG[rownames(gcount_sig) %in% DOWN_DG] = "DOWN"
gcount_sig$`CT SUB`[rownames(gcount_sig) %in% UP_CTSUB] = "UP"
gcount_sig$`CT SUB`[rownames(gcount_sig) %in% DOWN_CTSUB] = "DOWN"
gcount_sig$`NP SUB`[rownames(gcount_sig) %in% UP_NPSUB] = "UP"
gcount_sig$`NP SUB`[rownames(gcount_sig) %in% DOWN_NPSUB] = "DOWN"
gcount_sig$`SUB-ProS`[rownames(gcount_sig) %in% UP_SUBPROS] = "UP"
gcount_sig$`SUB-ProS`[rownames(gcount_sig) %in% DOWN_SUBPROS] = "DOWN"

UP = unname(apply(gcount_sig,1,function (x) paste(names(x[x == "UP"]), collapse = ",")))
count_UP = apply(gcount_sig,1,function (x) length(x[x == "UP"]))
coinc_UP = data.frame(count_UP,UP)
colnames(coinc_UP) = c("Number of Classes","Classes")
coinc_UP = coinc_UP[coinc_UP$`Number of Classes` !=0,]

DOWN = unname(apply(gcount_sig,1,function (x) paste(names(x[x == "DOWN"]), collapse = ",")))
count_DOWN = apply(gcount_sig,1,function (x) length(x[x == "DOWN"]))
coinc_DOWN = data.frame(count_DOWN,DOWN)
colnames(coinc_DOWN) = c("Number of Classes","Classes")
coinc_DOWN = coinc_DOWN[coinc_DOWN$`Number of Classes` !=0,]

png("Coinc_UP_2DEG.png",width = 800,height = 600)
par(mar = c(10,4,5,1))
barplot(table(coinc_UP[coinc_UP$`Number of Classes` == 2,length(coinc_UP)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()

png("Coinc_UP_3DEG.png",width = 800,height = 600)
par(mar = c(12.5,4,5,1))
barplot(table(coinc_UP[coinc_UP$`Number of Classes` == 3,length(coinc_UP)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()

png("Coinc_DOWN_2DEG.png",width = 800,height = 600)
par(mar = c(10,4,5,1))
barplot(table(coinc_DOWN[coinc_DOWN$`Number of Classes` == 2,length(coinc_DOWN)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()

png("Coinc_DOWN_3DEG.png",width = 800,height = 600)
par(mar = c(12.5,4,5,1))
barplot(table(coinc_DOWN[coinc_DOWN$`Number of Classes` == 3,length(coinc_DOWN)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()


pattern = c("CA1$","CA1,")
png("Coinc_UP_CA1.png")
barplot(table(coinc_UP$`Number of Classes`[grep(paste(pattern,collapse = "|"),coinc_UP$Classes)]),
        ylab = "Number of genes",
        xlab = "Num. of classes with which the DEG UP is shared (1 = CA1)",col = "lightgreen")
dev.off()

coinc_UP_CA1 = coinc_UP[grep(paste(pattern,collapse = "|"),coinc_UP$Classes),]
coinc_UP_CA1_unique = rownames(coinc_UP_CA1[coinc_UP_CA1$`Number of Classes` == 1,])
coinc_UP_CA1_d = rownames(coinc_UP_CA1[coinc_UP_CA1$`Number of Classes` == 2,])
coinc_UP_CA1_t = rownames(coinc_UP_CA1[coinc_UP_CA1$`Number of Classes` == 3,])
coinc_UP_CA1_q = rownames(coinc_UP_CA1[coinc_UP_CA1$`Number of Classes` == 4,])

pattern = c("CA1$","CA1,")
png("Coinc_DOWN_CA1.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep(paste(pattern,collapse = "|"),coinc_DOWN$Classes)]),
        ylab = "Number of genes",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = CA1)",col = "lightgreen")
dev.off()

coinc_DOWN_CA1 = coinc_DOWN[grep(paste(pattern,collapse = "|"),coinc_DOWN$Classes),]
coinc_DOWN_CA1_unique = rownames(coinc_DOWN_CA1[coinc_DOWN_CA1$`Number of Classes` == 1,])
coinc_DOWN_CA1_d = rownames(coinc_DOWN_CA1[coinc_DOWN_CA1$`Number of Classes` == 2,])
coinc_DOWN_CA1_t = rownames(coinc_DOWN_CA1[coinc_DOWN_CA1$`Number of Classes` == 3,])
coinc_DOWN_CA1_q = rownames(coinc_DOWN_CA1[coinc_DOWN_CA1$`Number of Classes` == 4,])
coinc_DOWN_CA1_c = rownames(coinc_DOWN_CA1[coinc_DOWN_CA1$`Number of Classes` == 5,])

png("Coinc_UP_CA1PROS.png")
barplot(table(coinc_UP$`Number of Classes`[grep("CA1-ProS",coinc_UP$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG UP is shared (1 = CA1-ProS)")
dev.off()

coinc_UP_CA1PROS = coinc_UP[grep("CA1-ProS",coinc_UP$Classes),]
coinc_UP_CA1PROS_unique = rownames(coinc_UP_CA1PROS[coinc_UP_CA1PROS$`Number of Classes` == 1,])
coinc_UP_CA1PROS_d = rownames(coinc_UP_CA1PROS[coinc_UP_CA1PROS$`Number of Classes` == 2,])
coinc_UP_CA1PROS_t = rownames(coinc_UP_CA1PROS[coinc_UP_CA1PROS$`Number of Classes` == 3,])
coinc_UP_CA1PROS_q = rownames(coinc_UP_CA1PROS[coinc_UP_CA1PROS$`Number of Classes` == 4,])

png("Coinc_DOWN_CA1PROS.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep("CA1-ProS",coinc_DOWN$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = CA1-ProS)")
dev.off()

coinc_DOWN_CA1PROS = coinc_DOWN[grep("CA1-ProS",coinc_DOWN$Classes),]
coinc_DOWN_CA1PROS_unique = rownames(coinc_DOWN_CA1PROS[coinc_DOWN_CA1PROS$`Number of Classes` == 1,])
coinc_DOWN_CA1PROS_d = rownames(coinc_DOWN_CA1PROS[coinc_DOWN_CA1PROS$`Number of Classes` == 2,])
coinc_DOWN_CA1PROS_t = rownames(coinc_DOWN_CA1PROS[coinc_DOWN_CA1PROS$`Number of Classes` == 3,])
coinc_DOWN_CA1PROS_q = rownames(coinc_DOWN_CA1PROS[coinc_DOWN_CA1PROS$`Number of Classes` == 4,])
coinc_DOWN_CA1PROS_c = rownames(coinc_DOWN_CA1PROS[coinc_DOWN_CA1PROS$`Number of Classes` == 5,])

png("Coinc_UP_CA2.png")
barplot(table(coinc_UP$`Number of Classes`[grep("CA2",coinc_UP$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG UP is shared (1 = CA2)")
dev.off()

coinc_UP_CA2 = coinc_UP[grep("CA2",coinc_UP$Classes),]
coinc_UP_CA2_unique = rownames(coinc_UP_CA2[coinc_UP_CA2$`Number of Classes` == 1,])
coinc_UP_CA2_d = rownames(coinc_UP_CA2[coinc_UP_CA2$`Number of Classes` == 2,])
coinc_UP_CA2_t = rownames(coinc_UP_CA2[coinc_UP_CA2$`Number of Classes` == 3,])
coinc_UP_CA2_q = rownames(coinc_UP_CA2[coinc_UP_CA2$`Number of Classes` == 4,])

png("Coinc_DOWN_CA2.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep("CA2",coinc_DOWN$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = CA2)")
dev.off()

coinc_DOWN_CA2 = coinc_DOWN[grep("CA2",coinc_DOWN$Classes),]
coinc_DOWN_CA2_unique = rownames(coinc_DOWN_CA2[coinc_DOWN_CA2$`Number of Classes` == 1,])
coinc_DOWN_CA2_d = rownames(coinc_DOWN_CA2[coinc_DOWN_CA2$`Number of Classes` == 2,])
coinc_DOWN_CA2_t = rownames(coinc_DOWN_CA2[coinc_DOWN_CA2$`Number of Classes` == 3,])
coinc_DOWN_CA2_q = rownames(coinc_DOWN_CA2[coinc_DOWN_CA2$`Number of Classes` == 4,])
coinc_DOWN_CA2_c = rownames(coinc_DOWN_CA2[coinc_DOWN_CA2$`Number of Classes` == 5,])

png("Coinc_UP_CA3.png")
barplot(table(coinc_UP$`Number of Classes`[grep("CA3",coinc_UP$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG UP is shared (1 = CA3)")
dev.off()

coinc_UP_CA3 = coinc_UP[grep("CA3",coinc_UP$Classes),]
coinc_UP_CA3_unique = rownames(coinc_UP_CA3[coinc_UP_CA3$`Number of Classes` == 1,])
coinc_UP_CA3_d = rownames(coinc_UP_CA3[coinc_UP_CA3$`Number of Classes` == 2,])
coinc_UP_CA3_t = rownames(coinc_UP_CA3[coinc_UP_CA3$`Number of Classes` == 3,])
coinc_UP_CA3_q = rownames(coinc_UP_CA3[coinc_UP_CA3$`Number of Classes` == 3,])

png("Coinc_DOWN_CA3.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep("CA3",coinc_DOWN$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = CA3)")
dev.off()

coinc_DOWN_CA3 = coinc_DOWN[grep("CA3",coinc_DOWN$Classes),]
coinc_DOWN_CA3_unique = rownames(coinc_DOWN_CA3[coinc_DOWN_CA3$`Number of Classes` == 1,])
coinc_DOWN_CA3_d = rownames(coinc_DOWN_CA3[coinc_DOWN_CA3$`Number of Classes` == 2,])
coinc_DOWN_CA3_t = rownames(coinc_DOWN_CA3[coinc_DOWN_CA3$`Number of Classes` == 3,])
coinc_DOWN_CA3_q = rownames(coinc_DOWN_CA3[coinc_DOWN_CA3$`Number of Classes` == 4,])
coinc_DOWN_CA3_c = rownames(coinc_DOWN_CA3[coinc_DOWN_CA3$`Number of Classes` == 5,])

png("Coinc_UP_DG.png")
barplot(table(coinc_UP$`Number of Classes`[grep("DG",coinc_UP$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG UP is shared (1 = DG)")
dev.off()

coinc_UP_DG = coinc_UP[grep("DG",coinc_UP$Classes),]
coinc_UP_DG_unique = rownames(coinc_UP_DG[coinc_UP_DG$`Number of Classes` == 1,])
coinc_UP_DG_d = rownames(coinc_UP_DG[coinc_UP_DG$`Number of Classes` == 2,])
coinc_UP_DG_t = rownames(coinc_UP_DG[coinc_UP_DG$`Number of Classes` == 3,])
coinc_UP_DG_q = rownames(coinc_UP_DG[coinc_UP_DG$`Number of Classes` == 4,])

png("Coinc_DOWN_DG.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep("DG",coinc_DOWN$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = DG)")
dev.off()

coinc_DOWN_DG = coinc_DOWN[grep("DG",coinc_DOWN$Classes),]
coinc_DOWN_DG_unique = rownames(coinc_DOWN_DG[coinc_DOWN_DG$`Number of Classes` == 1,])
coinc_DOWN_DG_d = rownames(coinc_DOWN_DG[coinc_DOWN_DG$`Number of Classes` == 2,])
coinc_DOWN_DG_t = rownames(coinc_DOWN_DG[coinc_DOWN_DG$`Number of Classes` == 3,])
coinc_DOWN_DG_q = rownames(coinc_DOWN_DG[coinc_DOWN_DG$`Number of Classes` == 4,])
coinc_DOWN_DG_c = rownames(coinc_DOWN_DG[coinc_DOWN_DG$`Number of Classes` == 5,])

png("Coinc_UP_CTSUB.png")
barplot(table(coinc_UP$`Number of Classes`[grep("CT SUB",coinc_UP$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG UP is shared (1 = CT SUB)")
dev.off()

coinc_UP_CTSUB = coinc_UP[grep("CT SUB",coinc_UP$Classes),]
coinc_UP_CTSUB_unique = rownames(coinc_UP_CTSUB[coinc_UP_CTSUB$`Number of Classes` == 1,])
coinc_UP_CTSUB_d = rownames(coinc_UP_CTSUB[coinc_UP_CTSUB$`Number of Classes` == 2,])
coinc_UP_CTSUB_t = rownames(coinc_UP_CTSUB[coinc_UP_CTSUB$`Number of Classes` == 3,])
coinc_UP_CTSUB_q = rownames(coinc_UP_CTSUB[coinc_UP_CTSUB$`Number of Classes` == 4,])

png("Coinc_DOWN_CTSUB.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep("CT SUB",coinc_DOWN$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = CT SUB)")
dev.off()

coinc_DOWN_CTSUB = coinc_DOWN[grep("CT SUB",coinc_DOWN$Classes),]
coinc_DOWN_CTSUB_unique = rownames(coinc_DOWN_CTSUB[coinc_DOWN_CTSUB$`Number of Classes` == 1,])
coinc_DOWN_CTSUB_d = rownames(coinc_DOWN_CTSUB[coinc_DOWN_CTSUB$`Number of Classes` == 2,])
coinc_DOWN_CTSUB_t = rownames(coinc_DOWN_CTSUB[coinc_DOWN_CTSUB$`Number of Classes` == 3,])
coinc_DOWN_CTSUB_q = rownames(coinc_DOWN_CTSUB[coinc_DOWN_CTSUB$`Number of Classes` == 4,])
coinc_DOWN_CTSUB_c = rownames(coinc_DOWN_CTSUB[coinc_DOWN_CTSUB$`Number of Classes` == 5,])

png("Coinc_UP_NPSUB.png")
barplot(table(coinc_UP$`Number of Classes`[grep("NP SUB",coinc_UP$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG UP is shared (1 = NP SUB)")
dev.off()

coinc_UP_NPSUB = coinc_UP[grep("NP SUB",coinc_UP$Classes),]
coinc_UP_NPSUB_unique = rownames(coinc_UP_NPSUB[coinc_UP_NPSUB$`Number of Classes` == 1,])
coinc_UP_NPSUB_d = rownames(coinc_UP_NPSUB[coinc_UP_NPSUB$`Number of Classes` == 2,])
coinc_UP_NPSUB_t = rownames(coinc_UP_NPSUB[coinc_UP_NPSUB$`Number of Classes` == 3,])
coinc_UP_NPSUB_q = rownames(coinc_UP_NPSUB[coinc_UP_NPSUB$`Number of Classes` == 4,])

png("Coinc_DOWN_NPSUB.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep("NP SUB",coinc_DOWN$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = NP SUB)")
dev.off()

coinc_DOWN_NPSUB = coinc_DOWN[grep("NP SUB",coinc_DOWN$Classes),]
coinc_DOWN_NPSUB_unique = rownames(coinc_DOWN_NPSUB[coinc_DOWN_NPSUB$`Number of Classes` == 1,])
coinc_DOWN_NPSUB_d = rownames(coinc_DOWN_NPSUB[coinc_DOWN_NPSUB$`Number of Classes` == 2,])
coinc_DOWN_NPSUB_t = rownames(coinc_DOWN_NPSUB[coinc_DOWN_NPSUB$`Number of Classes` == 3,])
coinc_DOWN_NPSUB_q = rownames(coinc_DOWN_NPSUB[coinc_DOWN_NPSUB$`Number of Classes` == 4,])
coinc_DOWN_NPSUB_c = rownames(coinc_DOWN_NPSUB[coinc_DOWN_NPSUB$`Number of Classes` == 5,])

png("Coinc_UP_SUBPROS.png")
barplot(table(coinc_UP$`Number of Classes`[grep("SUB-ProS",coinc_UP$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG UP is shared (1 = SUB-ProS)")
dev.off()

coinc_UP_SUBPROS = coinc_UP[grep("SUB-ProS",coinc_UP$Classes),]
coinc_UP_SUBPROS_unique = rownames(coinc_UP_SUBPROS[coinc_UP_SUBPROS$`Number of Classes` == 1,])
coinc_UP_SUBPROS_d = rownames(coinc_UP_SUBPROS[coinc_UP_SUBPROS$`Number of Classes` == 2,])
coinc_UP_SUBPROS_t = rownames(coinc_UP_SUBPROS[coinc_UP_SUBPROS$`Number of Classes` == 3,])
coinc_UP_SUBPROS_q = rownames(coinc_UP_SUBPROS[coinc_UP_SUBPROS$`Number of Classes` == 4,])

png("Coinc_DOWN_SUBPROS.png")
barplot(table(coinc_DOWN$`Number of Classes`[grep("SUB-ProS",coinc_DOWN$Classes)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of Classes with which the DEG DOWN is shared (1 = SUB-ProS)")
dev.off()

coinc_DOWN_SUBPROS = coinc_DOWN[grep("SUB-ProS",coinc_DOWN$Classes),]
coinc_DOWN_SUBPROS_unique = rownames(coinc_DOWN_SUBPROS[coinc_DOWN_SUBPROS$`Number of Classes` == 1,])
coinc_DOWN_SUBPROS_d = rownames(coinc_DOWN_SUBPROS[coinc_DOWN_SUBPROS$`Number of Classes` == 2,])
coinc_DOWN_SUBPROS_t = rownames(coinc_DOWN_SUBPROS[coinc_DOWN_SUBPROS$`Number of Classes` == 3,])
coinc_DOWN_SUBPROS_q = rownames(coinc_DOWN_SUBPROS[coinc_DOWN_SUBPROS$`Number of Classes` == 4,])
coinc_DOWN_SUBPROS_c = rownames(coinc_DOWN_SUBPROS[coinc_DOWN_SUBPROS$`Number of Classes` == 5,])

unique_up = c(length(coinc_UP_CA1_unique),length(coinc_UP_CA1PROS_unique),length(coinc_UP_CA2_unique),
             length(coinc_UP_CA3_unique),length(coinc_UP_DG_unique),length(coinc_UP_CTSUB_unique),
             length(coinc_UP_NPSUB_unique),length(coinc_UP_SUBPROS_unique))
two_up = c(length(coinc_UP_CA1_d),length(coinc_UP_CA1PROS_d),length(coinc_UP_CA2_d),
           length(coinc_UP_CA3_d),length(coinc_UP_DG_d),length(coinc_UP_CTSUB_d),
           length(coinc_UP_NPSUB_d),length(coinc_UP_SUBPROS_d))
three_up = c(length(coinc_UP_CA1_t),length(coinc_UP_CA1PROS_t),length(coinc_UP_CA2_t),
           length(coinc_UP_CA3_t),length(coinc_UP_DG_t),length(coinc_UP_CTSUB_t),
           length(coinc_UP_NPSUB_t),length(coinc_UP_SUBPROS_t))
four_up = c(length(coinc_UP_CA1_q),length(coinc_UP_CA1PROS_q),length(coinc_UP_CA2_q),
               length(coinc_UP_CA3_q),length(coinc_UP_DG_q),length(coinc_UP_CTSUB_q),
               length(coinc_UP_NPSUB_q),length(coinc_UP_SUBPROS_q))
tsummary_up = data.frame(unique_up,two_up,three_up,four_up)
colnames(tsummary_up) = c("Unique","Two Classes","Three Classes","Four Classes")
rownames(tsummary_up) = c("CA1","CA1-ProS","CA2","CA3","DG","CT SUB","NP SUB","SUB-ProS")

unique_down = c(length(coinc_DOWN_CA1_unique),length(coinc_DOWN_CA1PROS_unique),
               length(coinc_DOWN_CA2_unique),length(coinc_DOWN_CA3_unique),length(coinc_DOWN_DG_unique),
               length(coinc_DOWN_CTSUB_unique),length(coinc_DOWN_NPSUB_unique),
               length(coinc_DOWN_SUBPROS_unique))
two_down = c(length(coinc_DOWN_CA1_d),length(coinc_DOWN_CA1PROS_d),length(coinc_DOWN_CA2_d),
           length(coinc_DOWN_CA3_d),length(coinc_DOWN_DG_d),length(coinc_DOWN_CTSUB_d),
           length(coinc_DOWN_NPSUB_d),length(coinc_DOWN_SUBPROS_d))
three_down = c(length(coinc_DOWN_CA1_t),length(coinc_DOWN_CA1PROS_t),length(coinc_DOWN_CA2_t),
             length(coinc_DOWN_CA3_t),length(coinc_DOWN_DG_t),length(coinc_DOWN_CTSUB_t),
             length(coinc_DOWN_NPSUB_t),length(coinc_DOWN_SUBPROS_t))
four_down = c(length(coinc_DOWN_CA1_q),length(coinc_DOWN_CA1PROS_q),length(coinc_DOWN_CA2_q),
               length(coinc_DOWN_CA3_q),length(coinc_DOWN_DG_q),length(coinc_DOWN_CTSUB_q),
               length(coinc_DOWN_NPSUB_q),length(coinc_DOWN_SUBPROS_q))
five_down = c(length(coinc_DOWN_CA1_c),length(coinc_DOWN_CA1PROS_c),length(coinc_DOWN_CA2_c),
               length(coinc_DOWN_CA3_c),length(coinc_DOWN_DG_c),length(coinc_DOWN_CTSUB_c),
               length(coinc_DOWN_NPSUB_c),length(coinc_DOWN_SUBPROS_c))
tsummary_down = data.frame(uniques_down,dos_down,tres_down,quatre_down,five_down)
colnames(tsummary_down) = c("Unique","Two Classes","Three Classes","Four Classes","Five Classes")
rownames(tsummary_down) = c("CA1","CA1-ProS","CA2","CA3","DG","CT SUB","NP SUB","SUB-ProS")



write.csv(x = genecount_df,file = "DE_Genes_Class.csv")

write.csv(gens_class_CA1,"List_Genes_Sig_CA1.csv")
write.csv(gens_class_CA1PROS,"List_Genes_Sig_CA1PROS.csv")
write.csv(gens_class_CA2,"List_Genes_Sig_CA2.csv")
write.csv(gens_class_CA3,"List_Genes_Sig_CA3.csv")
write.csv(gens_class_DG,"List_Genes_Sig_DG.csv")
write.csv(gens_class_CTSUB,"List_Genes_Sig_CTSUB.csv")
write.csv(gens_class_NPSUB,"List_Genes_Sig_NPSUB.csv")
write.csv(gens_class_SUBPROS,"List_Genes_Sig_SUBPROS.csv")

write.csv(summary_table,"summary_table.csv")

write.csv(pattern_t_CA1,"pattern_Table_CA1.csv")
write.csv(pattern_t_CA1PROS,"pattern_Table_CA1PROS.csv")
write.csv(pattern_t_CA2,"pattern_Table_CA2.csv")
write.csv(pattern_t_CA3,"pattern_Table_CA3.csv")
write.csv(pattern_t_DG,"pattern_Table_DG.csv")
write.csv(pattern_t_CTSUB,"pattern_Table_CTSUB.csv")
write.csv(pattern_t_NPSUB,"pattern_Table_NPSUB.csv")
write.csv(pattern_t_SUBPROS,"pattern_Table_SUBPROS.csv")

write.csv(thres_genes_up,"Threshold_Genes_Total_UP.csv")
write.csv(thres_genes_down,"Threshold_Genes_Total_DOWN.csv")

write.csv(coinc_UP_CA1_unique,"DEGsUP_uniques_CA1.csv")
write.csv(coinc_UP_CA1_d,"DEGsUP_Shared2_CA1.csv")
write.csv(coinc_UP_CA1_t,"DEGsUP_Shared3_CA1.csv")
write.csv(coinc_UP_CA1PROS_unique,"DEGsUP_uniques_CA1PROS.csv")
write.csv(coinc_UP_CA1PROS_d,"DEGsUP_Shared2_CA1PROS.csv")
write.csv(coinc_UP_CA1PROS_t,"DEGsUP_Shared3_CA1PROS.csv")
write.csv(coinc_UP_CA2_unique,"DEGsUP_uniques_CA2.csv")
write.csv(coinc_UP_CA2_d,"DEGsUP_Shared2_CA2.csv")
write.csv(coinc_UP_CA2_t,"DEGsUP_Shared3_CA2.csv")
write.csv(coinc_UP_CA3_unique,"DEGsUP_uniques_CA3.csv")
write.csv(coinc_UP_CA3_d,"DEGsUP_Shared2_CA3.csv")
write.csv(coinc_UP_CA3_t,"DEGsUP_Shared3_CA3.csv")
write.csv(coinc_UP_DG_unique,"DEGsUP_uniques_DG.csv")
write.csv(coinc_UP_DG_d,"DEGsUP_Shared2_DG.csv")
write.csv(coinc_UP_DG_t,"DEGsUP_Shared3_DG.csv")
write.csv(coinc_UP_CTSUB_unique,"DEGsUP_uniques_CTSUB.csv")
write.csv(coinc_UP_CTSUB_d,"DEGsUP_Shared2_CTSUB.csv")
write.csv(coinc_UP_CTSUB_t,"DEGsUP_Shared3_CTSUB.csv")
write.csv(coinc_UP_NPSUB_unique,"DEGsUP_uniques_NPSUB.csv")
write.csv(coinc_UP_NPSUB_d,"DEGsUP_Shared2_NPSUB.csv")
write.csv(coinc_UP_NPSUB_t,"DEGsUP_Shared3_NPSUB.csv")
write.csv(coinc_UP_SUBPROS_unique,"DEGsUP_uniques_SUBPROS.csv")
write.csv(coinc_UP_SUBPROS_d,"DEGsUP_Shared2_SUBPROS.csv")
write.csv(coinc_UP_SUBPROS_t,"DEGsUP_Shared3_SUBPROS.csv")

write.csv(coinc_DOWN_CA1_unique,"DEGsDOWN_uniques_CA1.csv")
write.csv(coinc_DOWN_CA1_d,"DEGsDOWN_Shared2_CA1.csv")
write.csv(coinc_DOWN_CA1_t,"DEGsDOWN_Shared3_CA1.csv")
write.csv(coinc_DOWN_CA1PROS_unique,"DEGsDOWN_uniques_CA1PROS.csv")
write.csv(coinc_DOWN_CA1PROS_d,"DEGsDOWN_Shared2_CA1PROS.csv")
write.csv(coinc_DOWN_CA1PROS_t,"DEGsDOWN_Shared3_CA1PROS.csv")
write.csv(coinc_DOWN_CA2_unique,"DEGsDOWN_uniques_CA2.csv")
write.csv(coinc_DOWN_CA2_d,"DEGsDOWN_Shared2_CA2.csv")
write.csv(coinc_DOWN_CA2_t,"DEGsDOWN_Shared3_CA2.csv")
write.csv(coinc_DOWN_CA3_unique,"DEGsDOWN_uniques_CA3.csv")
write.csv(coinc_DOWN_CA3_d,"DEGsDOWN_Shared2_CA3.csv")
write.csv(coinc_DOWN_CA3_t,"DEGsDOWN_Shared3_CA3.csv")
write.csv(coinc_DOWN_DG_unique,"DEGsDOWN_uniques_DG.csv")
write.csv(coinc_DOWN_DG_d,"DEGsDOWN_Shared2_DG.csv")
write.csv(coinc_DOWN_DG_t,"DEGsDOWN_Shared3_DG.csv")
write.csv(coinc_DOWN_CTSUB_unique,"DEGsDOWN_uniques_CTSUB.csv")
write.csv(coinc_DOWN_CTSUB_d,"DEGsDOWN_Shared2_CTSUB.csv")
write.csv(coinc_DOWN_CTSUB_t,"DEGsDOWN_Shared3_CTSUB.csv")
write.csv(coinc_DOWN_NPSUB_unique,"DEGsDOWN_uniques_NPSUB.csv")
write.csv(coinc_DOWN_NPSUB_d,"DEGsDOWN_Shared2_NPSUB.csv")
write.csv(coinc_DOWN_NPSUB_t,"DEGsDOWN_Shared3_NPSUB.csv")
write.csv(coinc_DOWN_SUBPROS_unique,"DEGsDOWN_uniques_SUBPROS.csv")
write.csv(coinc_DOWN_SUBPROS_d,"DEGsDOWN_Shared2_SUBPROS.csv")
write.csv(coinc_DOWN_SUBPROS_t,"DEGsDOWN_Shared3_SUBPROS.csv")

write.csv(tsummary_up,"Summary_table_Coincidence_UP.csv")
write.csv(tsummary_down,"Summary_table_Coincidence_DOWN.csv")
