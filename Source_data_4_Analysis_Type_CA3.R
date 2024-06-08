##File to analyze differentially expressed genes between neuronal types of the CA3 Class.
library(dplyr)
library(stringr)

results_list = readRDS("Output_Types.RDS")
results_list = unlist(results_list,recursive = FALSE)

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
  gens_filtre = filter(x,(avg_log2FC > 0.6) | (avg_log2FC < -0.6) & (p_val_bh < 0.05))
  
  final_gene_set = rbind(final_gene_set,gens_filtre)
  return(final_gene_set)
}

sig_df = lapply(results_list,filtration)

n = length(sig_df)
count_df = data.frame(matrix(ncol = 10,nrow = n/10))
types = c("349_CA3","350_CA3","351_CA3","352_CA3","353_CA3","354_CA3","355_CA3",
         "356_CA3","357_CA3","358_CA3")
colnames(count_df) = types

CA3_349_sig = sig_df[grep("^349_Mossy-",names(sig_df))]
CA3_350_sig = sig_df[grep("^350_Mossy-",names(sig_df))]
CA3_351_sig = sig_df[grep("^351_CA3-",names(sig_df))]
CA3_352_sig = sig_df[grep("^352_CA3-",names(sig_df))]
CA3_353_sig = sig_df[grep("^353_CA3-",names(sig_df))]
CA3_354_sig = sig_df[grep("^354_CA3-",names(sig_df))]
CA3_355_sig = sig_df[grep("^355_CA3-",names(sig_df))]
CA3_356_sig = sig_df[grep("^356_CA3-",names(sig_df))]
CA3_357_sig = sig_df[grep("^357_CA3-",names(sig_df))]
CA3_358_sig = sig_df[grep("^358_CA3-",names(sig_df))]

sig_gens_CA3_349 = unique(unlist(sapply(CA3_349_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_350 = unique(unlist(sapply(CA3_350_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_351 = unique(unlist(sapply(CA3_351_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_352 = unique(unlist(sapply(CA3_352_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_353 = unique(unlist(sapply(CA3_353_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_354 = unique(unlist(sapply(CA3_354_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_355 = unique(unlist(sapply(CA3_355_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_356 = unique(unlist(sapply(CA3_356_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_357 = unique(unlist(sapply(CA3_357_sig,function(x) row.names(x))),use.names=FALSE)
sig_gens_CA3_358 = unique(unlist(sapply(CA3_358_sig,function(x) row.names(x))),use.names=FALSE)

sig_gens_total = unique(c(sig_gens_CA3_349,sig_gens_CA3_350,sig_gens_CA3_351,
                          sig_gens_CA3_352,sig_gens_CA3_353,sig_gens_CA3_354,sig_gens_CA3_355,
                          sig_gens_CA3_356,sig_gens_CA3_357,sig_gens_CA3_358))

gens_totals = unique(unlist(sapply(results_list,function(x) row.names(x))),use.names=FALSE)
gens_never_sig = gens_totals[!(gens_totals%in%sig_gens_total)]

count_df$`349_CA3` = sapply(CA3_349_sig,function(x) nrow(x))
count_df$`350_CA3` = sapply(CA3_350_sig,function(x) nrow(x))
count_df$`351_CA3` = sapply(CA3_351_sig,function(x) nrow(x))
count_df$`352_CA3` = sapply(CA3_352_sig,function(x) nrow(x))
count_df$`353_CA3` = sapply(CA3_353_sig,function(x) nrow(x))
count_df$`354_CA3` = sapply(CA3_354_sig,function(x) nrow(x))
count_df$`355_CA3` = sapply(CA3_355_sig,function(x) nrow(x))
count_df$`356_CA3` = sapply(CA3_356_sig,function(x) nrow(x))
count_df$`357_CA3` = sapply(CA3_357_sig,function(x) nrow(x))
count_df$`358_CA3` = sapply(CA3_358_sig,function(x) nrow(x))

genecount_df = data.frame(matrix(data = 0,ncol = 10,nrow = length(sig_gens_total)))
row.names(genecount_df) = sig_gens_total
colnames(genecount_df) = types

for (i in 1:length(sig_gens_total)) {
  genecount_df$`349_CA3`[i] = sum(sapply(CA3_349_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`350_CA3`[i] = sum(sapply(CA3_350_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`351_CA3`[i] = sum(sapply(CA3_351_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`352_CA3`[i] = sum(sapply(CA3_352_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`353_CA3`[i] = sum(sapply(CA3_353_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`354_CA3`[i] = sum(sapply(CA3_354_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`355_CA3`[i] = sum(sapply(CA3_355_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`356_CA3`[i] = sum(sapply(CA3_356_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`357_CA3`[i] = sum(sapply(CA3_357_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
  genecount_df$`358_CA3`[i] = sum(sapply(CA3_358_sig,
                                         function(x)  sig_gens_total[i]%in% row.names(x)))
}

iterations = nrow(count_df)
threshold = 0.9 * iterations
gens_threshold_CA3_349 = rownames(genecount_df[which(genecount_df$`349_CA3` >= threshold),])
gens_threshold_CA3_350 = rownames(genecount_df[which(genecount_df$`350_CA3` >= threshold),])
gens_threshold_CA3_351 = rownames(genecount_df[which(genecount_df$`351_CA3` >= threshold),])
gens_threshold_CA3_352 = rownames(genecount_df[which(genecount_df$`352_CA3` >= threshold),])
gens_threshold_CA3_353 = rownames(genecount_df[which(genecount_df$`353_CA3` >= threshold),])
gens_threshold_CA3_354 = rownames(genecount_df[which(genecount_df$`354_CA3` >= threshold),])
gens_threshold_CA3_355 = rownames(genecount_df[which(genecount_df$`355_CA3` >= threshold),])
gens_threshold_CA3_356 = rownames(genecount_df[which(genecount_df$`356_CA3` >= threshold),])
gens_threshold_CA3_357 = rownames(genecount_df[which(genecount_df$`357_CA3` >= threshold),])
gens_threshold_CA3_358 = rownames(genecount_df[which(genecount_df$`358_CA3` >= threshold),])

# CA3_349.
class_gens_CA3_349 = function(y,gens_CA3_349 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  
  ex = as.data.frame(sapply(CA3_349_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  
  g_mean = mean(df)
  
  mitpv = mean(pv)
  
  
  if (g_mean > 0) {
    gens_CA3_349[nrow(gens_CA3_349) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_349[nrow(gens_CA3_349) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_349)
}

gens_class_CA3_349 = lapply(gens_threshold_CA3_349,class_gens_CA3_349)

gens_class_CA3_349 = data.frame(matrix(unlist(gens_class_CA3_349),
                                       nrow = length(gens_class_CA3_349),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_349) = c("Gene","Express","Log2FC","P-value")


CA3_349_up = length(gens_class_CA3_349$Express[gens_class_CA3_349$Express == "UP"])
CA3_349_down = length(gens_class_CA3_349$Express[gens_class_CA3_349$Express == "DOWN"])

# CA3_350.
class_gens_CA3_350 = function(y,gens_CA3_350 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  
  ex = as.data.frame(sapply(CA3_350_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  
  g_mean = mean(df)
  
  mitpv = mean(pv)
  
  if (g_mean > 0) {
    gens_CA3_350[nrow(gens_CA3_350) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_350[nrow(gens_CA3_350) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_350)
}

gens_class_CA3_350 = lapply(gens_threshold_CA3_350,class_gens_CA3_350)

gens_class_CA3_350 = data.frame(matrix(unlist(gens_class_CA3_350),
                                       nrow = length(gens_class_CA3_350),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_350) = c("Gene","Express","Log2FC","P-value")


CA3_350_up = length(gens_class_CA3_350$Express[gens_class_CA3_350$Express == "UP"])
CA3_350_down = length(gens_class_CA3_350$Express[gens_class_CA3_350$Express == "DOWN"])

# CA3_351.
class_gens_CA3_351 = function(y,gens_CA3_351 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  
  ex = as.data.frame(sapply(CA3_351_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  
  g_mean = mean(df)
  
  mitpv = mean(pv)
  
  if (g_mean > 0) {
    gens_CA3_351[nrow(gens_CA3_351) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_351[nrow(gens_CA3_351) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_351)
}

gens_class_CA3_351 = lapply(gens_threshold_CA3_351,class_gens_CA3_351)

gens_class_CA3_351 = data.frame(matrix(unlist(gens_class_CA3_351),
                                       nrow = length(gens_class_CA3_351),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_351) = c("Gene","Express","Log2FC","P-value")


CA3_351_up = length(gens_class_CA3_351$Express[gens_class_CA3_351$Express == "UP"])
CA3_351_down = length(gens_class_CA3_351$Express[gens_class_CA3_351$Express == "DOWN"])

# CA3_352.
class_gens_CA3_352 = function(y,gens_CA3_352 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){

  ex = as.data.frame(sapply(CA3_352_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)

  g_mean = mean(df)

  mitpv = mean(pv)

  if (g_mean > 0) {
    gens_CA3_352[nrow(gens_CA3_352) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_352[nrow(gens_CA3_352) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_352)
}

gens_class_CA3_352 = lapply(gens_threshold_CA3_352,class_gens_CA3_352)

gens_class_CA3_352 = data.frame(matrix(unlist(gens_class_CA3_352),
                                       nrow = length(gens_class_CA3_352),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_352) = c("Gene","Express","Log2FC","P-value")


CA3_352_up = length(gens_class_CA3_352$Express[gens_class_CA3_352$Express == "UP"])
CA3_352_down = length(gens_class_CA3_352$Express[gens_class_CA3_352$Express == "DOWN"])

# CA3_353.
class_gens_CA3_353 = function(y,gens_CA3_353 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){

  ex = as.data.frame(sapply(CA3_353_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)

  g_mean = mean(df)

  mitpv = mean(pv)

  if (g_mean > 0) {
    gens_CA3_353[nrow(gens_CA3_353) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_353[nrow(gens_CA3_353) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_353)
}

gens_class_CA3_353 = lapply(gens_threshold_CA3_353,class_gens_CA3_353)

gens_class_CA3_353 = data.frame(matrix(unlist(gens_class_CA3_353),
                                       nrow = length(gens_class_CA3_353),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_353) = c("Gene","Express","Log2FC","P-value")


CA3_353_up = length(gens_class_CA3_353$Express[gens_class_CA3_353$Express == "UP"])
CA3_353_down = length(gens_class_CA3_353$Express[gens_class_CA3_353$Express == "DOWN"])

# CA3_354.
class_gens_CA3_354 = function(y,gens_CA3_354 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){

  ex = as.data.frame(sapply(CA3_354_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)

  g_mean = mean(df)

  mitpv = mean(pv)

  if (g_mean > 0) {
    gens_CA3_354[nrow(gens_CA3_354) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_354[nrow(gens_CA3_354) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_354)
}

gens_class_CA3_354 = lapply(gens_threshold_CA3_354,class_gens_CA3_354)

gens_class_CA3_354 = data.frame(matrix(unlist(gens_class_CA3_354),
                                       nrow = length(gens_class_CA3_354),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_354) = c("Gene","Express","Log2FC","P-value")


CA3_354_up = length(gens_class_CA3_354$Express[gens_class_CA3_354$Express == "UP"])
CA3_354_down = length(gens_class_CA3_354$Express[gens_class_CA3_354$Express == "DOWN"])

# CA3_355.
class_gens_CA3_355 = function(y,gens_CA3_355 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  
  ex = as.data.frame(sapply(CA3_355_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  
  g_mean = mean(df)
  
  mitpv = mean(pv)
  
  if (g_mean > 0) {
    gens_CA3_355[nrow(gens_CA3_355) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_355[nrow(gens_CA3_355) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_355)
}

gens_class_CA3_355 = lapply(gens_threshold_CA3_355,class_gens_CA3_355)

gens_class_CA3_355 = data.frame(matrix(unlist(gens_class_CA3_355),
                                       nrow = length(gens_class_CA3_355),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_355) = c("Gene","Express","Log2FC","P-value")


CA3_355_up = length(gens_class_CA3_355$Express[gens_class_CA3_355$Express == "UP"])
CA3_355_down = length(gens_class_CA3_355$Express[gens_class_CA3_355$Express == "DOWN"])

# CA3_356.
class_gens_CA3_356 = function(y,gens_CA3_356 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){

  ex = as.data.frame(sapply(CA3_356_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  
  g_mean = mean(df)
  
  mitpv = mean(pv)
  
  if (g_mean > 0) {
    gens_CA3_356[nrow(gens_CA3_356) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_356[nrow(gens_CA3_356) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_356)
}

gens_class_CA3_356 = lapply(gens_threshold_CA3_356,class_gens_CA3_356)

gens_class_CA3_356 = data.frame(matrix(unlist(gens_class_CA3_356),
                                       nrow = length(gens_class_CA3_356),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_356) = c("Gene","Express","Log2FC","P-value")


CA3_356_up = length(gens_class_CA3_356$Express[gens_class_CA3_356$Express == "UP"])
CA3_356_down = length(gens_class_CA3_356$Express[gens_class_CA3_356$Express == "DOWN"])

# CA3_357.
class_gens_CA3_357 = function(y,gens_CA3_357 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  
  ex = as.data.frame(sapply(CA3_357_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  
  g_mean = mean(df)
  
  mitpv = mean(pv)
  
  if (g_mean > 0) {
    gens_CA3_357[nrow(gens_CA3_357) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_357[nrow(gens_CA3_357) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_357)
}

gens_class_CA3_357 = lapply(gens_threshold_CA3_357,class_gens_CA3_357)

gens_class_CA3_357 = data.frame(matrix(unlist(gens_class_CA3_357),
                                       nrow = length(gens_class_CA3_357),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_357) = c("Gen","Express","Log2FC","P-value")


CA3_357_up = length(gens_class_CA3_357$Express[gens_class_CA3_357$Express == "UP"])
CA3_357_down = length(gens_class_CA3_357$Express[gens_class_CA3_357$Express == "DOWN"])

# CA3_358.
class_gens_CA3_358 = function(y,gens_CA3_358 = data.frame(Gene = character(),
                                                          Express = character(),
                                                          FC = numeric(),
                                                          PV = numeric(),
                                                          stringsAsFactors = F)){
  
  ex = as.data.frame(sapply(CA3_358_sig,function (x) x[row.names(x) == y,]))
  df = ex[2,]
  df = unlist(df)
  pv = ex[6,]
  pv = unlist(pv)
  
  g_mean = mean(df)
  
  mitpv = mean(pv)
  
  if (g_mean > 0) {
    gens_CA3_358[nrow(gens_CA3_358) + 1,] = c(y,"UP",g_mean,mitpv)
  } else {
    gens_CA3_358[nrow(gens_CA3_358) + 1,] = c(y,"DOWN",g_mean,mitpv)
  }
  return(gens_CA3_358)
}

gens_class_CA3_358 = lapply(gens_threshold_CA3_358,class_gens_CA3_358)

gens_class_CA3_358 = data.frame(matrix(unlist(gens_class_CA3_358),
                                       nrow = length(gens_class_CA3_358),
                                       byrow = TRUE),stringsAsFactors = F)
colnames(gens_class_CA3_358) = c("Gene","Express","Log2FC","P-value")


CA3_358_up = length(gens_class_CA3_358$Express[gens_class_CA3_358$Express == "UP"])
CA3_358_down = length(gens_class_CA3_358$Express[gens_class_CA3_358$Express == "DOWN"])


rc_totals = c(length(sig_gens_CA3_349),length(sig_gens_CA3_350),length(sig_gens_CA3_351),
              length(sig_gens_CA3_352),length(sig_gens_CA3_353),length(sig_gens_CA3_354),
              length(sig_gens_CA3_355),length(sig_gens_CA3_356),length(sig_gens_CA3_357),
              length(sig_gens_CA3_358))
rc_threshold = c(length(gens_threshold_CA3_349),length(gens_threshold_CA3_350),
               length(gens_threshold_CA3_351),length(gens_threshold_CA3_352),
               length(gens_threshold_CA3_353),length(gens_threshold_CA3_354),
               length(gens_threshold_CA3_355),length(gens_threshold_CA3_356),
               length(gens_threshold_CA3_357),length(gens_threshold_CA3_358))
rc_up = c(CA3_349_up,CA3_350_up,CA3_351_up,CA3_352_up,CA3_353_up,CA3_354_up,CA3_355_up,CA3_356_up,
          CA3_357_up,CA3_358_up)
rc_down = c(CA3_349_down,CA3_350_down,CA3_351_down,CA3_352_down,CA3_353_down,CA3_354_down,
            CA3_355_down,CA3_356_down,CA3_357_down,CA3_358_down)
summary_table = data.frame(colnames(genecount_df),rc_totals,rc_threshold,rc_up,rc_down)
colnames(summary_table) = c("Types","Significant Genes","Threshold Genes","Threshold Genes UP",
                          "Threshold Genes DOWN")


llin_gens_up = unique(c(gens_class_CA3_349$Gen[gens_class_CA3_349$Express == "UP"],
                        gens_class_CA3_350$Gen[gens_class_CA3_350$Express == "UP"],
                        gens_class_CA3_351$Gen[gens_class_CA3_351$Express == "UP"],
                        gens_class_CA3_352$Gen[gens_class_CA3_352$Express == "UP"],
                        gens_class_CA3_353$Gen[gens_class_CA3_353$Express == "UP"],
                        gens_class_CA3_354$Gen[gens_class_CA3_354$Express == "UP"],
                        gens_class_CA3_355$Gen[gens_class_CA3_355$Express == "UP"],
                        gens_class_CA3_356$Gen[gens_class_CA3_356$Express == "UP"],
                        gens_class_CA3_357$Gen[gens_class_CA3_357$Express == "UP"],
                        gens_class_CA3_358$Gen[gens_class_CA3_358$Express == "UP"]))
llin_gens_down = unique(c(gens_class_CA3_349$Gen[gens_class_CA3_349$Express == "DOWN"],
                          gens_class_CA3_350$Gen[gens_class_CA3_350$Express == "DOWN"],
                          gens_class_CA3_351$Gen[gens_class_CA3_351$Express == "DOWN"],
                          gens_class_CA3_352$Gen[gens_class_CA3_352$Express == "DOWN"],
                          gens_class_CA3_353$Gen[gens_class_CA3_353$Express == "DOWN"],
                          gens_class_CA3_354$Gen[gens_class_CA3_354$Express == "DOWN"],
                          gens_class_CA3_355$Gen[gens_class_CA3_355$Express == "DOWN"],
                          gens_class_CA3_356$Gen[gens_class_CA3_356$Express == "DOWN"],
                          gens_class_CA3_357$Gen[gens_class_CA3_357$Express == "DOWN"],
                          gens_class_CA3_358$Gen[gens_class_CA3_358$Express == "DOWN"]))

png(filename = "Pattern_349_Mossy.png",width = 800,height = 600)
pattern_349_CA3 = hist(genecount_df$`349_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[35],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_349_CA3 = transform(table(cut(genecount_df$`349_CA3`,breaks = 10)))

png(filename = "Pattern_350_Mossy.png",width = 800,height = 600)
pattern_350_CA3 = hist(genecount_df$`350_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[36],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_350_CA3 = transform(table(cut(genecount_df$`350_CA3`,breaks = 10)))

png(filename = "Pattern_351_CA3-ve.png",width = 800,height = 600)
pattern_351_CA3 = hist(genecount_df$`351_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[37],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_351_CA3 = transform(table(cut(genecount_df$`351_CA3`,breaks = 10)))

png(filename = "Pattern_352_CA3-ve.png",width = 800,height = 600)
pattern_352_CA3 = hist(genecount_df$`352_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[38],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_352_CA3 = transform(table(cut(genecount_df$`352_CA3`,breaks = 10)))

png(filename = "Pattern_353_CA3-ve.png",width = 800,height = 600)
pattern_353_CA3 = hist(genecount_df$`353_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[39],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_353_CA3 = transform(table(cut(genecount_df$`353_CA3`,breaks = 10)))

png(filename = "Pattern_354_CA3-ve.png",width = 800,height = 600)
pattern_354_CA3 = hist(genecount_df$`354_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[40],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_354_CA3 = transform(table(cut(genecount_df$`354_CA3`,breaks = 10)))

png(filename = "Pattern_355_CA3-ve.png",width = 800,height = 600)
pattern_355_CA3 = hist(genecount_df$`355_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[41],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_355_CA3 = transform(table(cut(genecount_df$`355_CA3`,breaks = 10)))

png(filename = "Pattern_356_CA3-do.png",width = 800,height = 600)
pattern_356_CA3 = hist(genecount_df$`356_CA3`,breaks = 10,main = "Expression Pattern",
                     xlab = colnames(genecount_df)[42],ylab = "Frequency",col = "lightblue")
dev.off()
pattern_t_356_CA3 = transform(table(cut(genecount_df$`356_CA3`,breaks = 10)))

png(filename = "Pattern_357_CA3-do.png",width = 800,height = 600)
pattern_357_CA3 = hist(genecount_df$`357_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[43],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_357_CA3 = transform(table(cut(genecount_df$`357_CA3`,breaks = 10)))

png(filename = "Pattern_358_CA3-do.png",width = 800,height = 600)
pattern_358_CA3 = hist(genecount_df$`358_CA3`,breaks = 10,main = "Expression Pattern",
                     ylim = c(0,1000),xlab = colnames(genecount_df)[44],ylab = "Frequency",
                     col = "lightblue")
dev.off()
pattern_t_358_CA3 = transform(table(cut(genecount_df$`358_CA3`,breaks = 10)))


UP_CA3_349 = gens_class_CA3_349$Gen[gens_class_CA3_349$Express == "UP"]
DOWN_CA3_349 = gens_class_CA3_349$Gen[gens_class_CA3_349$Express == "DOWN"]
UP_CA3_350 = gens_class_CA3_350$Gen[gens_class_CA3_350$Express == "UP"]
DOWN_CA3_350 = gens_class_CA3_350$Gen[gens_class_CA3_350$Express == "DOWN"]
UP_CA3_351 = gens_class_CA3_351$Gen[gens_class_CA3_351$Express == "UP"]
DOWN_CA3_351 = gens_class_CA3_351$Gen[gens_class_CA3_351$Express == "DOWN"]
UP_CA3_352 = gens_class_CA3_352$Gen[gens_class_CA3_352$Express == "UP"]
DOWN_CA3_352 = gens_class_CA3_352$Gen[gens_class_CA3_352$Express == "DOWN"]
UP_CA3_353 = gens_class_CA3_353$Gen[gens_class_CA3_353$Express == "UP"]
DOWN_CA3_353 = gens_class_CA3_353$Gen[gens_class_CA3_353$Express == "DOWN"]
UP_CA3_354 = gens_class_CA3_354$Gen[gens_class_CA3_354$Express == "UP"]
DOWN_CA3_354 = gens_class_CA3_354$Gen[gens_class_CA3_354$Express == "DOWN"]
UP_CA3_355 = gens_class_CA3_355$Gen[gens_class_CA3_355$Express == "UP"]
DOWN_CA3_355 = gens_class_CA3_355$Gen[gens_class_CA3_355$Express == "DOWN"]
UP_CA3_356 = gens_class_CA3_356$Gen[gens_class_CA3_356$Express == "UP"]
DOWN_CA3_356 = gens_class_CA3_356$Gen[gens_class_CA3_356$Express == "DOWN"]
UP_CA3_357 = gens_class_CA3_357$Gen[gens_class_CA3_357$Express == "UP"]
DOWN_CA3_357 = gens_class_CA3_357$Gen[gens_class_CA3_357$Express == "DOWN"]
UP_CA3_358 = gens_class_CA3_358$Gen[gens_class_CA3_358$Express == "UP"]
DOWN_CA3_358 = gens_class_CA3_358$Gen[gens_class_CA3_358$Express == "DOWN"]

gcount_sig = data.frame(genecount_df >= threshold)
gcount_sig = gcount_sig[rowSums(gcount_sig) != 0,]
colnames(gcount_sig) = c("349_CA3","350_CA3","351_CA3","352_CA3","353_CA3","354_CA3",
                         "355_CA3","356_CA3","357_CA3","358_CA3")

gcount_sig$`349_CA3`[rownames(gcount_sig) %in% UP_CA3_349] = "UP"
gcount_sig$`349_CA3`[rownames(gcount_sig) %in% DOWN_CA3_349] = "DOWN"
gcount_sig$`350_CA3`[rownames(gcount_sig) %in% UP_CA3_350] = "UP"
gcount_sig$`350_CA3`[rownames(gcount_sig) %in% DOWN_CA3_350] = "DOWN"
gcount_sig$`351_CA3`[rownames(gcount_sig) %in% UP_CA3_351] = "UP"
gcount_sig$`351_CA3`[rownames(gcount_sig) %in% DOWN_CA3_351] = "DOWN"
gcount_sig$`352_CA3`[rownames(gcount_sig) %in% UP_CA3_352] = "UP"
gcount_sig$`352_CA3`[rownames(gcount_sig) %in% DOWN_CA3_352] = "DOWN"
gcount_sig$`353_CA3`[rownames(gcount_sig) %in% UP_CA3_353] = "UP"
gcount_sig$`353_CA3`[rownames(gcount_sig) %in% DOWN_CA3_353] = "DOWN"
gcount_sig$`354_CA3`[rownames(gcount_sig) %in% UP_CA3_354] = "UP"
gcount_sig$`354_CA3`[rownames(gcount_sig) %in% DOWN_CA3_354] = "DOWN"
gcount_sig$`355_CA3`[rownames(gcount_sig) %in% UP_CA3_355] = "UP"
gcount_sig$`355_CA3`[rownames(gcount_sig) %in% DOWN_CA3_355] = "DOWN"
gcount_sig$`356_CA3`[rownames(gcount_sig) %in% UP_CA3_356] = "UP"
gcount_sig$`356_CA3`[rownames(gcount_sig) %in% DOWN_CA3_356] = "DOWN"
gcount_sig$`357_CA3`[rownames(gcount_sig) %in% UP_CA3_357] = "UP"
gcount_sig$`357_CA3`[rownames(gcount_sig) %in% DOWN_CA3_357] = "DOWN"
gcount_sig$`358_CA3`[rownames(gcount_sig) %in% UP_CA3_358] = "UP"
gcount_sig$`358_CA3`[rownames(gcount_sig) %in% DOWN_CA3_358] = "DOWN"

UP = unname(apply(gcount_sig,1,function (x) paste(names(x[x == "UP"]), collapse = ",")))
count_UP = apply(gcount_sig,1,function (x) length(x[x == "UP"]))
coinc_UP = data.frame(count_UP,UP)
colnames(coinc_UP) = c("Types number","Types")
coinc_UP = coinc_UP[coinc_UP$`Nombre Types` !=0,]

DOWN = unname(apply(gcount_sig,1,function (x) paste(names(x[x == "DOWN"]), collapse = ",")))
count_DOWN = apply(gcount_sig,1,function (x) length(x[x == "DOWN"]))
coinc_DOWN = data.frame(count_DOWN,DOWN)
colnames(coinc_DOWN) = c("Types number","Types")
coinc_DOWN = coinc_DOWN[coinc_DOWN$`Nombre Types` !=0,]

png("Coinc_UP_2DEG.png",width = 800,height = 600)
par(mar = c(10,4,5,1))
barplot(table(coinc_UP[coinc_UP$`Types number` == 2,length(coinc_UP)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()

png("Coinc_UP_3DEG.png",width = 800,height = 600)
par(mar = c(12.5,4,5,1))
barplot(table(coinc_UP[coinc_UP$`Types number` == 3,length(coinc_UP)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()

png("Coinc_DOWN_2DEG.png",width = 800,height = 600)
par(mar = c(10,4,5,1))
barplot(table(coinc_DOWN[coinc_DOWN$`Types number` == 2,length(coinc_DOWN)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()

png("Coinc_DOWN_3DEG.png",width = 800,height = 600)
par(mar = c(12.5,4,5,1))
barplot(table(coinc_DOWN[coinc_DOWN$`Types number` == 3,length(coinc_DOWN)]),col = "lightblue",
        ylab = "Number of genes",las = 2)
abline(h = 5,col = "red",lwd = 2,lty = 2)
dev.off()

png("Coinc_UP_349_Mossy.png")
barplot(table(coinc_UP$`Types number`[grep("349_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 349_Mossy)")
dev.off()

coinc_UP_CA3_349 = coinc_UP[grep("349_CA3",coinc_UP$Types),]
coinc_UP_CA3_349_unique = rownames(coinc_UP_CA3_349[coinc_UP_CA3_349$`Types number` == 1,])
coinc_UP_CA3_349_d = rownames(coinc_UP_CA3_349[coinc_UP_CA3_349$`Types number` == 2,])
coinc_UP_CA3_349_t = rownames(coinc_UP_CA3_349[coinc_UP_CA3_349$`Types number` == 3,])
coinc_UP_CA3_349_q = rownames(coinc_UP_CA3_349[coinc_UP_CA3_349$`Types number` == 4,])
coinc_UP_CA3_349_c = rownames(coinc_UP_CA3_349[coinc_UP_CA3_349$`Types number` == 5,])
coinc_UP_CA3_349_s = rownames(coinc_UP_CA3_349[coinc_UP_CA3_349$`Types number` == 6,])


png("Coinc_DOWN_349_Mossy.png")
barplot(table(coinc_DOWN$`Types number`[grep("349_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 349_Mossy)")
dev.off()

coinc_DOWN_CA3_349 = coinc_DOWN[grep("349_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_349_unique = rownames(coinc_DOWN_CA3_349[coinc_DOWN_CA3_349$`Types number` == 1,])
coinc_DOWN_CA3_349_d = rownames(coinc_DOWN_CA3_349[coinc_DOWN_CA3_349$`Types number` == 2,])
coinc_DOWN_CA3_349_t = rownames(coinc_DOWN_CA3_349[coinc_DOWN_CA3_349$`Types number` == 3,])
coinc_DOWN_CA3_349_q = rownames(coinc_DOWN_CA3_349[coinc_DOWN_CA3_349$`Types number` == 4,])
coinc_DOWN_CA3_349_c = rownames(coinc_DOWN_CA3_349[coinc_DOWN_CA3_349$`Types number` == 5,])
coinc_DOWN_CA3_349_s = rownames(coinc_DOWN_CA3_349[coinc_DOWN_CA3_349$`Types number` == 6,])


png("Coinc_UP_350_Mossy.png")
barplot(table(coinc_UP$`Types number`[grep("350_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 350_Mossy)")
dev.off()

coinc_UP_CA3_350 = coinc_UP[grep("350_CA3",coinc_UP$Types),]
coinc_UP_CA3_350_unique = rownames(coinc_UP_CA3_350[coinc_UP_CA3_350$`Types number` == 1,])
coinc_UP_CA3_350_d = rownames(coinc_UP_CA3_350[coinc_UP_CA3_350$`Types number` == 2,])
coinc_UP_CA3_350_t = rownames(coinc_UP_CA3_350[coinc_UP_CA3_350$`Types number` == 3,])
coinc_UP_CA3_350_q = rownames(coinc_UP_CA3_350[coinc_UP_CA3_350$`Types number` == 4,])
coinc_UP_CA3_350_c = rownames(coinc_UP_CA3_350[coinc_UP_CA3_350$`Types number` == 5,])
coinc_UP_CA3_350_s = rownames(coinc_UP_CA3_350[coinc_UP_CA3_350$`Types number` == 6,])


png("Coinc_DOWN_350_Mossy.png")
barplot(table(coinc_DOWN$`Types number`[grep("350_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 350_Mossy)")
dev.off()

coinc_DOWN_CA3_350 = coinc_DOWN[grep("350_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_350_unique = rownames(coinc_DOWN_CA3_350[coinc_DOWN_CA3_350$`Types number` == 1,])
coinc_DOWN_CA3_350_d = rownames(coinc_DOWN_CA3_350[coinc_DOWN_CA3_350$`Types number` == 2,])
coinc_DOWN_CA3_350_t = rownames(coinc_DOWN_CA3_350[coinc_DOWN_CA3_350$`Types number` == 3,])
coinc_DOWN_CA3_350_q = rownames(coinc_DOWN_CA3_350[coinc_DOWN_CA3_350$`Types number` == 4,])
coinc_DOWN_CA3_350_c = rownames(coinc_DOWN_CA3_350[coinc_DOWN_CA3_350$`Types number` == 5,])
coinc_DOWN_CA3_350_s = rownames(coinc_DOWN_CA3_350[coinc_DOWN_CA3_350$`Types number` == 6,])


png("Coinc_UP_351_CA3-ve.png")
barplot(table(coinc_UP$`Types number`[grep("351_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 351_CA3-ve)")
dev.off()

coinc_UP_CA3_351 = coinc_UP[grep("351_CA3",coinc_UP$Types),]
coinc_UP_CA3_351_unique = rownames(coinc_UP_CA3_351[coinc_UP_CA3_351$`Types number` == 1,])
coinc_UP_CA3_351_d = rownames(coinc_UP_CA3_351[coinc_UP_CA3_351$`Types number` == 2,])
coinc_UP_CA3_351_t = rownames(coinc_UP_CA3_351[coinc_UP_CA3_351$`Types number` == 3,])
coinc_UP_CA3_351_q = rownames(coinc_UP_CA3_351[coinc_UP_CA3_351$`Types number` == 4,])
coinc_UP_CA3_351_c = rownames(coinc_UP_CA3_351[coinc_UP_CA3_351$`Types number` == 5,])
coinc_UP_CA3_351_s = rownames(coinc_UP_CA3_351[coinc_UP_CA3_351$`Types number` == 6,])


png("Coinc_DOWN_351_CA3-ve.png")
barplot(table(coinc_DOWN$`Types number`[grep("351_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 351_CA3-ve)")
dev.off()

coinc_DOWN_CA3_351 = coinc_DOWN[grep("351_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_351_unique = rownames(coinc_DOWN_CA3_351[coinc_DOWN_CA3_351$`Types number` == 1,])
coinc_DOWN_CA3_351_d = rownames(coinc_DOWN_CA3_351[coinc_DOWN_CA3_351$`Types number` == 2,])
coinc_DOWN_CA3_351_t = rownames(coinc_DOWN_CA3_351[coinc_DOWN_CA3_351$`Types number` == 3,])
coinc_DOWN_CA3_351_q = rownames(coinc_DOWN_CA3_351[coinc_DOWN_CA3_351$`Types number` == 4,])
coinc_DOWN_CA3_351_c = rownames(coinc_DOWN_CA3_351[coinc_DOWN_CA3_351$`Types number` == 5,])
coinc_DOWN_CA3_351_s = rownames(coinc_DOWN_CA3_351[coinc_DOWN_CA3_351$`Types number` == 6,])


png("Coinc_UP_352_CA3-ve.png")
barplot(table(coinc_UP$`Types number`[grep("352_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 352_CA3-ve)")
dev.off()

coinc_UP_CA3_352 = coinc_UP[grep("352_CA3",coinc_UP$Types),]
coinc_UP_CA3_352_unique = rownames(coinc_UP_CA3_352[coinc_UP_CA3_352$`Types number` == 1,])
coinc_UP_CA3_352_d = rownames(coinc_UP_CA3_352[coinc_UP_CA3_352$`Types number` == 2,])
coinc_UP_CA3_352_t = rownames(coinc_UP_CA3_352[coinc_UP_CA3_352$`Types number` == 3,])
coinc_UP_CA3_352_q = rownames(coinc_UP_CA3_352[coinc_UP_CA3_352$`Types number` == 4,])
coinc_UP_CA3_352_c = rownames(coinc_UP_CA3_352[coinc_UP_CA3_352$`Types number` == 5,])
coinc_UP_CA3_352_s = rownames(coinc_UP_CA3_352[coinc_UP_CA3_352$`Types number` == 6,])


png("Coinc_DOWN_352_CA3-ve.png")
barplot(table(coinc_DOWN$`Types number`[grep("352_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 352_CA3-ve)")
dev.off()

coinc_DOWN_CA3_352 = coinc_DOWN[grep("352_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_352_unique = rownames(coinc_DOWN_CA3_352[coinc_DOWN_CA3_352$`Types number` == 1,])
coinc_DOWN_CA3_352_d = rownames(coinc_DOWN_CA3_352[coinc_DOWN_CA3_352$`Types number` == 2,])
coinc_DOWN_CA3_352_t = rownames(coinc_DOWN_CA3_352[coinc_DOWN_CA3_352$`Types number` == 3,])
coinc_DOWN_CA3_352_q = rownames(coinc_DOWN_CA3_352[coinc_DOWN_CA3_352$`Types number` == 4,])
coinc_DOWN_CA3_352_c = rownames(coinc_DOWN_CA3_352[coinc_DOWN_CA3_352$`Types number` == 5,])
coinc_DOWN_CA3_352_s = rownames(coinc_DOWN_CA3_352[coinc_DOWN_CA3_352$`Types number` == 6,])


png("Coinc_UP_353_CA3-ve.png")
barplot(table(coinc_UP$`Types number`[grep("353_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 353_CA3)")
dev.off()

coinc_UP_CA3_353 = coinc_UP[grep("353_CA3",coinc_UP$Types),]
coinc_UP_CA3_353_unique = rownames(coinc_UP_CA3_353[coinc_UP_CA3_353$`Types number` == 1,])
coinc_UP_CA3_353_d = rownames(coinc_UP_CA3_353[coinc_UP_CA3_353$`Types number` == 2,])
coinc_UP_CA3_353_t = rownames(coinc_UP_CA3_353[coinc_UP_CA3_353$`Types number` == 3,])
coinc_UP_CA3_353_q = rownames(coinc_UP_CA3_353[coinc_UP_CA3_353$`Types number` == 4,])
coinc_UP_CA3_353_c = rownames(coinc_UP_CA3_353[coinc_UP_CA3_353$`Types number` == 5,])
coinc_UP_CA3_353_s = rownames(coinc_UP_CA3_353[coinc_UP_CA3_353$`Types number` == 6,])


png("Coinc_DOWN_353_CA3-ve.png")
barplot(table(coinc_DOWN$`Types number`[grep("353_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 353_CA3-ve)")
dev.off()

coinc_DOWN_CA3_353 = coinc_DOWN[grep("353_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_353_unique = rownames(coinc_DOWN_CA3_353[coinc_DOWN_CA3_353$`Types number` == 1,])
coinc_DOWN_CA3_353_d = rownames(coinc_DOWN_CA3_353[coinc_DOWN_CA3_353$`Types number` == 2,])
coinc_DOWN_CA3_353_t = rownames(coinc_DOWN_CA3_353[coinc_DOWN_CA3_353$`Types number` == 3,])
coinc_DOWN_CA3_353_q = rownames(coinc_DOWN_CA3_353[coinc_DOWN_CA3_353$`Types number` == 4,])
coinc_DOWN_CA3_353_c = rownames(coinc_DOWN_CA3_353[coinc_DOWN_CA3_353$`Types number` == 5,])
coinc_DOWN_CA3_353_s = rownames(coinc_DOWN_CA3_353[coinc_DOWN_CA3_353$`Types number` == 6,])


png("Coinc_UP_354_CA3-ve.png")
barplot(table(coinc_UP$`Types number`[grep("354_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 354_CA3-ve)")
dev.off()

coinc_UP_CA3_354 = coinc_UP[grep("354_CA3",coinc_UP$Types),]
coinc_UP_CA3_354_unique = rownames(coinc_UP_CA3_354[coinc_UP_CA3_354$`Types number` == 1,])
coinc_UP_CA3_354_d = rownames(coinc_UP_CA3_354[coinc_UP_CA3_354$`Types number` == 2,])
coinc_UP_CA3_354_t = rownames(coinc_UP_CA3_354[coinc_UP_CA3_354$`Types number` == 3,])
coinc_UP_CA3_354_q = rownames(coinc_UP_CA3_354[coinc_UP_CA3_354$`Types number` == 4,])
coinc_UP_CA3_354_c = rownames(coinc_UP_CA3_354[coinc_UP_CA3_354$`Types number` == 5,])
coinc_UP_CA3_354_s = rownames(coinc_UP_CA3_354[coinc_UP_CA3_354$`Types number` == 6,])


png("Coinc_DOWN_354_CA3-ve.png")
barplot(table(coinc_DOWN$`Types number`[grep("354_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 354_CA3-ve)")
dev.off()

coinc_DOWN_CA3_354 = coinc_DOWN[grep("354_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_354_unique = rownames(coinc_DOWN_CA3_354[coinc_DOWN_CA3_354$`Types number` == 1,])
coinc_DOWN_CA3_354_d = rownames(coinc_DOWN_CA3_354[coinc_DOWN_CA3_354$`Types number` == 2,])
coinc_DOWN_CA3_354_t = rownames(coinc_DOWN_CA3_354[coinc_DOWN_CA3_354$`Types number` == 3,])
coinc_DOWN_CA3_354_q = rownames(coinc_DOWN_CA3_354[coinc_DOWN_CA3_354$`Types number` == 4,])
coinc_DOWN_CA3_354_c = rownames(coinc_DOWN_CA3_354[coinc_DOWN_CA3_354$`Types number` == 5,])
coinc_DOWN_CA3_354_s = rownames(coinc_DOWN_CA3_354[coinc_DOWN_CA3_354$`Types number` == 6,])


png("Coinc_UP_355_CA3-ve.png")
barplot(table(coinc_UP$`Types number`[grep("355_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 355_CA3-ve)")
dev.off()

coinc_UP_CA3_355 = coinc_UP[grep("355_CA3",coinc_UP$Types),]
coinc_UP_CA3_355_unique = rownames(coinc_UP_CA3_355[coinc_UP_CA3_355$`Types number` == 1,])
coinc_UP_CA3_355_d = rownames(coinc_UP_CA3_355[coinc_UP_CA3_355$`Types number` == 2,])
coinc_UP_CA3_355_t = rownames(coinc_UP_CA3_355[coinc_UP_CA3_355$`Types number` == 3,])
coinc_UP_CA3_355_q = rownames(coinc_UP_CA3_355[coinc_UP_CA3_355$`Types number` == 4,])
coinc_UP_CA3_355_c = rownames(coinc_UP_CA3_355[coinc_UP_CA3_355$`Types number` == 5,])
coinc_UP_CA3_355_s = rownames(coinc_UP_CA3_355[coinc_UP_CA3_355$`Types number` == 6,])


png("Coinc_DOWN_355_CA3-ve.png")
barplot(table(coinc_DOWN$`Types number`[grep("355_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 355_CA3-ve)")
dev.off()

coinc_DOWN_CA3_355 = coinc_DOWN[grep("355_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_355_unique = rownames(coinc_DOWN_CA3_355[coinc_DOWN_CA3_355$`Types number` == 1,])
coinc_DOWN_CA3_355_d = rownames(coinc_DOWN_CA3_355[coinc_DOWN_CA3_355$`Types number` == 2,])
coinc_DOWN_CA3_355_t = rownames(coinc_DOWN_CA3_355[coinc_DOWN_CA3_355$`Types number` == 3,])
coinc_DOWN_CA3_355_q = rownames(coinc_DOWN_CA3_355[coinc_DOWN_CA3_355$`Types number` == 4,])
coinc_DOWN_CA3_355_c = rownames(coinc_DOWN_CA3_355[coinc_DOWN_CA3_355$`Types number` == 5,])
coinc_DOWN_CA3_355_s = rownames(coinc_DOWN_CA3_355[coinc_DOWN_CA3_355$`Types number` == 6,])


png("Coinc_UP_356_CA3-do.png")
barplot(table(coinc_UP$`Types number`[grep("356_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 356_CA3-do)")
dev.off()

coinc_UP_CA3_356 = coinc_UP[grep("356_CA3",coinc_UP$Types),]
coinc_UP_CA3_356_unique = rownames(coinc_UP_CA3_356[coinc_UP_CA3_356$`Types number` == 1,])
coinc_UP_CA3_356_d = rownames(coinc_UP_CA3_356[coinc_UP_CA3_356$`Types number` == 2,])
coinc_UP_CA3_356_t = rownames(coinc_UP_CA3_356[coinc_UP_CA3_356$`Types number` == 3,])
coinc_UP_CA3_356_q = rownames(coinc_UP_CA3_356[coinc_UP_CA3_356$`Types number` == 4,])
coinc_UP_CA3_356_c = rownames(coinc_UP_CA3_356[coinc_UP_CA3_356$`Types number` == 5,])
coinc_UP_CA3_356_s = rownames(coinc_UP_CA3_356[coinc_UP_CA3_356$`Types number` == 6,])


png("Coinc_DOWN_356_CA3-do.png")
barplot(table(coinc_DOWN$`Types number`[grep("356_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 356_CA3-do)")
dev.off()

coinc_DOWN_CA3_356 = coinc_DOWN[grep("356_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_356_unique = rownames(coinc_DOWN_CA3_356[coinc_DOWN_CA3_356$`Types number` == 1,])
coinc_DOWN_CA3_356_d = rownames(coinc_DOWN_CA3_356[coinc_DOWN_CA3_356$`Types number` == 2,])
coinc_DOWN_CA3_356_t = rownames(coinc_DOWN_CA3_356[coinc_DOWN_CA3_356$`Types number` == 3,])
coinc_DOWN_CA3_356_q = rownames(coinc_DOWN_CA3_356[coinc_DOWN_CA3_356$`Types number` == 4,])
coinc_DOWN_CA3_356_c = rownames(coinc_DOWN_CA3_356[coinc_DOWN_CA3_356$`Types number` == 5,])
coinc_DOWN_CA3_356_s = rownames(coinc_DOWN_CA3_356[coinc_DOWN_CA3_356$`Types number` == 6,])


png("Coinc_UP_357_CA3-do.png")
barplot(table(coinc_UP$`Types number`[grep("357_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 357_CA3-do)")
dev.off()

coinc_UP_CA3_357 = coinc_UP[grep("357_CA3",coinc_UP$Types),]
coinc_UP_CA3_357_unique = rownames(coinc_UP_CA3_357[coinc_UP_CA3_357$`Types number` == 1,])
coinc_UP_CA3_357_d = rownames(coinc_UP_CA3_357[coinc_UP_CA3_357$`Types number` == 2,])
coinc_UP_CA3_357_t = rownames(coinc_UP_CA3_357[coinc_UP_CA3_357$`Types number` == 3,])
coinc_UP_CA3_357_q = rownames(coinc_UP_CA3_357[coinc_UP_CA3_357$`Types number` == 4,])
coinc_UP_CA3_357_c = rownames(coinc_UP_CA3_357[coinc_UP_CA3_357$`Types number` == 5,])
coinc_UP_CA3_357_s = rownames(coinc_UP_CA3_357[coinc_UP_CA3_357$`Types number` == 6,])


png("Coinc_DOWN_357_CA3-do.png")
barplot(table(coinc_DOWN$`Types number`[grep("357_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 357_CA3-do)")
dev.off()

coinc_DOWN_CA3_357 = coinc_DOWN[grep("357_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_357_unique = rownames(coinc_DOWN_CA3_357[coinc_DOWN_CA3_357$`Types number` == 1,])
coinc_DOWN_CA3_357_d = rownames(coinc_DOWN_CA3_357[coinc_DOWN_CA3_357$`Types number` == 2,])
coinc_DOWN_CA3_357_t = rownames(coinc_DOWN_CA3_357[coinc_DOWN_CA3_357$`Types number` == 3,])
coinc_DOWN_CA3_357_q = rownames(coinc_DOWN_CA3_357[coinc_DOWN_CA3_357$`Types number` == 4,])
coinc_DOWN_CA3_357_c = rownames(coinc_DOWN_CA3_357[coinc_DOWN_CA3_357$`Types number` == 5,])
coinc_DOWN_CA3_357_s = rownames(coinc_DOWN_CA3_357[coinc_DOWN_CA3_357$`Types number` == 6,])


png("Coinc_UP_358_CA3-do.png")
barplot(table(coinc_UP$`Types number`[grep("358_CA3",coinc_UP$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG UP is shared (1 = 358_CA3-do)")
dev.off()

coinc_UP_CA3_358 = coinc_UP[grep("358_CA3",coinc_UP$Types),]
coinc_UP_CA3_358_unique = rownames(coinc_UP_CA3_358[coinc_UP_CA3_358$`Types number` == 1,])
coinc_UP_CA3_358_d = rownames(coinc_UP_CA3_358[coinc_UP_CA3_358$`Types number` == 2,])
coinc_UP_CA3_358_t = rownames(coinc_UP_CA3_358[coinc_UP_CA3_358$`Types number` == 3,])
coinc_UP_CA3_358_q = rownames(coinc_UP_CA3_358[coinc_UP_CA3_358$`Types number` == 4,])
coinc_UP_CA3_358_c = rownames(coinc_UP_CA3_358[coinc_UP_CA3_358$`Types number` == 5,])
coinc_UP_CA3_358_s = rownames(coinc_UP_CA3_358[coinc_UP_CA3_358$`Types number` == 6,])


png("Coinc_DOWN_358_CA3-do.png")
barplot(table(coinc_DOWN$`Types number`[grep("358_CA3",coinc_DOWN$Types)]),
        ylab = "Number of genes",col = "lightgreen",
        xlab = "Num. of neuronal types with which the DEG DOWN is shared (1 = 358_CA3-do)")
dev.off()

coinc_DOWN_CA3_358 = coinc_DOWN[grep("358_CA3",coinc_DOWN$Types),]
coinc_DOWN_CA3_358_unique = rownames(coinc_DOWN_CA3_358[coinc_DOWN_CA3_358$`Types number` == 1,])
coinc_DOWN_CA3_358_d = rownames(coinc_DOWN_CA3_358[coinc_DOWN_CA3_358$`Types number` == 2,])
coinc_DOWN_CA3_358_t = rownames(coinc_DOWN_CA3_358[coinc_DOWN_CA3_358$`Types number` == 3,])
coinc_DOWN_CA3_358_q = rownames(coinc_DOWN_CA3_358[coinc_DOWN_CA3_358$`Types number` == 4,])
coinc_DOWN_CA3_358_c = rownames(coinc_DOWN_CA3_358[coinc_DOWN_CA3_358$`Types number` == 5,])
coinc_DOWN_CA3_358_s = rownames(coinc_DOWN_CA3_358[coinc_DOWN_CA3_358$`Types number` == 6,])


unique_up = c(length(coinc_UP_CA3_349_unique),length(coinc_UP_CA3_350_unique),
             length(coinc_UP_CA3_351_unique),length(coinc_UP_CA3_352_unique),
             length(coinc_UP_CA3_353_unique),length(coinc_UP_CA3_354_unique),
             length(coinc_UP_CA3_355_unique),length(coinc_UP_CA3_356_unique),
             length(coinc_UP_CA3_357_unique),length(coinc_UP_CA3_358_unique))
two_up = c(length(coinc_UP_CA3_349_d),length(coinc_UP_CA3_350_d),length(coinc_UP_CA3_351_d),
           length(coinc_UP_CA3_352_d),length(coinc_UP_CA3_353_d),length(coinc_UP_CA3_354_d),
           length(coinc_UP_CA3_355_d),length(coinc_UP_CA3_356_d),length(coinc_UP_CA3_357_d),
           length(coinc_UP_CA3_358_d))
three_up = c(length(coinc_UP_CA3_349_t),length(coinc_UP_CA3_350_t),length(coinc_UP_CA3_351_t),
            length(coinc_UP_CA3_352_t),length(coinc_UP_CA3_353_t),length(coinc_UP_CA3_354_t),
            length(coinc_UP_CA3_355_t),length(coinc_UP_CA3_356_t),length(coinc_UP_CA3_357_t),
            length(coinc_UP_CA3_358_t))
four_up = c(length(coinc_UP_CA3_349_q),length(coinc_UP_CA3_350_q),length(coinc_UP_CA3_351_q),
              length(coinc_UP_CA3_352_q),length(coinc_UP_CA3_353_q),length(coinc_UP_CA3_354_q),
              length(coinc_UP_CA3_355_q),length(coinc_UP_CA3_356_q),length(coinc_UP_CA3_357_q),
              length(coinc_UP_CA3_358_q))
five_up = c(length(coinc_UP_CA3_349_c),length(coinc_UP_CA3_350_c),length(coinc_UP_CA3_351_c),
            length(coinc_UP_CA3_352_c),length(coinc_UP_CA3_353_c),length(coinc_UP_CA3_354_c),
            length(coinc_UP_CA3_355_c),length(coinc_UP_CA3_356_c),length(coinc_UP_CA3_357_c),
            length(coinc_UP_CA3_358_c))
six_up = c(length(coinc_UP_CA3_349_s),length(coinc_UP_CA3_350_s),length(coinc_UP_CA3_351_s),
           length(coinc_UP_CA3_352_s),length(coinc_UP_CA3_353_s),length(coinc_UP_CA3_354_s),
           length(coinc_UP_CA3_355_s),length(coinc_UP_CA3_356_s),length(coinc_UP_CA3_357_s),
           length(coinc_UP_CA3_358_s))
tsummary_up = data.frame(unique_up,two_up,three_up,four_up,five_up,six_up)
colnames(tsummary_up) = c("Unique","Two types","Three types","Four types","Five types","Six types")
rownames(tsummary_up) = c("349_Mossy","350_Mossy","351_CA3-ve","352_CA3-ve","353_CA3-ve","354_CA3-ve",
                        "355_CA3-ve","356_CA3-do","357_CA3-do","358_CA3-do")


unique_down = c(length(coinc_DOWN_CA3_349_unique),length(coinc_DOWN_CA3_350_unique),
             length(coinc_DOWN_CA3_351_unique),length(coinc_DOWN_CA3_352_unique),
             length(coinc_DOWN_CA3_353_unique),length(coinc_DOWN_CA3_354_unique),
             length(coinc_DOWN_CA3_355_unique),length(coinc_DOWN_CA3_356_unique),
             length(coinc_DOWN_CA3_357_unique),length(coinc_DOWN_CA3_358_unique))
two_down = c(length(coinc_DOWN_CA3_349_d),length(coinc_DOWN_CA3_350_d),length(coinc_DOWN_CA3_351_d),
           length(coinc_DOWN_CA3_352_d),length(coinc_DOWN_CA3_353_d),length(coinc_DOWN_CA3_354_d),
           length(coinc_DOWN_CA3_355_d),length(coinc_DOWN_CA3_356_d),length(coinc_DOWN_CA3_357_d),
           length(coinc_DOWN_CA3_358_d))
three_down = c(length(coinc_DOWN_CA3_349_t),length(coinc_DOWN_CA3_350_t),length(coinc_DOWN_CA3_351_t),
            length(coinc_DOWN_CA3_352_t),length(coinc_DOWN_CA3_353_t),length(coinc_DOWN_CA3_354_t),
            length(coinc_DOWN_CA3_355_t),length(coinc_DOWN_CA3_356_t),length(coinc_DOWN_CA3_357_t),
            length(coinc_DOWN_CA3_358_t))
four_down = c(length(coinc_DOWN_CA3_349_q),length(coinc_DOWN_CA3_350_q),
                length(coinc_DOWN_CA3_351_q),length(coinc_DOWN_CA3_352_q),
                length(coinc_DOWN_CA3_353_q),length(coinc_DOWN_CA3_354_q),
                length(coinc_DOWN_CA3_355_q),length(coinc_DOWN_CA3_356_q),
                length(coinc_DOWN_CA3_357_q),length(coinc_DOWN_CA3_358_q))
five_down = c(length(coinc_DOWN_CA3_349_c),length(coinc_DOWN_CA3_350_c),length(coinc_DOWN_CA3_351_c),
            length(coinc_DOWN_CA3_352_c),length(coinc_DOWN_CA3_353_c),length(coinc_DOWN_CA3_354_c),
            length(coinc_DOWN_CA3_355_c),length(coinc_DOWN_CA3_356_c),length(coinc_DOWN_CA3_357_c),
            length(coinc_DOWN_CA3_358_c))
six_down = c(length(coinc_DOWN_CA3_349_s),length(coinc_DOWN_CA3_350_s),length(coinc_DOWN_CA3_351_s),
           length(coinc_DOWN_CA3_352_s),length(coinc_DOWN_CA3_353_s),length(coinc_DOWN_CA3_354_s),
           length(coinc_DOWN_CA3_355_s),length(coinc_DOWN_CA3_356_s),length(coinc_DOWN_CA3_357_s),
           length(coinc_DOWN_CA3_358_s))
tsummary_down = data.frame(unique_down,two_down,three_down,four_down,five_down,six_down)
colnames(tsummary_down) = c("Unique","Two types","Three types","Four types","Five types","Six types")
rownames(tsummary_down) = c("349_Mossy","350_Mossy","351_CA3-ve","352_CA3-ve","353_CA3-ve",
                          "354_CA3-ve","355_CA3-ve","356_CA3-do","357_CA3-do","358_CA3-do")

write.csv(x = genecount_df,file = "DE_Genes_Types_CA3.csv")

write.csv(gens_class_CA3_349,"List_Genes_Sig_349_Mossy.csv")
write.csv(gens_class_CA3_350,"List_Genes_Sig_350_Mossy.csv")
write.csv(gens_class_CA3_351,"List_Genes_Sig_351_CA3-ve.csv")
write.csv(gens_class_CA3_352,"List_Genes_Sig_352_CA3-ve.csv")
write.csv(gens_class_CA3_353,"List_Genes_Sig_353_CA3-ve.csv")
write.csv(gens_class_CA3_354,"List_Genes_Sig_354_CA3-ve.csv")
write.csv(gens_class_CA3_355,"List_Genes_Sig_355_CA3-ve.csv")
write.csv(gens_class_CA3_356,"List_Genes_Sig_356_CA3-do.csv")
write.csv(gens_class_CA3_357,"List_Genes_Sig_357_CA3-do.csv")
write.csv(gens_class_CA3_358,"List_Genes_Sig_358_CA3-do.csv")

write.csv(summary_table,"summary_table.csv")

write.csv(pattern_t_349_CA3,"pattern_Table_349_Mossy.csv")
write.csv(pattern_t_350_CA3,"pattern_Table_350_Mossy.csv")
write.csv(pattern_t_351_CA3,"pattern_Table_351_CA3-ve.csv")
write.csv(pattern_t_352_CA3,"pattern_Table_352_CA3-ve.csv")
write.csv(pattern_t_353_CA3,"pattern_Table_353_CA3-ve.csv")
write.csv(pattern_t_354_CA3,"pattern_Table_354_CA3-ve.csv")
write.csv(pattern_t_355_CA3,"pattern_Table_355_CA3-ve.csv")
write.csv(pattern_t_356_CA3,"pattern_Table_356_CA3-do.csv")
write.csv(pattern_t_357_CA3,"pattern_Table_357_CA3-do.csv")
write.csv(pattern_t_358_CA3,"pattern_Table_358_CA3-do.csv")

write.csv(llin_gens_up,"Genes_threshold_Total_UP.csv")
write.csv(llin_gens_down,"Genes_threshold_Total_DOWN.csv")

write.csv(coinc_UP_CA3_349_unique,"DEGsUP_uniques_349_Mossy.csv")
write.csv(coinc_UP_CA3_349_d,"DEGsUP_Shared2_349_Mossy.csv")
write.csv(coinc_UP_CA3_349_t,"DEGsUP_Shared3_349_Mossy.csv")
write.csv(coinc_UP_CA3_350_unique,"DEGsUP_uniques_350_Mossy.csv")
write.csv(coinc_UP_CA3_350_d,"DEGsUP_Shared2_350_Mossy.csv")
write.csv(coinc_UP_CA3_350_t,"DEGsUP_Shared3_350_Mossy.csv")
write.csv(coinc_UP_CA3_351_unique,"DEGsUP_uniques_351_CA3-ve.csv")
write.csv(coinc_UP_CA3_351_d,"DEGsUP_Shared2_351_CA3-ve.csv")
write.csv(coinc_UP_CA3_351_t,"DEGsUP_Shared3_351_CA3-ve.csv")
write.csv(coinc_UP_CA3_352_unique,"DEGsUP_uniques_352_CA3-ve.csv")
write.csv(coinc_UP_CA3_352_d,"DEGsUP_Shared2_352_CA3-ve.csv")
write.csv(coinc_UP_CA3_352_t,"DEGsUP_Shared3_352_CA3-ve.csv")
write.csv(coinc_UP_CA3_353_unique,"DEGsUP_uniques_353_CA3-ve.csv")
write.csv(coinc_UP_CA3_353_d,"DEGsUP_Shared2_353_CA3-ve.csv")
write.csv(coinc_UP_CA3_353_t,"DEGsUP_Shared3_353_CA3-ve.csv")
write.csv(coinc_UP_CA3_354_unique,"DEGsUP_uniques_354_CA3-ve.csv")
write.csv(coinc_UP_CA3_354_d,"DEGsUP_Shared2_354_CA3-ve.csv")
write.csv(coinc_UP_CA3_354_t,"DEGsUP_Shared3_354_CA3-ve.csv")
write.csv(coinc_UP_CA3_355_unique,"DEGsUP_uniques_355_CA3-ve.csv")
write.csv(coinc_UP_CA3_355_d,"DEGsUP_Shared2_355_CA3-ve.csv")
write.csv(coinc_UP_CA3_355_t,"DEGsUP_Shared3_355_CA3-ve.csv")
write.csv(coinc_UP_CA3_356_unique,"DEGsUP_uniques_356_CA3-do.csv")
write.csv(coinc_UP_CA3_356_d,"DEGsUP_Shared2_356_CA3-do.csv")
write.csv(coinc_UP_CA3_356_t,"DEGsUP_Shared3_356_CA3-do.csv")
write.csv(coinc_UP_CA3_357_unique,"DEGsUP_uniques_357_CA3-do.csv")
write.csv(coinc_UP_CA3_357_d,"DEGsUP_Shared2_357_CA3-do.csv")
write.csv(coinc_UP_CA3_357_t,"DEGsUP_Shared3_357_CA3-do.csv")
write.csv(coinc_UP_CA3_358_unique,"DEGsUP_uniques_358_CA3-do.csv")
write.csv(coinc_UP_CA3_358_d,"DEGsUP_Shared2_358_CA3-do.csv")
write.csv(coinc_UP_CA3_358_t,"DEGsUP_Shared3_358_CA3-do.csv")

write.csv(coinc_DOWN_CA3_349_unique,"DEGsDOWN_uniques_349_Mossy.csv")
write.csv(coinc_DOWN_CA3_349_d,"DEGsDOWN_Shared2_349_Mossy.csv")
write.csv(coinc_DOWN_CA3_349_t,"DEGsDOWN_Shared3_349_Mossy.csv")
write.csv(coinc_DOWN_CA3_350_unique,"DEGsDOWN_uniques_350_Mossy.csv")
write.csv(coinc_DOWN_CA3_350_d,"DEGsDOWN_Shared2_350_Mossy.csv")
write.csv(coinc_DOWN_CA3_350_t,"DEGsDOWN_Shared3_350_Mossy.csv")
write.csv(coinc_DOWN_CA3_351_unique,"DEGsDOWN_uniques_351_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_351_d,"DEGsDOWN_Shared2_351_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_351_t,"DEGsDOWN_Shared3_351_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_352_unique,"DEGsDOWN_uniques_352_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_352_d,"DEGsDOWN_Shared2_352_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_352_t,"DEGsDOWN_Shared3_352_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_353_unique,"DEGsDOWN_uniques_353_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_353_d,"DEGsDOWN_Shared2_353_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_353_t,"DEGsDOWN_Shared3_353_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_354_unique,"DEGsDOWN_uniques_354_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_354_d,"DEGsDOWN_Shared2_354_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_354_t,"DEGsDOWN_Shared3_354_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_355_unique,"DEGsDOWN_uniques_355_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_355_d,"DEGsDOWN_Shared2_355_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_355_t,"DEGsDOWN_Shared3_355_CA3-ve.csv")
write.csv(coinc_DOWN_CA3_356_unique,"DEGsDOWN_uniques_356_CA3-do.csv")
write.csv(coinc_DOWN_CA3_356_d,"DEGsDOWN_Shared2_356_CA3-do.csv")
write.csv(coinc_DOWN_CA3_356_t,"DEGsDOWN_Shared3_356_CA3-do.csv")
write.csv(coinc_DOWN_CA3_357_unique,"DEGsDOWN_uniques_357_CA3-do.csv")
write.csv(coinc_DOWN_CA3_357_d,"DEGsDOWN_Shared2_357_CA3-do.csv")
write.csv(coinc_DOWN_CA3_357_t,"DEGsDOWN_Shared3_357_CA3-do.csv")
write.csv(coinc_DOWN_CA3_358_unique,"DEGsDOWN_uniques_358_CA3-do.csv")
write.csv(coinc_DOWN_CA3_358_d,"DEGsDOWN_Shared2_358_CA3-do.csv")
write.csv(coinc_DOWN_CA3_358_t,"DEGsDOWN_Shared3_358_CA3-do.csv")

write.csv(tresum_up,"Summary_table_Coincidence_UP.csv")
write.csv(tresum_down,"Summary_table_Coincidence_DOWN.csv")
