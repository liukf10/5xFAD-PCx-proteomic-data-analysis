
#quantWorkFlow preprocess############ 
# extract data
input = "20210219_Maxquant_proteinGroups.txt"
output = "PCx"
blast_path = "ncbi-blast-2.7.1+/bin/"
db1_name = "human-uniprot-20181022.fasta"
db2_path = "mouse-all-20210202-APP added.fasta"
source("proteomic_data_extract_R")
#remove human APP expression
rm(list = ls())
load("PCx_LFQ.RData")
pos = which(rawdata$inf$ori.ID == "P05067"); #401
rawdata$inf <- rawdata$inf[-pos,];
rawdata$intensity <- rawdata$intensity[-pos,];
loginfo <- c(loginfo,remove_APP = TRUE)
save(rawdata,loginfo, file = "PCx_LFQ.RData")
#remove outlier and normalization
rm(list = ls())
input = "PCx_LFQ.RData"
output = "PCx"
funcdir = "D:/func"
step = 123
source("D:/func/step1_pro_R")
#covariance adjust
rm(list=ls())
quantfile="PCx_rmOut_impute.RData"
covar <- read.csv("5FAD_PCx_covar.csv");
covar <- covar[,-1]
output="PCx"
keepvar = "age,disease"
funcdir = "D:/func"
source("D:/func/step2_covariate_adj_R")

rm(list = ls())
library(DDPNA)
load("PCx_adj_final.RData")
load("inf.RData")
save(quant_final, covar_hid, covar,inf,
     file = paste0("PCx_adj_final.RData"))

#DDPNA############
#DEP 
rm(list = ls())
funcdir = "D:/func/"
input = "PCx_adj_final.RData"
output = "PCx"
source(paste0(funcdir,"DEP_func.R"))
load(input);
covar$diagnose <- paste0(covar$disease,covar$age)
library(openxlsx)
stat2way <- function(quant_final,covar, Egrp, Cgrp, 
                     student = TRUE, anova = TRUE){
  disease = covar$disease; age  = covar$age;
  fc <- fc_func(quant_final, covar$diagnose,Egrp,Cgrp);
  if(student) tp <- StudentT(quant_final,covar$diagnose,Egrp,Cgrp,var.equal = TRUE);
  if(!student) tp <- StudentT(quant_final,covar$diagnose,Egrp,Cgrp,var.equal = FALSE);
  if(anova) ap <- twoway_anova(quant = quant_final, group1 = disease,
                               group2 = age);
  if(anova)
    statinf <- data.frame(fc,ttestP = tp, ttestQ = p.adjust(tp,"fdr"),
                          ap, diseaseQ = p.adjust(ap[,1],"fdr"),
                          ageQ = p.adjust(ap[,2],"fdr"),
                          interactionQ = p.adjust(ap[,3],"fdr"))
  if(!anova) 
    statinf <- data.frame(fc,ttestP = tp, ttestQ = p.adjust(tp,"fdr"))
  statinf
}
info1 <- stat2way(quant_final,covar,Egrp = "AD3",Cgrp = "WT3",anova = FALSE)
info1 <- info1[,c(2,7,3,8,4,9,5,10,1,6,11,12,13)];
colnames(info1)[11:13] <- paste0("mo3_",c("log2fc","P","Q"));
info2 <- stat2way(quant_final,covar,Egrp = "AD4",Cgrp = "WT4",anova = FALSE)
info2 <- info2[,11:13]; colnames(info2)
colnames(info2)[1:3] <- paste0("mo4_",c("log2fc","P","Q"));
info3 <- stat2way(quant_final,covar,Egrp = "AD6",Cgrp = "WT6",anova = FALSE)
info3 <- info3[,11:13];colnames(info3)
colnames(info3)[1:3] <- paste0("mo6_",c("log2fc","P","Q"));
info4 <- stat2way(quant_final,covar,Egrp = "AD9",Cgrp = "WT9",anova = FALSE)
info4 <- info4[,11:13];colnames(info4)
colnames(info4)[1:3] <- paste0("mo9_",c("log2fc","P","Q"));
info5 <- stat2way(quant_final,covar,Egrp = "AD11",Cgrp = "WT11")
info5 <- info5[,11:19];colnames(info5)
colnames(info5)[1:3] <- paste0("mo11_",c("log2fc","P","Q"));
statinf <- data.frame(info1,info2,info3,info4,info5)
write.csv(statinf, file = paste0(output,"_stat.csv"))

#WGCNA
rm(list = ls())
library(WGCNA)
library(DDPNA)
input = "PCx_adj_final.RData";  
output ="PCx";  
power = 4; deepSplit = 2; minModuleSize = 30; minKMEtoStay = 0.3;
load(input);
data = quant_final; 
wgcnadata <- t(data)
CoExpNet = blockwiseModules(wgcnadata, power = power, maxBlockSize = 6000,
                            TOMType = "unsigned", deepSplit = deepSplit,
                            minModuleSize = minModuleSize,
                            minKMEtoStay = minKMEtoStay,
                            reassignThreshold = 0.05, mergeCutHeight = 0.15,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            verbose = 3)
moduleinf <- Module_inf(CoExpNet,inf)
moduleinf <- data.frame(moduleinf,GN = inf$GN)
MEinf <- ME_inf(CoExpNet$MEs,data,intensity.type = "none")
moduleinf = data.frame(moduleinf,inf[match(inf$ori.ID,moduleinf$ori.ID),5:6])
save(CoExpNet, MEinf, moduleinf, file = paste0(output,"_WGCNAnet.Rdata"))

#DDPNA-DEP-mod-enrich
rm(list = ls())
library(DDPNA)
output = "PCx"; 
funcdir = "D:/func/"
load(paste0(output,"_WGCNAnet.Rdata"))
stat <- read.csv(paste0(output,"_stat.csv"))
source(paste0(funcdir,"DEP_Mod_func.R"))
stat3 <- data.frame(ori.ID = stat$X,fc = 2^stat$mo3_log2fc, p = stat$mo3_P)
mo3 <- FCS(stat3,moduleinf,ordername = "fc",fc = 1.2)
stat4 <- data.frame(ori.ID = stat$X,fc = 2^stat$mo4_log2fc, p = stat$mo4_P)
mo4 <- FCS(stat4,moduleinf,ordername = "fc",fc = 1.2)
stat6 <- data.frame(ori.ID = stat$X,fc = 2^stat$mo6_log2fc, p = stat$mo6_P)
mo6 <- FCS(stat6,moduleinf,ordername = "fc",fc = 1.2)
stat9 <- data.frame(ori.ID = stat$X,fc = 2^stat$mo9_log2fc, p = stat$mo9_P)
mo9 <- FCS(stat9,moduleinf,ordername = "fc",fc = 1.2)
stat11 <- data.frame(ori.ID = stat$X,fc = 2^stat$mo11_log2fc, p = stat$mo11_P)
mo11 <- FCS(stat11,moduleinf,ordername = "fc",fc = 1.2)
save(mo3,mo4,mo6,mo9,mo11,file = paste0(output,"_DEP_ModEnrich.RData"))

#DDPNA-DAP
rm(list = ls())
library(DDPNA)
input = "PCx"; output = "PCx";
funcdir = "D:/func/"
source(paste0(funcdir,"DAP_func.R"))
load(paste0(input,"_DEP_ModEnrich.RData"))
load(paste0(input,"_WGCNAnet.Rdata"))
load(paste0(input,"_adj_final.RData"))
stat <- read.csv(paste0(input,"_stat.csv"))

DAP_3M <- DAPbyMod(quant_final, moduleinf, 13, mo3$DEPinf$ori.ID,
                   IDname="DEP_3M", cor.adj="fdr", coln = "ori.ID");

DAP_4M <- DAPbyMod(quant_final, moduleinf, 3, mo4$DEPinf$ori.ID,
                   IDname="DEP_4M", cor.adj="fdr", coln = "ori.ID");

DAP_6M <- DAPbyMod(quant_final, moduleinf, c(3,4), mo6$DEPinf$ori.ID,
                   IDname="DEP_6M",cor.adj="fdr", coln = "ori.ID");

DAP_9M <- DAPbyMod(quant_final, moduleinf, 3, mo9$DEPinf$ori.ID,
                   IDname="DEP_9M",cor.adj="fdr", coln = "ori.ID");

DAP_11M <- DAPbyMod(quant_final, moduleinf, c(2,3,13), mo11$DEPinf$ori.ID,
                    IDname="DEP_11M",cor.adj="fdr", coln = "ori.ID");

save(DAP_3M, DAP_4M,DAP_6M,DAP_9M,DAP_11M, file = paste0(output,"_DAP.RData"))

#DDPNA-DAPcore
rm(list = ls())
load("PCx_DAP.RData")
load("PCx_adj_final.RData")
DAPcore <- intersect(DAP_4M,DAP_6M);
DAPcore <- intersect(DAPcore, DAP_9M); DAPcore <- intersect(DAPcore, DAP_11M);
DAPcoreGN <- unique(inf$GNmouse[inf$ori.ID %in% DAPcore])
write.table(DAPcoreGN, file = "PCx_DAPcore.txt", quote = F, row.names = F, col.names = F)

#DDPNA-DAPnet
rm(list = ls())
library(DDPNA)
input = "PCx"; output = "PCx";
funcdir = "D:/func/"
source(paste0(funcdir,"DAP_func.R"))
load(paste0(input,"_DEP_ModEnrich.RData"))
load(paste0(input,"_WGCNAnet.Rdata"))
load(paste0(input,"_adj_final.RData"))
stat <- read.csv(paste0(input,"_stat.csv"))

netDAP_4M <- DAPextract2(quant_final, moduleinf, 3, mo4$DEPinf$ori.ID,
                         IDname="DEP_4M",cor.adj="fdr",coln = "ori.ID");
netDAP_6M <- DAPextract2(quant_final, moduleinf, 3, mo6$DEPinf$ori.ID,
                         IDname="DEP_6M",cor.adj="fdr",coln = "ori.ID");
netDAP_9M <- DAPextract2(quant_final, moduleinf, 3, mo9$DEPinf$ori.ID,
                         IDname="DEP_9M",cor.adj="fdr",coln = "ori.ID");
netDAP_11M <- DAPextract2(quant_final, moduleinf, 3, mo11$DEPinf$ori.ID,
                          IDname="DEP_11M",cor.adj="fdr",coln = "ori.ID");
save(netDAP_4M,netDAP_6M,netDAP_9M,netDAP_11M,file=paste0(output,"_DAPnet.RData"))

#enrichment #####
#DEP enrichment
rm(list = ls())
load("PCx_DEP_ModEnrich.RData")
load("PCx_WGCNAnet.Rdata")
library(clusterProfiler)
library(org.Mm.eg.db)
x = moduleinf$GNmouse[moduleinf$ori.ID %in% mo3$DEPinf$ori.ID]
ego_DEP3= enrichGO(x,"org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",
                   qvalueCutoff = 1,pvalueCutoff = 0.05)
x = moduleinf$GNmouse[moduleinf$ori.ID %in% mo4$DEPinf$ori.ID]
ego_DEP4= enrichGO(x,"org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",
                   qvalueCutoff = 1,pvalueCutoff = 0.05)
x = moduleinf$GNmouse[moduleinf$ori.ID %in% mo6$DEPinf$ori.ID]
ego_DEP6= enrichGO(x,"org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",
                   qvalueCutoff = 1,pvalueCutoff = 0.05)
x = moduleinf$GNmouse[moduleinf$ori.ID %in% mo9$DEPinf$ori.ID]
ego_DEP9= enrichGO(x,"org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",
                   qvalueCutoff = 1,pvalueCutoff = 0.05)
x = moduleinf$GNmouse[moduleinf$ori.ID %in% mo11$DEPinf$ori.ID]
ego_DEP11= enrichGO(x,"org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",
                    qvalueCutoff = 1,pvalueCutoff = 0.05)
save(ego_DEP3,ego_DEP4,ego_DEP6,ego_DEP9,ego_DEP11, file = "ego_DEP.RData")

#DAPcore and M3 enrichment
rm(list = ls())
load("PCx_DAPnet.RData")
load("moduleinf.RData")
gene <- intersect(netDAP_4M$netgene,netDAP_6M$netgene)
gene <- intersect(gene,netDAP_9M$netgene)
gene <- intersect(gene,netDAP_11M$netgene)
x <- moduleinf$GNmouse[moduleinf$ori.ID %in% gene]
library(clusterProfiler)
ego_DAPcore2= enrichGO(x,"org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",
                       qvalueCutoff = 1,pvalueCutoff = 0.05, universe = moduleinf$GNmouse)
save(ego_DAPcore2,file = "ego_DAPcore.RData")
x<- moduleinf$GNmouse[moduleinf$moduleNum == 3]
ego_M3_2= enrichGO(x,"org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",
                   qvalueCutoff = 1,pvalueCutoff = 0.05, universe = moduleinf$GNmouse)
save(ego_M3_2,file = "egoM3.RData")

##AD_LFQ_adj_batch_remove#########
rm(list = ls())
input = c("ACT", "Banner","MSBB", "BLSA_FC","BLSA_PC")
for ( i in input) {
  load(paste0(gsub("_(PC|FC)","",i),"/",i,"_adj_final.RData"))
  if(i == input[1]) {
    quant = matrix(nrow = nrow(quant_final),ncol = 0)
    meta_all = matrix(nrow = 0,ncol = ncol(covar))
    inf <- data_rmOut_imp$inf;} 
  pos <- match(inf$GN,data_rmOut_imp$inf$GN)
  quant <- quant[!is.na(pos),];
  inf <- inf[!is.na(pos),];
  quant <- cbind(quant,quant_final[pos[!is.na(pos)],]);
  quant <- quant[!duplicated(inf$GN),];
  inf <- inf[!duplicated(inf$GN),];
  covar <- data.frame(SampleName = rownames(covar), 
                      diagnose = covar$diagnose, batch = i)
  covar$diagnose[covar$diagnose == "asym"] <- "MCI";
  covar$diagnose[covar$diagnose == "ad"] <- "AD";
  covar$diagnose[covar$diagnose == "ctl"] <- "CTL";
  meta_all <- rbind(meta_all,covar);
  print(nrow(quant))
}

#ACT 3438   3284   #3103  #3057 #3022 #3022
covar <- meta_all; 
loginfo <- list(mergedata = list(dataset = input,type = "LFQ",time = Sys.time()))
save(covar,quant,inf,loginfo,file = "AD_LFQ_adj_merge.RData")

rm(list = ls())
quantfile = "AD_LFQ_adj_merge.RData";
output = "AD_LFQ_adj"
funcdir = "D:/func"
keepvar = "diagnose"
source("D:/func/v1.2/step2_batch_adj_R")

###AD_LFQ_5FAD_feature.RData###############
rm(list = ls())
load("PCx_adj_final.RData")
covar$diagnose <- paste0(covar$age, covar$disease)
PCx <- list(covar = covar, inf = inf, quant = quant_final)
rm(covar, inf, quant_final)
load("AD_LFQ_adj_batch_remove.RData")
ADlfq <- list(covar = covar, inf = inf, quant = quant)
rm(covar, inf, quant)

pos <- match(PCx$inf$GN,ADlfq$inf$GN)
quant <- PCx$quant[!is.na(pos),]
inf <- PCx$inf[!is.na(pos),];
quant <- cbind(quant,ADlfq$quant[pos[!is.na(pos)],]);
inf <- cbind(inf, ori.ID_hs=ADlfq$inf[1][pos[!is.na(pos)],])
inf$GN[duplicated(inf$GN)]
quant <- quant[!duplicated(inf$GN),]
inf <- inf[!duplicated(inf$GN),]

meta_all <- data.frame(SampleName = rownames(PCx$covar), 
                       diagnose = PCx$covar$diagnose,
                       batch = "5xFAD")	
meta_all <- rbind(meta_all ,
                  data.frame(SampleName = rownames(ADlfq$covar), 
                             diagnose = ADlfq$covar$diagnose,
                             batch = ADlfq$covar$batch));
covar = meta_all
save(covar,quant, inf, file = "AD_LFQ_5FAD_adj_merge.RData")

rm(list = ls())
quantfile = "AD_LFQ_5FAD_adj_merge.RData";
output = "AD_LFQ_5FAD_adj"
funcdir = "D:/func/v1.1"
source("D:/func/v1.2/step2_batch_adj_R")

load("AD_LFQ_5FAD_adj_batch_remove.RData")
oneway_anova <- function(quant, group, NAnum = NA,
                         var.equal = NA) {
  if(!is.na(NAnum)) quant[quant == NAnum] <- NA;
  if(is.na(var.equal))p <- anova_p(quant,group);
  if(isTRUE(var.equal)) {
    p = NULL
    for(i in 1:nrow(quant)) {
      p = c(p,oneway.test(unlist(quant[i,]) ~ group)$p.value);
      if(i == trunc(0.25*nrow(quant)) & i > 2000)
        print("....25% completed ......")
      if(i == trunc(0.5*nrow(quant)) & i > 4000)
        print("....50% completed ......")  
      if(i == trunc(0.75*nrow(quant)) & i > 4000)
        print("....75% completed ......")
    }
  }
  if(isFALSE(var.equal)) {
    p = NULL
    for(i in 1:nrow(quant)) {
      p = c(p,oneway.test(unlist(quant[i,]) ~ group,var.equal = FALSE)$p.value);
      if(i == trunc(0.25*nrow(quant)) & i > 2000)
        print("....25% completed ......")
      if(i == trunc(0.5*nrow(quant)) & i > 4000)
        print("....50% completed ......")  
      if(i == trunc(0.75*nrow(quant)) & i > 4000)
        print("....75% completed ......")
    }
  }
  print("....100% completed ......")
  p
}
group = covar$diagnose[covar$batch != "5xFAD"]
table(group); group <- factor(group,levels = c("CTL","MCI","AD"))
quant2 = quant[,covar$batch != "5xFAD"]
p <- oneway_anova(quant2, group)
pos2 = p.adjust(p,method = "fdr") <0.05
quant_filter_aov <- quant[pos2,]
save(quant_filter_aov,covar,
     file = "AD_LFQ_5FAD_feature.RData")
