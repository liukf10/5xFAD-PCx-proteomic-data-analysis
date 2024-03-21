##Fig.2b PCA######
rm(list = ls())
load("PCx_adj_final.RData")

library(ggprism)
library(ggplot2)
data.pca <- prcomp(t(quant_final),scale = TRUE)
importance <- summary(data.pca)$importance;
p1 <- as.numeric(sprintf("%.3f",importance[2,1]))*100;
p2 <- as.numeric(sprintf("%.3f",importance[2,2]))*100;
pcadata <- as.data.frame(data.pca$x);

plotdata <- data.frame(pcadata[,1:2], covar)
plotdata$disease <- factor(plotdata$disease, levels = c("WT","AD"),labels = c("WT","5xFAD"))
plotdata$age <- as.factor(plotdata$age)
colnames(plotdata)[4] = "group"
values = c("#00A600FF","#E6E600FF","#EDB48EFF","#FF6000FF","#FF0000FF")

ggplot(plotdata,aes(PC1,PC2)) + 
  theme_bw() +
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  geom_point(aes(color = age, shape = group),size = 3) +
  scale_color_manual(values =values) +
  labs(x=paste0("PC1(",p1,"%)"),y=paste0("PC2(",p2,"%)"))

##Fig.2c ########
rm(list = ls())
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)
#Core DEPs were manual selected accroding to our work and also fucntion of these proteins
#Core DEPs includes unique DEP in different grounp and also overlapping DEPs of different group
Core <- read.csv("core.csv",header = T,fill=NA)
Core$Group <- as.factor(Core$Group)
Core$Symbol <- factor(Core$Symbol,
                      levels = c(unique(Core$Symbol)))

p <-ggplot(data = Core,mapping = aes(x=Group,y=Symbol))
p+geom_point(aes(size=Log2FC,color=p))+theme_bw()+
  scale_colour_gradient(low="#ff0000", high="#ffcccc")+
  theme(panel.grid=element_blank())+
  theme(text = element_text(size = 20),axis.text = element_text(size = 20))+
  labs(x=NULL,y=NULL)
##Fig.2d #######
rm(list = ls())
load('ego_DEP.RData')
load("compareCluster_template.RData")
xx@compareClusterResult = rbind(data.frame(Cluster = "DEP03M",as.data.frame(ego_DEP3)),
                                data.frame(Cluster = "DEP04M",as.data.frame(ego_DEP4)),
                                data.frame(Cluster = "DEP06M",as.data.frame(ego_DEP6)),
                                data.frame(Cluster = "DEP09M",as.data.frame(ego_DEP9)),
                                data.frame(Cluster = "DEP11M",as.data.frame(ego_DEP11)))
xx@compareClusterResult$Cluster <- as.factor(xx@compareClusterResult$Cluster)
xx@geneClusters <- list(ego_DEP3@gene,ego_DEP4@gene,ego_DEP6@gene,ego_DEP9@gene,ego_DEP11@gene)
names(xx@geneClusters) <- c("DEP03M","DEP04M","DEP06M","DEP09M",'DEP11M')

dotplot(xx,label_format = 65,font.size = 16)
##Supplementary Fig.3 ########
rm(list = ls())
library(RRHO)
stat <- read.csv("PCx_stat.csv")

source("trend_func.R")

mo3 <- stat[,c("X","mo3_log2fc")];
mo4 <- stat[,c("X","mo4_log2fc")];
mo6 <- stat[,c("X","mo6_log2fc")];
mo9 <- stat[,c("X","mo9_log2fc")];
mo11 <- stat[,c("X","mo11_log2fc")];
rrho3vs4 <- RRHO(mo3,mo4,alternative = "enrichment",
                 labels = c("mo3","mo4"))
rrho4vs6 <- RRHO(mo4,mo6,alternative = "enrichment",
                 labels = c("mo4","mo6"))
rrho6vs9 <- RRHO(mo6,mo9,alternative = "enrichment",
                 labels = c("mo6","mo9"))
rrho9vs11 <- RRHO(mo9,mo11,alternative = "enrichment",
                  labels = c("mo9","mo11"))
library(ggheatmap)
rrho_plot(rrho3vs4,labels = c("mo3","mo4"), color_range = c(0,110))
rrho_plot(rrho4vs6,labels = c("mo4","mo6"), color_range = c(0,110))
rrho_plot(rrho6vs9,labels = c("mo6","mo9"), color_range = c(0,110))
rrho_plot(rrho9vs11,labels = c("mo9","mo11"), color_range = c(0,110))

####Fig.3a-b  WGCNA dendrogram##########
rm(list = ls())
input = "PCx"; output = "PCx"; 
metaflie = "5FAD_PCx_covar.csv";

library(WGCNA)
load(paste0(input,"_WGCNAnet.Rdata"))
#meta info extract
meta <- read.csv(metaflie);
for (i in 1:ncol(meta)) {
  if(all(rownames(MEinf) %in% meta[,i]))
    break
  if(i == ncol(meta)) i = 0;
}
pos <- match(rownames(MEinf),meta[,i]);
meta <- meta[pos,-i]; #meta & MEs are matched
rownames(meta) <- rownames(MEinf);
if(ncol(meta)!=1) 
  meta <- meta[,c(which(colnames(meta) == "diagnose"),which(colnames(meta) != "diagnose"))]
#plot1
moduleLabels = CoExpNet$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(CoExpNet$dendrograms[[1]], 
                    moduleColors[CoExpNet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#plot2
#Module and meta relationship
MEs =  orderMEs(MEinf)
for(i in 1:ncol(meta))  meta[,i]<-as.numeric(factor(meta[,i]));
modmodCor = cor(MEs, meta, use = "p")
modmodP = corPvalueStudent(modmodCor, nrow(MEs));
modPadj = matrix(p.adjust(modmodP,method = "fdr"),
                 nrow = nrow(modmodP), ncol = ncol(modmodP));
rownames(modPadj) <- rownames(modmodP); colnames(modPadj) <- colnames(modmodP);

textMatrix = paste(signif(modmodCor, 2), "\n(", signif(modPadj, 1), ")", sep = "")
dim(textMatrix) = dim(modmodCor)
textMatrix[modPadj>0.05]<-""
par(mar = c(5, 6, 3, 4));
labeledHeatmap(Matrix = t(modmodCor), xLabels = colnames(MEs), 
               yLabels = colnames(meta), 
               cex.lab = 0.7, xLabelsAngle = 0, xLabelsAdj = 0.5,
               ySymbols = colnames(meta), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = t(textMatrix), setStdMargins = FALSE, 
               cex.text = 0.6, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

####Fig.3c part#########
rm(list = ls())
library(DDPNA)
load("PCx_DEP_ModEnrich.RData")
load("PCx_WGCNAnet.Rdata")
source("DEP_Mod_func.R")
MO3 <- enrichmod(moduleinf, mo3, coln = "ori.ID", FCS = F)
MO4 <- enrichmod(moduleinf, mo4, coln = "ori.ID", FCS = F)
MO6 <- enrichmod(moduleinf, mo6, coln = "ori.ID", FCS = F)
MO9 <- enrichmod(moduleinf, mo9, coln = "ori.ID", FCS = F)
MO11 <- enrichmod(moduleinf, mo11, coln = "ori.ID", FCS = F)

MO3_up <- enrichmod(moduleinf, mo3, IDtype = "up",coln = "ori.ID", FCS = F)
MO4_up <- enrichmod(moduleinf, mo4, IDtype = "up", coln = "ori.ID", FCS = F)
MO6_up <- enrichmod(moduleinf, mo6, IDtype = "up", coln = "ori.ID", FCS = F)
MO9_up <- enrichmod(moduleinf, mo9, IDtype = "up", coln = "ori.ID", FCS = F)
MO11_up <- enrichmod(moduleinf, mo11, IDtype = "up", coln = "ori.ID", FCS = F)

MO3_down <- enrichmod(moduleinf, mo3, IDtype = "down", coln = "ori.ID", FCS = F)
MO4_down <- enrichmod(moduleinf, mo4, IDtype = "down", coln = "ori.ID", FCS = F)
MO6_down <- enrichmod(moduleinf, mo6, IDtype = "down", coln = "ori.ID", FCS = F)
MO9_down <- enrichmod(moduleinf, mo9, IDtype = "down", coln = "ori.ID", FCS = F)
MO11_down <- enrichmod(moduleinf, mo11, IDtype = "down", coln = "ori.ID", FCS = F)

DEP_Modall <- list('3M' = MO3, '3M_up' = MO3_up,'3M_down' = MO3_down,
                   '4M' = MO4, '4M_up' = MO4_up,'4M_down' = MO4_down,
                   '6M' = MO6, '6M_up' = MO6_up,'6M_down' = MO6_down,
                   '9M' = MO9, '9M_up' = MO9_up,'9M_down' = MO9_down, 
                   '11M' = MO11, '11M_up' = MO11_up,'11M_down' = MO11_down)

heatMapInfall <- DEP_Mod_HeatMap(DEP_Modall, xlab = "Mod",filter = "p",
                                 xLabelsAngle = 0);



##Fig.3d ####
rm(list = ls())
library(reshape2)
library(ggplot2)
library("cowplot")
data <- read.csv("cell type enrichment.csv",header = T, stringsAsFactors = F)
rownames(data) <- data$X
data <- data[,-1]
data2 <- -log10(data)
dep4 <- data2[2,];dep6 <- data2[3,];dep9 <- data2[4,];dep11 <- data2[5,]
dap4 <- data2[7,];dap6 <- data2[8,];dap9 <- data2[9,];dap11 <- data2[10,]
dep4 <- melt(dep4);dep6 <- melt(dep6);dep9 <- melt(dep9);dep11 <- melt(dep11) 
dap4 <- melt(dap4);dap6 <- melt(dap6);dap9 <- melt(dap9);dap11 <- melt(dap11) 

pdep4 <- ggplot(data = dep4, aes(x=variable,y=value,fill=variable))
dep.4 <- pdep4 +geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

pdep6 <- ggplot(data = dep6, aes(x=variable,y=value,fill=variable))
dep.6 <- pdep6+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

pdep9 <- ggplot(data = dep9, aes(x=variable,y=value,fill=variable))
dep.9 <- pdep9+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

pdep11 <- ggplot(data = dep11, aes(x=variable,y=value,fill=variable))
dep.11 <- pdep11+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

pdap4 <- ggplot(data = dap4, aes(x=variable,y=value,fill=variable))
dap.4 <- pdap4+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

pdap6 <- ggplot(data = dap6, aes(x=variable,y=value,fill=variable))
dap.6 <- pdap6+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

pdap9 <- ggplot(data = dap9, aes(x=variable,y=value,fill=variable))
dap.9 <- pdap9+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

pdap11<- ggplot(data = dap11, aes(x=variable,y=value,fill=variable))
dap.11 <- pdap11+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

plot_grid(dep.4,dap.4,dep.6,dap.6,dep.9,dap.9,dep.11,dap.11,nrow = 4,ncol = 2)

##Fig.3e-f #########
rm(list = ls())
load("ego_DAPcore.RData")
ego_DAPcore2@result <- ego_DAPcore2@result[c(1,3,4,7,18,19,20,21,26,73),]
ego_DAPcore2@pvalueCutoff <- 0.1
dotplot(ego_DAPcore2,label_format = 40)

load("egoM3.RData")
ego_M3_2@result <- ego_M3_2@result[c(1,3,4,5,11,16,18,29,36,49),]
barplot(ego_M3_2,label_format = 45,showCategory = 10)

##Fig.5a heatmap#########
rm(list = ls())
load("PCx_DAP.RData")
load("moduleinf.RData")
load("PCx_adj_final.RData")
covar$diagnose <- paste0(covar$age,covar$disease)
mean <- apply(quant_final,1,function(x) tapply(x, covar$diagnose, mean))
mean <-  data.frame(t(mean));
mean <- mean[,c(4,3,6,5,8,7,10,9,2,1)]
colnames(mean) = c("WT3M","AD3M","WT4M","AD4M","WT6M","AD6M","WT9M","AD9M","WT11M","AD11M")
GN = moduleinf$GNmouse
while(sum(duplicated(GN))>0)
  GN[duplicated(GN)] = paste0(GN[duplicated(GN)]," ")
rownames(mean) = GN

library(Mfuzz)
library(pheatmap)
x = moduleinf$GNmouse[moduleinf$ori.ID %in% DAP_4M & moduleinf$celltype == "Astro"]
pick <- mean[rownames(mean)%in%x,]
pick <- new('ExpressionSet',exprs = as.matrix(pick))
pick <- standardise(pick)@assayData$exprs
pheatmap(pick,cluster_cols = FALSE)


x = GN[moduleinf$ori.ID %in% DAP_4M & moduleinf$celltype == "Micro"]
pick <- mean[rownames(mean)%in%x,]
pick <- new('ExpressionSet',exprs = as.matrix(pick))
pick <- standardise(pick)@assayData$exprs
pheatmap(pick,cluster_cols = FALSE)

## Fig.5d  ###########
rm(list = ls())
load("AD_LFQ_5FAD_feature.RData")
quant <- quant_filter_aov
pca <- function(exprM, scale = TRUE,...){
  sum <- apply(exprM,1,sum);
  exprM <- exprM[!is.na(sum),];
  data <- t(as.matrix(exprM));
  # do PCA 
  data.pca <- prcomp(data,scale = scale,...)
  # fetch the proportion of PC1 and PC2
  pc <- as.data.frame(data.pca$x);
}
pcadata <- pca(quant)
pcadata <- data.frame(pcadata, covar)
pcadata <- data.frame(pcadata[,1:2], covar)

plotdata <- pcadata
plotdata <- rbind(plotdata[plotdata$diagnose == "CTL",],
                  plotdata[plotdata$diagnose == "MCI",],
                  plotdata[plotdata$diagnose == "AD",],
                  plotdata[plotdata$batch == "5xFAD",])
plotdata$group <- plotdata$batch
plotdata$group[grepl("WT",plotdata$diagnose)] = "WT"
plotdata$group[grepl("AD",plotdata$diagnose)] = "5xFAD"
plotdata$group[plotdata$batch != "5xFAD"] = "human"
plotdata$group <- factor(plotdata$group, levels = c("human","WT","5xFAD"))
plotdata$size <- 1
plotdata$size[plotdata$batch == "5xFAD"] <- 2
plotdata$mouse <- plotdata$diagnose
plotdata$diagnose[plotdata$batch == "5xFAD"] <- gsub("[A-Z]+","",plotdata$diagnose[plotdata$batch == "5xFAD"]) 
plotdata$diagnose <- factor(plotdata$diagnose, levels = c("CTL",'MCI','AD',"3","4","6","9","11"))

library(ggprism)
plotdata$y <- NA
plotdata$y[plotdata$group =="WT"]  <- 0.02
plotdata$y[plotdata$group =="5xFAD"] <- 0.03
levels(plotdata$diagnose) <- c("CTL","MCI","SAD","03","04","06","09","11")

#seperate mouse data
plotdata2 <- pcadata[pcadata$batch == "5xFAD",]
plotdata2$group <- gsub("[0-9]+","",plotdata2$diagnose)
plotdata2$age <- gsub("[A-Z]+","",plotdata2$diagnose)
plotdata2$age <- factor(as.numeric(plotdata2$age))
plotdata2$y <- -0.01 #0.02
plotdata2$y[plotdata2$group =="AD"] <- -0.005 #0.03
plotdata2$group <- factor(plotdata2$group,levels = c("WT","AD"), labels = c("WT","5xFAD"))

ggplot(plotdata[plotdata$batch != "5xFAD"&plotdata$PC1>-20,], aes(x=PC1,fill=diagnose)) +
  geom_density(alpha=.25) +#facet_wrap(~diagnose,ncol = 1) +
  theme_prism() +
  geom_point(mapping=aes(x = PC1, y = y, shape = group, color = age),fill=NA,
             data = plotdata2,size = 1) +
  scale_color_manual(values = c("#00A600FF","#E6E600FF","#EDB48EFF","#FF6000FF","#FF0000FF")) +
  scale_fill_manual(values = c("#FAEBD755","#00F5FF55","#008B8B55"),
                    labels = c("CTL","MCI","AD")) +
  ylab("density")

## Fig.5e-g #########
rm(list = ls())
load("ROSMAP_impute_DAP4M.RData")
library("Hmisc")
library(ggplot2)
library(ggprism)
gn <- unique(info$GN[duplicated(info$GN)])
pos = NULL
for(i in gn) {
  keep <- which(info$ori == gene$new.ID[gene$GN == i]);
  removepos <- which(info$GN %in% i)
  if(length(keep) == 0) 
    removepos <- removepos[-1] else  
      removepos <- removepos[removepos!=keep]
    pos <- c(pos,removepos)
}
info <- info[-pos,]
quant <- quant[-pos,]
ast <- info$ori.ID[info$GN %in% gene$GN[gene$celltype =="Astro"]]
mic <- info$ori.ID[info$GN %in% gene$GN[gene$celltype =="Micro"]]

sig_cor<-function(cor_matrix,p=0.05,r=NA,type=1){
  cor_matrix$P[upper.tri(cor_matrix$P,diag = TRUE)]<-NA;
  cor_matrix$r[upper.tri(cor_matrix$P,diag = TRUE)]<-NA;
  if(!is.na(p)){
    logi<-cor_matrix$P>p;
    cor_matrix$P[logi]<-NA;
    cor_matrix$r[logi]<-NA;
  }
  if(!is.na(r)){
    logi2<-abs(cor_matrix$r)<r
    cor_matrix$P[logi2]<-NA;
    cor_matrix$r[logi2]<-NA;
  }
  if(!is.na(p)|!is.na(r)){
    posit<-which(!is.na(cor_matrix$P), arr.ind=TRUE);
    if(type==2){
      #对这些位置排序，提取出方阵
      column<-sort(unique(c(posit[,1],posit[,2])));
      cor_matrix$P<-cor_matrix$P[column,column];
      cor_matrix$r<-cor_matrix$r[column,column];
    }
    else if(type==1){
      #对这些位置排序提取矩阵
      roworder<-sort(unique(posit[,1]));
      column<-sort(unique(posit[,2]));
      cor_matrix$P<-cor_matrix$P[roworder,column];
      cor_matrix$r<-cor_matrix$r[roworder,column];
    }
  }
  cor_matrix
}
orderCor <- function(cormat, pmat) {
  postion<-which(!is.na(cormat),arr.ind = TRUE);
  pos<-which(!is.na(cormat));
  ro<-attr(cormat,"dimnames")[[1]][postion[,1]];
  co<-attr(cormat,"dimnames")[[2]][postion[,2]];
  r<-cormat[pos];
  p<-pmat[pos];
  flat<-data.frame(row=ro,column=co,cor=r,p=p);
  orderCor<-flat[order(abs(flat$cor),decreasing = TRUE),];
  rownames(orderCor)<-order(abs(orderCor$cor),decreasing = TRUE);
  orderCor
}
corlist <- function(group = "diagnose", name, times = 10,quant) {
  n <- min(table(metainfo[,group])[table(metainfo[,group])>9])
  n =80
  quant2 <- t(quant[,metainfo[,group] == name])
  if(nrow(quant2) > n) {
    for(j in 1:times) {
      set.seed(j)
      pickSample <- sample(1:nrow(quant2),n)
      quant3 <- quant2[pickSample,]
      x1 <- rcorr(quant3,type = "pearson")
      x1 <-sig_cor(x1,p=1,r=NA,type = 2)
      x1<-orderCor(x1$r,x1$P)
      for( i in which(x1$p == 0)){
        fit <- lm(quant3[,x1$row[i]]~quant3[,x1$column[i]]);
        tstats <- coef(fit)/sqrt(diag(vcov(fit)));
        x1$p[i] <- 2*pt(abs(tstats),df=df.residual(fit),lower.tail = FALSE)[2]
      }
      if(j == 1) {x <- x1; r1 <- x1$cor} else r1 <- r1+x1$cor
    }
    x$cor <- r1/times
    x$p <- 2 * (1 - pt(abs(x$cor) * sqrt(80-2)/sqrt(pmax(1 - x$cor * x$cor, 0)), 80-2))
  } else {
    x <- rcorr(quant2,type = "pearson")
    x <-sig_cor(x,p=1,r=NA,type = 2)
    x<-orderCor(x$r,x$P)
    for( i in which(x$p == 0)){
      fit <- lm(quant2[,x$row[i]]~quant2[,x$column[i]]);
      tstats <- coef(fit)/sqrt(diag(vcov(fit)));
      x$p[i] <- 2*pt(abs(tstats),df=df.residual(fit),lower.tail = FALSE)[2]
    }
  }
  x$p.adj <- p.adjust(x$p,method = "fdr") 
  x
}
linkfilter <- function(x){
  pos1 = x[,1] %in% ast & x[,2] %in% mic
  pos2 = x[,1] %in% mic & x[,2] %in% ast
  x[pos1|pos2,]
}
ctl <- corlist(name = "ctl",times=10,quant=quant)
mci <- corlist(name = "mci",times=10,quant=quant)
ad <- corlist(name = "ad",times=10,quant=quant)
ctl_c <- linkfilter(ctl);
mci_c <- linkfilter(mci);
ad_c <- linkfilter(ad);
#Fig.5g write graphml file and then used Cytoscape to modifiy
{
mergelink <- function(x) {
  net <- x[[1]]
  net$cor[net$p > 0.05] <- NA
  colnames(net)[3] <- names(x)[1];
  net <- net[,-c(4:ncol(net))]
  POS <- is.na(net[,3])
  for(i in 2:length(x)) {
    net$x <- NA; colnames(net)[i+2] <- names(x)[i]
    pos = x[[i]]$p < 0.05
    pos1 = paste0(x[[i]][pos,1],x[[i]][pos,2])
    net[match(pos1, paste0(net[,1],net[,2])),i+2] <- x[[i]][pos,3]
    POS <- c(POS&is.na(net[,i+2]))
  }
  net <- net[!POS,]
  net[,1] <- info$GN[match(net[,1],info$ori)]
  net[,2] <- info$GN[match(net[,2],info$ori)]
  net
}
cellinteraction <- mergelink(list(ad=ad_c,mci=mci_c,ctl = ctl_c))
library(igraph)
link <- cellinteraction
g <- graph.data.frame(link)
netgene <- names(V(g)); 
genecolor <- gene$celltype[match(netgene,gene$GN)]
genecolor[genecolor== ""] = "gray"
genecolor[genecolor== "Astro"] = "cyan"
genecolor[genecolor== "Micro"] = "green"
genecolor[genecolor== "Oligo"] = "yellow"
genecolor[genecolor== "Neuro"] = "red"
V(g)$color <- genecolor
linkcolor <- rep("black",nrow(link))
linkcolor[is.na(link$mci)&is.na(link$ctl)] <- "red"
linkcolor[is.na(link$ad)&is.na(link$ctl)] <- "#FF4500"
linkcolor[is.na(link$ad)&is.na(link$mci)] <- "yellow"
linkcolor[(!is.na(link$ad))&(!is.na(link$mci))] <- "gray11"
linkcolor[(!is.na(link$mci))&(!is.na(link$ctl))] <- "gray22"
linkcolor[(!is.na(link$ad))&(!is.na(link$ctl))] <- "darkblue"
table(linkcolor)
E(g)$color <- linkcolor
E(g)$lwd <- apply(link[,3:5],1,mean,na.rm=TRUE)
#write graphml file and then used Cytoscape to modifiy
#write_graph(g,"Fig.5g.graphml","graphml")
}

#Fig.5e
plotdfunc <- function(x,name = "cor"){
  plotdata <- NULL
  for( i in 1:length(x)) {
    plotdata <- rbind(plotdata,data.frame(group = names(x)[i],value = abs(x[[i]][,name])))
  }
  plotdata$group <- factor(plotdata$group,levels = c("ctl","mci","ad"));
  p1 <- t.test(value~group,plotdata,group %in% c("ctl","mci"))$p.value
  p2 <- t.test(value~group,plotdata,group %in% c("ctl","ad"))$p.value
  p3 <- t.test(value~group,plotdata,group %in% c("mci","ad"))$p.value
  print(paste0("t.test:Pctlvsmci: ",p1,"; Pctlvsad: ",p2,"; Pmcivsad: ",p3))
  p1 <- wilcox.test(value~group,plotdata,group %in% c("ctl","mci"))$p.value
  p2 <- wilcox.test(value~group,plotdata,group %in% c("ctl","ad"))$p.value
  p3 <- wilcox.test(value~group,plotdata,group %in% c("mci","ad"))$p.value
  print(paste0("wilcox:Pctlvsmci: ",p1,"; Pctlvsad: ",p2,"; Pmcivsad: ",p3))
  plotdata
}
VBplot_1way <- function(data_m, coln ,violin = FALSE,boxplot = TRUE,
                        outlier = TRUE,point = FALSE,
                        meanpoint = TRUE,prism = TRUE, prismpalette = "floral",
                        vcolor = TRUE, vfill = FALSE, scaleviolin = FALSE,
                        bcolor = TRUE, bfill = FALSE, 
                        title = NULL, y = NULL, boxw = 0.2, 
                        dotsize = 1, binwidth = 1.5,
                        ylim, ...){
  variable = factor(data_m[,1]);
  data_m <- data.frame(variable = variable, value = data_m[,coln]);
  if(length(vcolor) > 1 | length(bcolor) > 1) prismpalette = FALSE
  plotparameter = list(vcolor=vcolor,vfill=vfill,bcolor=bcolor,bfill=bfill)
  for(i in seq_along(plotparameter)) {
    x = plotparameter[[i]]
    if(isTRUE(x)) 
      x = variable else if(length(x) == nlevels(variable))
        x = as.character(factor(variable,labels = x)) else if(length(x) == length(variable))
          x = x else x = 1;
          plotparameter[[i]] = x;
  }
  vcolor = plotparameter[[1]]; vfill = plotparameter[[2]];
  bcolor = plotparameter[[3]]; bfill = plotparameter[[4]];
  if(vcolor[1] == 1 & length(vcolor) == 1) vcolor = NULL; 
  if(vfill[1] == 1 & length(vfill) == 1) vfill = NULL;
  if(bcolor[1] == 1 & length(bcolor) == 1) bcolor = NULL;
  if(bfill[1] == 1 & length(bfill) == 1) bfill = NULL;
  p <- ggplot(data_m, aes(x=variable,y=value)) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(title= title, x= NULL,y = y) +
    theme(axis.text.x=element_text(hjust=0.5,vjust=0.5),
          plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 
  if(violin & scaleviolin) 
    p <- p + geom_violin(aes(color = vcolor, fill = vfill), scale = "width");
  if(violin & !scaleviolin) 
    p <- p + geom_violin(aes(color = vcolor, fill = vfill));
  if(boxplot & !outlier) 
    p <- p + geom_boxplot(aes(color = bcolor, fill = bfill),width=boxw, outlier.size = -1) 
  if (boxplot & outlier)
    p <- p + geom_boxplot(aes(color = bcolor, fill = bfill),width=boxw)  
  if (point)
    p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
                          dotsize = dotsize, binwidth = binwidth)
  if(meanpoint) 
    p <- p + geom_point(stat = "summary",
                        fun = "mean",shape = 1) #fun.y 变成了fun
  if(!missing(ylim)) p <- p + ylim(ylim);
  if(prism) {
    p<- p + theme_prism(base_size = 11) +
      scale_shape_prism() 
    if(!isFALSE(prismpalette)) 
      p <- p + scale_colour_prism(palette = prismpalette) + 
        scale_fill_prism(palette = prismpalette);
  }
  p
}
pd <- plotdfunc(list(ad=ad_c, mci=mci_c,ctl=ctl_c))
VBplot_1way(pd,2,violin = FALSE, outlier = FALSE) + ylim(0,0.2)

#Fig.5f chisq-test
kafang <- function(total,sig){
  tab <- as.table(rbind(total-sig, sig))  
  dimnames(tab) <- list(c("unsig", "sig"), compare)
  tab
  odds <- (tab[2,1]/tab[1,1])/(tab[2,2]/tab[1,2]);
  if(min(apply(tab,2,sum))*min(apply(tab,1,sum))/sum(tab) >= 5) 
    chisq <- chisq.test(tab, correct = FALSE)
  if(min(apply(tab,2,sum))*min(apply(tab,1,sum))/sum(tab) < 5)
    chisq <- chisq.test(tab, correct = TRUE)
  if(min(apply(tab,2,sum))*min(apply(tab,1,sum))/sum(tab) < 1) 
    chisq <- fisher.test(tab)
  p <- chisq$p.value
  c(odds,p,chisq$statistic)
}
total <- c(169,169); compare <- c("mci", "ctl")
kafang(total,c(6,18))
kafang(total,c(14,18))
library(gmodels)
M <- matrix(c(18,6,14,151,163,155),nrow=3,byrow = F,
            dimnames = list(group = c("ctl","mci","ad"),
                            effect = c("sig","unsig")))
chisq.test(M)



##Supplementary Fig.5a ############
#write graphml file and then used Cytoscape to modifiy
load("PCx_DAPnet.RData")
load("moduleinf.RData")

library(igraph)
netgene <- netDAP_4M$netgene;
g <- netDAP_4M$PMFG;
hub = netDAP_4M$hub

gene <- names(V(g));
gn = moduleinf$GNmouse[match(names(V(g)),moduleinf$ori.ID)]
while(sum(duplicated(gn))>0)
  gn[duplicated(gn)] = paste0(gn[duplicated(gn)]," ") 
V(g)$name = gn
degreenet = degree(g);
degreenet = degreenet[degreenet != 1];
degreenet = sort(degreenet,decreasing = T)
link = as_data_frame(g)
write.csv(link, file = "supplementary table3.csv")
genecolor = moduleinf$celltype[match(gene, moduleinf$ori.ID)];
genecolor[genecolor== ""] = "gray"
genecolor[genecolor== "Astro"] = "cyan"
genecolor[genecolor== "Micro"] = "green"
genecolor[genecolor== "Oligo"] = "yellow"
genecolor[genecolor== "Neuro"] = "red"
V(g)$color <- genecolor;
V(g)$shape = "circle"
V(g)$shape[gene %in% hub] <- "square"
V(g)$size = 5;   V(g)$label = NA;
write_graph(g,"SI Fig.5a.graphml","graphml")

## Supplementary Fig.5b-c ##########
rm(list = ls())
load("moduleinf.RData")
load("PCx_DAPnet.RData")
g <- netDAP_4M$PMFG;
link = as_data_frame(g)
gene <- names(V(g));
gn = moduleinf$GNmouse[match(names(V(g)),moduleinf$ori.ID)]
genecolor = moduleinf$celltype[match(gene, moduleinf$ori.ID)];
x <- data.frame(gene,gn,genecolor)
#Supplementary Fig.5b
{
x1 = subset(x,genecolor == "Astro")$gene
x2 = subset(x,genecolor == "Micro")$gene
linknum = sum(link$from %in% x1 & link$to %in% x2) + sum(link$to %in% x1 & link$from %in% x2)

xcell = subset(x,genecolor != "")$gene
xcell = x$gene
rlinkN = NULL
set.seed(1000)
for(i in 1:1000) {
  x1 = sample(xcell,14)
  x2 = sample(xcell[!xcell %in% x1],15)
  lN1 = sum(link$from %in% x1 & link$to %in% x2) + sum(link$to %in% x1 & link$from %in% x2)
  rlinkN = c(rlinkN,lN1)
}
library(ggplot2)
library(ggprism)
plotdata <- data.frame(rlinkN)
ggplot(plotdata,aes(x = rlinkN))+ theme_bw()+ 
  geom_density(fill = "gray",color = "black",
               adjust = 1.5) +
  theme_prism()

sum(rlinkN>linknum)
linknum
}
#Supplementary Fig.5c
{
source("stringfunction.R")
stringdb <- STRINGdb$new( version="11", species=10090,
                          score_threshold=200, input_directory="")
x$new.ID = moduleinf$new.ID[match(names(V(g)),moduleinf$ori.ID)]
allmap <- STRmap(x, colname = "gn")
link_STR <- STRlink(subset(x,genecolor %in% c("Astro","Micro"))$gn, 
                    allmap, input_type = "gn")
x1 = subset(x,genecolor == "Astro")$gn
x2 = subset(x,genecolor == "Micro")$gn
STRnum = sum(link_STR$from %in% x1 & link_STR$to %in% x2) + sum(link_STR$to %in% x1 & link_STR$from %in% x2)

rd = randSTRlink(fixedID = NULL, allmap, input_type = "gn",
                 randgenenum = 29)
rSTR_N = NULL
for(i in 1:1000) {
  STRnum1 = sum(rd[[i]]$from %in% x1 & rd[[i]]$to %in% x2) + sum(rd[[i]]$to %in% x1 & rd[[i]]$from %in% x2)
  rSTR_N = c(rSTR_N,STRnum1)
}
library(ggplot2)
library(ggprism)
plotdata <- data.frame(rSTR_N)
ggplot(plotdata,aes(x = rSTR_N))+ theme_bw()+ 
  geom_density(fill = "gray",color = "black",
               adjust = 3) +
  xlim(c(0,10)) +
  theme_prism()
sum(rSTR_N>STRnum)
STRnum 
}
## Supplementary Fig.6 a-c ####
rm(list = ls())
load("AD_LFQ_impute_DAP4M.RData")
rownames(quant) <- info$ori.ID
library(Hmisc)
library(ggplot2)
library(ggprism)
gn <- unique(info$GN[duplicated(info$GN)])
ast <- info$ori.ID[info$GN %in% gene$GN[gene$celltype =="Astro"]]
mic <- info$ori.ID[info$GN %in% gene$GN[gene$celltype =="Micro"]]

sig_cor<-function(cor_matrix,p=0.05,r=NA,type=1){
  cor_matrix$P[upper.tri(cor_matrix$P,diag = TRUE)]<-NA;
  cor_matrix$r[upper.tri(cor_matrix$P,diag = TRUE)]<-NA;
  if(!is.na(p)){
    logi<-cor_matrix$P>p;
    cor_matrix$P[logi]<-NA;
    cor_matrix$r[logi]<-NA;
  }
  if(!is.na(r)){
    logi2<-abs(cor_matrix$r)<r
    cor_matrix$P[logi2]<-NA;
    cor_matrix$r[logi2]<-NA;
  }
  if(!is.na(p)|!is.na(r)){
    posit<-which(!is.na(cor_matrix$P), arr.ind=TRUE);
    if(type==2){
      #对这些位置排序，提取出方阵
      column<-sort(unique(c(posit[,1],posit[,2])));
      cor_matrix$P<-cor_matrix$P[column,column];
      cor_matrix$r<-cor_matrix$r[column,column];
    }
    else if(type==1){
      #对这些位置排序提取矩阵
      roworder<-sort(unique(posit[,1]));
      column<-sort(unique(posit[,2]));
      cor_matrix$P<-cor_matrix$P[roworder,column];
      cor_matrix$r<-cor_matrix$r[roworder,column];
    }
  }
  cor_matrix
}
orderCor <- function(cormat, pmat) {
  postion<-which(!is.na(cormat),arr.ind = TRUE);
  pos<-which(!is.na(cormat));
  ro<-attr(cormat,"dimnames")[[1]][postion[,1]];
  co<-attr(cormat,"dimnames")[[2]][postion[,2]];
  r<-cormat[pos];
  p<-pmat[pos];
  flat<-data.frame(row=ro,column=co,cor=r,p=p);
  orderCor<-flat[order(abs(flat$cor),decreasing = TRUE),];
  rownames(orderCor)<-order(abs(orderCor$cor),decreasing = TRUE);
  orderCor
}
corlist <- function(group = "diagnose", name, times = 10,quant) {
  n <- min(table(metainfo[,group])[table(metainfo[,group])>9])
  n =50
  quant2 <- t(quant[,metainfo[,group] == name])
  if(nrow(quant2) > n) {
    for(j in 1:times) {
      set.seed(j)
      pickSample <- sample(1:nrow(quant2),n)
      quant3 <- quant2[pickSample,]
      x1 <- rcorr(quant3,type = "pearson")
      x1 <-sig_cor(x1,p=1,r=NA,type = 2)
      x1<-orderCor(x1$r,x1$P)
      for( i in which(x1$p == 0)){
        fit <- lm(quant3[,x1$row[i]]~quant3[,x1$column[i]]);
        tstats <- coef(fit)/sqrt(diag(vcov(fit)));
        x1$p[i] <- 2*pt(abs(tstats),df=df.residual(fit),lower.tail = FALSE)[2]
      }
      if(j == 1) {x <- x1; r1 <- x1$cor} else r1 <- r1+x1$cor
    }
    x$cor <- r1/times
    x$p <- 2 * (1 - pt(abs(x$cor) * sqrt(80-2)/sqrt(pmax(1 - x$cor * x$cor, 0)), 80-2))
  } else {
    x <- rcorr(quant2,type = "pearson")
    x <-sig_cor(x,p=1,r=NA,type = 2)
    x<-orderCor(x$r,x$P)
    for( i in which(x$p == 0)){
      fit <- lm(quant2[,x$row[i]]~quant2[,x$column[i]]);
      tstats <- coef(fit)/sqrt(diag(vcov(fit)));
      x$p[i] <- 2*pt(abs(tstats),df=df.residual(fit),lower.tail = FALSE)[2]
    }
  }
  x$p.adj <- p.adjust(x$p,method = "fdr") 
  x
}
linkfilter <- function(x){
  pos1 = x[,1] %in% ast & x[,2] %in% mic
  pos2 = x[,1] %in% mic & x[,2] %in% ast
  x[pos1|pos2,]
}
ctl <- corlist(name = "CTL",times=10,quant=quant)
mci <- corlist(name = "MCI",times=10,quant=quant)
ad <- corlist(name = "AD",times=10,quant=quant)
ctl_c <- linkfilter(ctl);
mci_c <- linkfilter(mci);
ad_c <- linkfilter(ad);
#Supplementar Fig.6c 
#write graphml file and then used Cytoscape to modifiy
{
  mergelink <- function(x) {
    net <- x[[1]]
    net$cor[net$p > 0.05] <- NA
    colnames(net)[3] <- names(x)[1];
    net <- net[,-c(4:ncol(net))]
    POS <- is.na(net[,3])
    for(i in 2:length(x)) {
      net$x <- NA; colnames(net)[i+2] <- names(x)[i]
      pos = x[[i]]$p < 0.05
      pos1 = paste0(x[[i]][pos,1],x[[i]][pos,2])
      net[match(pos1, paste0(net[,1],net[,2])),i+2] <- x[[i]][pos,3]
      POS <- c(POS&is.na(net[,i+2]))
    }
    net <- net[!POS,]
    net[,1] <- info$GN[match(net[,1],info$ori)]
    net[,2] <- info$GN[match(net[,2],info$ori)]
    net
  }
  cellinteraction <- mergelink(list(ad=ad_c,mci=mci_c,ctl = ctl_c))
  library(igraph)
  link <- cellinteraction
  g <- graph.data.frame(link)
  netgene <- names(V(g)); 
  genecolor <- gene$celltype[match(netgene,gene$GN)]
  genecolor[genecolor== ""] = "gray"
  genecolor[genecolor== "Astro"] = "cyan"
  genecolor[genecolor== "Micro"] = "green"
  genecolor[genecolor== "Oligo"] = "yellow"
  genecolor[genecolor== "Neuro"] = "red"
  V(g)$color <- genecolor
  linkcolor <- rep("black",nrow(link))
  linkcolor[is.na(link$mci)&is.na(link$ctl)] <- "red"
  linkcolor[is.na(link$ad)&is.na(link$ctl)] <- "#FF4500"
  linkcolor[is.na(link$ad)&is.na(link$mci)] <- "yellow"
  linkcolor[(!is.na(link$ad))&(!is.na(link$mci))] <- "gray11"
  linkcolor[(!is.na(link$mci))&(!is.na(link$ctl))] <- "gray22"
  linkcolor[(!is.na(link$ad))&(!is.na(link$ctl))] <- "darkblue"
  table(linkcolor)
  E(g)$color <- linkcolor
  E(g)$lwd <- apply(link[,3:5],1,mean,na.rm=TRUE)
  #write graphml file and then used Cytoscape to modifiy
  #write_graph(g,"SI Fig.6c.graphml","graphml")
}

#Supplementar Fig.6a
plotdfunc <- function(x,name = "cor"){
  plotdata <- NULL
  for( i in 1:length(x)) {
    plotdata <- rbind(plotdata,data.frame(group = names(x)[i],value = abs(x[[i]][,name])))
  }
  plotdata$group <- factor(plotdata$group,levels = c("ctl","mci","ad"));
  p1 <- t.test(value~group,plotdata,group %in% c("ctl","mci"))$p.value
  p2 <- t.test(value~group,plotdata,group %in% c("ctl","ad"))$p.value
  p3 <- t.test(value~group,plotdata,group %in% c("mci","ad"))$p.value
  print(paste0("t.test:Pctlvsmci: ",p1,"; Pctlvsad: ",p2,"; Pmcivsad: ",p3))
  p1 <- wilcox.test(value~group,plotdata,group %in% c("ctl","mci"))$p.value
  p2 <- wilcox.test(value~group,plotdata,group %in% c("ctl","ad"))$p.value
  p3 <- wilcox.test(value~group,plotdata,group %in% c("mci","ad"))$p.value
  print(paste0("wilcox:Pctlvsmci: ",p1,"; Pctlvsad: ",p2,"; Pmcivsad: ",p3))
  plotdata
}
VBplot_1way <- function(data_m, coln ,violin = FALSE,boxplot = TRUE,
                        outlier = TRUE,point = FALSE,
                        meanpoint = TRUE,prism = TRUE, prismpalette = "floral",
                        vcolor = TRUE, vfill = FALSE, scaleviolin = FALSE,
                        bcolor = TRUE, bfill = FALSE, 
                        title = NULL, y = NULL, boxw = 0.2, 
                        dotsize = 1, binwidth = 1.5,
                        ylim, ...){
  variable = factor(data_m[,1]);
  data_m <- data.frame(variable = variable, value = data_m[,coln]);
  if(length(vcolor) > 1 | length(bcolor) > 1) prismpalette = FALSE
  plotparameter = list(vcolor=vcolor,vfill=vfill,bcolor=bcolor,bfill=bfill)
  for(i in seq_along(plotparameter)) {
    x = plotparameter[[i]]
    if(isTRUE(x)) 
      x = variable else if(length(x) == nlevels(variable))
        x = as.character(factor(variable,labels = x)) else if(length(x) == length(variable))
          x = x else x = 1;
          plotparameter[[i]] = x;
  }
  vcolor = plotparameter[[1]]; vfill = plotparameter[[2]];
  bcolor = plotparameter[[3]]; bfill = plotparameter[[4]];
  if(vcolor[1] == 1 & length(vcolor) == 1) vcolor = NULL; 
  if(vfill[1] == 1 & length(vfill) == 1) vfill = NULL;
  if(bcolor[1] == 1 & length(bcolor) == 1) bcolor = NULL;
  if(bfill[1] == 1 & length(bfill) == 1) bfill = NULL;
  p <- ggplot(data_m, aes(x=variable,y=value)) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(title= title, x= NULL,y = y) +
    theme(axis.text.x=element_text(hjust=0.5,vjust=0.5),
          plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 
  if(violin & scaleviolin) 
    p <- p + geom_violin(aes(color = vcolor, fill = vfill), scale = "width");
  if(violin & !scaleviolin) 
    p <- p + geom_violin(aes(color = vcolor, fill = vfill));
  if(boxplot & !outlier) 
    p <- p + geom_boxplot(aes(color = bcolor, fill = bfill),width=boxw, outlier.size = -1) 
  if (boxplot & outlier)
    p <- p + geom_boxplot(aes(color = bcolor, fill = bfill),width=boxw)  
  if (point)
    p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
                          dotsize = dotsize, binwidth = binwidth)
  if(meanpoint) 
    p <- p + geom_point(stat = "summary",
                        fun = "mean",shape = 1) #fun.y 变成了fun
  if(!missing(ylim)) p <- p + ylim(ylim);
  if(prism) {
    p<- p + theme_prism(base_size = 11) +
      scale_shape_prism() 
    if(!isFALSE(prismpalette)) 
      p <- p + scale_colour_prism(palette = prismpalette) + 
        scale_fill_prism(palette = prismpalette);
  }
  p
}
pd <- plotdfunc(list(ad=ad_c, mci=mci_c,ctl=ctl_c))
VBplot_1way(pd,2,violin = FALSE, outlier = FALSE) + ylim(0,0.2)

#Supplementary Fig.6b chisq-test
sum(ad_c$p<0.05); sum(mci_c$p<0.05); sum(ctl_c$p<0.05)
total <- c(143,143); compare <- c("mci", "ctl")
kafang(total,c(7,4))
kafang(total,c(7,7))
library(gmodels)
M <- matrix(c(7,4,7,136,139,136),nrow=3,byrow = F,
            dimnames = list(group = c("ctl","mci","ad"),
                            effect = c("sig","unsig")))
chisq.test(M)


