#1.DEPs-associated protein(DAP) acquisition
{
library(DDPNA)
source("DAP_func.R")
load("5FAD_PCx_DEP_ModEnrich.RData")
load("5FAD_PCx_WGCNAnet.Rdata")
load("5FAD_PCx2_adj_final.RData")
stat <- read.csv("5FAD_PCx_stat.csv")

#DAP_3M_13 indicates DAP of 3-month-old DEPs in module 13 (Figure 3.c), naming below is all follow rules like this.
DAP_3M_13 <- DAPbyMod(quant_final, moduleinf, 13, mo3$DEPinf$ori.ID,
                      IDname="DEP_3M", cor.adj="fdr", coln = "ori.ID");
DAP_3M_2 <- DAPbyMod(quant_final, moduleinf, 2, mo3$DEPinf$ori.ID,
                     IDname="DAP_3M_extra", cor.adj="fdr", coln = "ori.ID");
DAP_3M_11 <- DAPbyMod(quant_final, moduleinf, 11, mo3$DEPinf$ori.ID,
                      IDname="DAP_3M_extra", cor.adj="fdr", coln = "ori.ID");

DAP_4M_3 <- DAPbyMod(quant_final, moduleinf, 3, mo4$DEPinf$ori.ID,
                     IDname="DEP_4M", cor.adj="fdr", coln = "ori.ID");
DAP_4M_10 <- DAPbyMod(quant_final, moduleinf, 10, mo4$DEPinf$ori.ID,
                      IDname="DAP_4M_extra", cor.adj="fdr", coln = "ori.ID");

DAP_6M_3 <- DAPbyMod(quant_final, moduleinf, 3, mo6$DEPinf$ori.ID,
                     IDname="DEP_6M",cor.adj="fdr", coln = "ori.ID");
DAP_6M_4 <- DAPbyMod(quant_final, moduleinf, 4, mo6$DEPinf$ori.ID,
                     IDname="DAP_6M_extra", cor.adj="fdr", coln = "ori.ID");
DAP_6M_2 <- DAPbyMod(quant_final, moduleinf, 2, mo6$DEPinf$ori.ID,
                     IDname="DEP_6M",cor.adj="fdr", coln = "ori.ID");
DAP_6M_13 <- DAPbyMod(quant_final, moduleinf, 13, mo6$DEPinf$ori.ID,
                      IDname="DEP_6M",cor.adj="fdr", coln = "ori.ID");

DAP_9M_3 <- DAPbyMod(quant_final, moduleinf, 3, mo9$DEPinf$ori.ID,
                   IDname="DEP_9M",cor.adj="fdr", coln = "ori.ID");


DAP_11M_2 <- DAPbyMod(quant_final, moduleinf, 2, mo11$DEPinf$ori.ID,
                      IDname="DEP_11M",cor.adj="fdr", coln = "ori.ID");
DAP_11M_3 <- DAPbyMod(quant_final, moduleinf, 3, mo11$DEPinf$ori.ID,
                      IDname="DAP_11M_extra", cor.adj="fdr", coln = "ori.ID");
DAP_11M_7 <- DAPbyMod(quant_final, moduleinf, 7, mo11$DEPinf$ori.ID,
                      IDname="DAP_11M_extra", cor.adj="fdr", coln = "ori.ID");
DAP_11M_13 <- DAPbyMod(quant_final, moduleinf, 13, mo11$DEPinf$ori.ID,
                       IDname="DAP_11M_extra", cor.adj="fdr", coln = "ori.ID");

#DAPcore is the combination of DAPs from 4, 6, 9, 11-month-old DEP only in module3
DAPcore <- intersect(DAP_4M_3,DAP_6M_3); 
DAPcore <- intersect(DAPcore, DAP_9M_3); 
DAPcore <- intersect(DAPcore, DAP_11M_3);
DAPcoreGN <- unique(inf$GN[inf$ori.ID %in% DAPcore])
}

#2.Customiszed different month-old DEPs overlapping plot, related to Figure 2.c
{
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
}

#3.Customized cell type enrichment plot, related to Figure
{
library(reshape2)
library(ggplot2)
library("cowplot")
data <- read.csv("cell type enrichment.csv",header = T, stringsAsFactors = F)
rownames(data) <- data$X
data <- data[,-1]
data2 <- -log10(data)
dep3 <- data2[1,]
dep4 <- data2[2,]
dep6 <- data2[3,]
dep9 <- data2[4,]
dep11 <- data2[5,]
dap4 <- data2[7,]
dap6 <- data2[8,]  
dap9 <- data2[9,]
dap11 <- data2[10,]
dap3_13 <- data2[6,]
dap11_13 <- data2[11,]
dep3 <- melt(dep3) 
dep4 <- melt(dep4) 
dep6 <- melt(dep6) 
dep9 <- melt(dep9) 
dep11 <- melt(dep11) 
dap3_13 <- melt(dap3_13) 
dap4 <- melt(dap4) 
dap6 <- melt(dap6) 
dap9 <- melt(dap9) 
dap11 <- melt(dap11) 
dap11_13 <- melt(dap11_13) 

pdep3 <- ggplot(data = dep3, aes(x=variable,y=value,fill=variable))
dep.3 <- pdep3 +geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(plot.title = element_text(hjust = 0.5,size = 15))+
  guides(fill="none")

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

pdap3_13 <- ggplot(data = dap3_13, aes(x=variable,y=value,fill=variable))
dap.3.13 <- pdap3_13 +geom_bar(stat = "summary",fun = mean,width = 1)+
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

pdap11_13 <- ggplot(data = dap11_13, aes(x=variable,y=value,fill=variable))
dap.11_13 <- pdap11_13+geom_bar(stat = "summary",fun = mean,width = 1)+
  ylim(0,6)+
  geom_hline(yintercept = -log10(0.05),color = "gray")+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x = NULL,axis.text.y = element_text(size=15))+
  guides(fill="none")

plot_grid(dep.4,dap.4,dep.6,dap.6,dep.9,dap.9,dep.11,dap.11,nrow = 4,ncol = 2)
}