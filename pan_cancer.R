
rm(list = ls())
setwd("D:/000Researching/020Yao's Research/pan-cancer/Cuprotosis/Data")
dir <- "D:/000Researching/020Yao's Research/pan-cancer/Cuprotosis/"
getwd()

# Load library ------------------------------------------------------------

# install.packages("xlsx")
# install.packages("rJava")
# BiocManager::install("limma")
# BiocManager::install("GenVisR")
# install.packages("forestploter")
# BiocManager::install("impute")
library(GenVisR) 
set.seed(383)
library(BiocManager)
library(rJava)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(xlsx)
library(ggsignif)
library(limma)
library(dplyr)
library(plyr)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(scales)
library(corrplot)
library(ggcorrplot)
library(Hmisc)
library(ggradar)
library(ggrepel)
library(ggsci)
library(forestploter)
library(RColorBrewer)

# Expression in normal tissue ---------------------------------------------
gene <- c("FDX1", "LIAS", "LIPT1", "DLAT", "DLD", "PDHA1", "PDHB", "SLC31A1")
Fdat <- read.table("../../HPA/rna_tissue_fantom.tsv", header = T, sep = "\t")
fdat <- Fdat[Fdat$Gene.name %in% gene,]
writexl::write_xlsx(fdat, "../../HPA/gene_fantom.xlsx")
unique(Fdat$Tissue)
p <- ggplot(fdat, aes(x=Tissue, y = Scaled.tags.per.million, fill= Gene.name))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_lancet(alpha = 0.5)+
  facet_grid(Gene.name~., scales = "free")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  ylab("Expression(Scaled tags per million)"); p
ggsave(filename =paste0(dir, "Graph/","Fantom","_normal.svg"), plot = p,width = 11, height = 5.5)
Hdat <- read.table("../../HPA/rna_tissue_hpa.tsv", header = T, sep = "\t")
hdat <- Hdat[Hdat$Gene.name %in% gene, ]
unique(hdat$Tissue) 
tis <- unique(hdat[hdat$Tissue %in% unique(Fdat$Tissue), ]$Tissue)
tis <- c(tis, "bone marrow")
hdat <- hdat[hdat$Tissue %in% tis,]
writexl::write_xlsx(hdat, "../../HPA/gene_hpa.xlsx")
p <- ggplot(hdat, aes(x=Tissue, y = TPM, fill= Gene.name))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_lancet(alpha = 0.5)+
  facet_grid(Gene.name~., scales = "free")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  ylab("Expression(mean of transcription per million)"); p
ggsave(filename =paste0(dir, "Graph/","HPA","_normal.svg"), plot = p,width = 10, height = 5.5)
Gdat <- read.table("../../HPA/rna_tissue_gtex.tsv", header = T, sep = "\t")
gdat <- Gdat[Gdat$Gene.name %in% gene, ]
unique(gdat$Tissue)
writexl::write_xlsx(gdat, "../../HPA/gene_GTEx.xlsx")
p <- ggplot(gdat, aes(x=Tissue, y = nTPM, fill= Gene.name))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_lancet(alpha = 0.5)+
  facet_grid(Gene.name~., scales = "free")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  ylab("Expression(transcription per million)"); p
ggsave(filename =paste0(dir, "Graph/","GTEx","_normal.svg"), plot = p,width = 10, height = 5.5)


diff <- readxl::read_xlsx("diff.xlsx", sheet = 1)
diff <- diff[,c(4,1,2,3)]
for (i in 2:8) {
  diff0 <- readxl::read_xlsx("diff.xlsx", sheet = i)
  diff0 <- diff0[,c(4,1,2,3)]
  diff <- merge(diff, diff0, all =T)
}
writexl::write_xlsx(diff, "Results.xlsx")

# diff <- readxl::read_xlsx("Results.xlsx", sheet=1)
diff$Sample <- as.factor(diff$Sample)
diff1 <- group_by(diff, Sample)
col <-c("#F0AD4E","#D9534F","#337AB7","#5CB85C")
name <- colnames(diff1)[3:10]
for (i in 4:11) {
  dat <- diff1[,c(2,3,i)]
  colnames(dat) <- c("Sample", "Group", "Gene")
  p <- ggplot(dat, aes(x=Sample, y=Gene, color = Group))+
    geom_boxplot()+theme_classic()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    ylab(paste0(colnames(diff1)[i]," Expression"))+
    scale_color_manual(values = col)+
    stat_compare_means(aes(group = Group),method="wilcox.test",hide.ns = T,
                       label="p.signif",vjust = 0.8)
  ggsave(filename =paste0(dir, "Graph/",colnames(diff1)[i],"_diff.svg"), plot = p,width = 13, height = 4) 
}
diff_p <- data.frame(1:34)
row.names(diff_p) <- t$Sample
for (i in 4:11) {
  form <- as.formula(paste0(colnames(diff1)[i], "~Group"))
  t <- compare_means(form, data = diff1, group.by = "Sample",method = "wilcox.test", 
                     paired = F) 
  diff_p <- cbind(diff_p, t$p.signif)
  colnames(diff_p)[i-2] <- colnames(diff1)[i]
}
diff_p <- diff_p[,-1]
diff_p <- diff_p[order(row.names(diff_p)), ]

diff1 <- diff1[,-1]
diff2 <- aggregate(.~Sample+Group, median, data = diff1)
diff2 <- diff2[order(diff2$Sample),]
write.xlsx(diff2, "diff_T_N.xlsx", sheetName = "Median", append = T)
diff2 <- readxl::read_xlsx("diff_T_N.xlsx", sheet = 4)
diff2 <- as.data.frame(t(diff2))
colnames(diff2) <- diff2[1,]
diff2 <- diff2[-1,]
diff2[,c(1:34)] <- lapply(diff2[,c(1:34)], as.numeric)
diff2 <- as.matrix(diff2)
for (i in seq(2,68,2)) {
  a <- diff2[i, ]-diff2[i-1, ]
  
  
}
diff1$Sample <- as.factor(diff1$Sample)
diff1$Group <- as.factor(diff$Group)
diff2 <- diff1 %>% group_by(Sample, Group) %>% summarise(FDX1=mean(FDX1), LIAS=mean(LIAS), LIPT1=mean(LIPT1), DLD=mean(DLD),
                                                    DLAT=mean(DLAT), PDHA1=mean(PDHA1), PDHB = mean(PDHB), SLC31A1=mean(SLC31A1))
sum <- diff2  


p <- ggplot(diff1, aes(x=Sample, y=FDX1, color = Group))+
    geom_boxplot()+theme_classic()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    ylab(paste0("FDX1"," Expression"))+
    scale_color_manual(values = col)+
    stat_compare_means(aes(group = Group),method="wilcox.test",hide.ns = T,
                       label="p.signif",vjust = 0.8)


diff_T <- diff[diff$Group == "Tumor", ]
diff_T <- diff_T[, c(2,4:11)]
diff_T <- gather(diff_T, key=Gene, value = Expression, -Sample)
mean <- aggregate(diff_T$Expression,by = list(Sample=diff_T$Sample, Gene=diff_T$Gene),FUN = mean, data = diff_T)
sd <- aggregate(diff_T$Expression,by = list(Sample=diff_T$Sample, Gene=diff_T$Gene),FUN = sd)
diff_T <- data.frame(mean, sd$x)
colnames(diff_T) <- c("Sample", "Gene", "Expression", "sd")
diff_T <- separate(diff_T, Sample, into = c("Sample", "dis"), sep="[(]")
p <- ggplot(diff_T, aes(x=Sample, y=Expression, fill = Gene))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = Expression-sd,ymax = Expression+sd), position = position_dodge(), width=.2)+
  facet_grid(Gene~., scales = "free")+
  scale_fill_lancet(alpha = 0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab( "Expression(FPKM)");p

ggsave(filename =paste0(dir, "Graph/","GTEx","_Tumor.svg"), plot = p,width = 10, height = 5)

exp <- readxl::read_xlsx("exp.xlsx", sheet = 1)
exp <- exp[,c(4,1,2,3)]
for (i in 2:8) {
  exp0 <- readxl::read_xlsx("exp.xlsx", sheet = i)
  exp0 <- exp0[,c(4,1,2,3)]
  exp <- merge(exp, exp0, all =T)
}
write.xlsx(exp, "exp_merge.xlsx")

res <- readxl::read_xlsx("exp_merge.xlsx")[,-1]

sam <- unique(res$Sample)
diff_exp <- data.frame()
logFC <- data.frame(1:8)
row.names(logFC) <- row.names(out)
avgexp <- data.frame(1:8)
row.names(avgexp) <- row.names(out)
for (i in 1:length(sam)) {
  res2 <- filter(res, Sample ==sam[i])
  res3 <- res2[,-c(2,3)]
  res3 <- as.data.frame(t(res3))
  colnames(res3) <- res3[1,]
  res3 <- res3[-1,]
  res3[,c(1:ncol(res3)) ] <- lapply(res1[,c(1:ncol(res3)) ], as.numeric)
  list <- factor(res2$Group)
  head(list)
  list <- model.matrix(~factor(list)+0)
  colnames(list) <- substr(colnames(list), 13, nchar(colnames(list)))
  fit <- lmFit(res3, list)
  contr <- makeContrasts(Normal-Tumor, levels = list)
  fit1 <- contrasts.fit(fit, contr)
  fit1 <- eBayes(fit)
  out <- topTable(fit1, n=Inf, coef = 1)
  logFC <- cbind(logFC, out$logFC)
  colnames(logFC)[i+1] <- sam[i]
  avgexp <- cbind(avgexp, out$AveExpr)
  colnames(avgexp)[i+1] <- sam[i]
  out$Sample <- rep(sam[i], 8)
  diff_exp <- rbind(diff_exp,out)
}

diff_exp$Gene <- row.names(diff_exp)
write.xlsx(diff_exp, "diff_exp.xlsx", sheetName = "diff_exp")
write.xlsx(logFC, "diff_exp.xlsx", sheetName = "logFC",append = T)
write.xlsx(avgexp, "diff_exp.xlsx", sheetName = "avgExp",append = T)

logFC <- logFC[,-1]
logMA <- as.matrix(logFC)
avgExp <- avgexp[,-1]
avgExp <- as.matrix(avgExp)

diff_p <- diff_p[, -c(1,10)]
diff_p1 <- as.data.frame(t(diff_p))
diff_p1 <- as.matrix(diff_p1)

diff3 <- diff1[,-2]
diff3 <- as.data.frame(t(diff3))
colnames(diff3) <- diff3[1,]
diff3 <- diff3[-1,]
diff3[,c(1:22046)] <- lapply(diff3[,c(1:22046)], as.numeric)
diff3 <- as.matrix(diff3)
bk <- c(seq(-1,-0.01,0.01), seq(0,1,0.01))
color <- c(colorRampPalette(colors = c("navy", "white"))(length(bk)/2),
           colorRampPalette(colors = c("white", "firebrick3"))(length(bk)/2))
color <- c(colorRampPalette(colors = c("#FDE725FF", "#1E9D89FF"))(length(bk)/2),
           colorRampPalette(colors = c("#1E9D89FF", "#3C508BFF"))(length(bk)/2))
pdf(paste0(dir, "/Graph/diff_heatplot.pdf"), width = 12, height = 5)
pheatmap(diff2,  cutree_cols = 10,legend_breaks = seq(-1,1,0.5),breaks = bk,col = color)
dev.off()

Heatmap(diff2, cell_fun = function(j, i, x, y, w, h, col) {
  grid.text(diff_p[i, j], x, y) })


a <- viridis_pal()(50)

waterfall(brcaMAF, mainRecurCutoff = 0.06)

readx <- function(x){
  dat <- readxl::read_xlsx(paste0(x,".xlsx"), sheet = 1)
  dat <- dat[,c(4,1,2,3)]
  for (i in 2:8) {
    dat0 <- readxl::read_xlsx(paste0(x,".xlsx"), sheet = i)
    dat0 <- dat0[,c(4,1,2,3)]
    dat <- merge(dat, dat0, all =T)
  }
  return(dat)
}
tmb <- readx("TMB")
writexl::write_xlsx(tmb, "TMB_merge.xlsx")
cordat <- corda(tmb)
A <- unique(cordat$Sample)
table <- data.frame()
corMerge <- function(x){
  for (i in 1:length(A)) {
    B <- A[i]
    dat0 <- cordat[which(cordat$Sample %in% B), ]
    cor0 <- rcorr(as.matrix(dat0[, c("TMB", x)]), type = "spearman")
    table0 <- data.frame(B, x,cor0$r[2], cor0$P[2])
    table <- rbind(table, table0)
  }
  colnames(table) <- c("Sample","Gene","R", "pvalue")
  return(table)
}
TMB_merge_cor <- lapply(colnames(cordat)[3:10], corMerge)
names(TMB_merge_cor) <- colnames(cordat)[3:10]
merge_cor <- as.data.frame(TMB_merge_cor)
write.xlsx(merge_cor, "TMB_merge_cor.xlsx")

redar <- function(x){
  table <- data.frame()
  for (i in 1:length(A)) {
    B <- A[i]
    dat0 <- cordat[which(cordat$Sample %in% B), ]
    cor0 <- rcorr(as.matrix(dat0[, c("TMB", x)]), type = "spearman")
    table0 <- data.frame(B, cor0$r[2], cor0$P[2])
    table <- rbind(table, table0)
  }
  colnames(table) <- c("Sample", paste0("R",x), "pvalue")
  table$pvalue <- paste("p=", round(table$pvalue, 4))
  data <- unite(table,"Sample", Sample, pvalue,  sep = ",")
  data <- as.data.frame(t(data))
  colnames(data) <- data[1,]
  data <- data[-1, ]
  data[,1:ncol(data)] <- as.numeric(data[,1:ncol(data)])
  data$group <- row.names(data)
  row.names(data) <- NULL
  data <- data[, c(ncol(data), 1:ncol(data)-1)]
  p <- ggradar(data,grid.min = -0.35, grid.mid = 0, grid.max = 0.5 ,values.radar = c("-0.35", "0", "0.5"),
               base.size = 14, axis.label.size = 4, grid.label.size = 6, axis.label.offset = 1.07,
               background.circle.transparency = 0.2,group.point.size = 3, group.colours = "#337AB7", plot.title = x)+
    theme(plot.title = element_text(hjust = 0.5,size = 16))
  
  ggsave(p, filename =paste0(dir, "Graph/TMB_redarplot_", x,".svg"), width = 10, height = 10)
}
sapply(colnames(cordat)[3:10], redar)


msi <- readx("MSI")
cordat <- corda(msi)
A <- unique(cordat$Sample)
writexl::write_xlsx(msi, "MSI_merge.xlsx")
table <- data.frame()
corMerge <- function(x){
  for (i in 1:length(A)) {
    B <- A[i]
    dat0 <- cordat[which(cordat$Sample %in% B), ]
    cor0 <- rcorr(as.matrix(dat0[, c("MSI", x)]), type = "spearman")
    table0 <- data.frame(B, x,cor0$r[2], cor0$P[2])
    table <- rbind(table, table0)
  }
  colnames(table) <- c("Sample","Gene","R", "pvalue")
  return(table)
}
MSI_merge_cor <- lapply(colnames(cordat)[3:10], corMerge)
names(MSI_merge_cor) <- colnames(cordat)[3:10]
merge_cor <- as.data.frame(MSI_merge_cor)
write.xlsx(merge_cor, "MSI_merge_cor.xlsx")


NEO <- readx("NEO")  
corda <- function(x){
  cordat <- x[, -1]
  cordat <- separate(cordat, Sample, into = c("Sample", "name"), sep = "[(]")
  cordat <- cordat[,-2]
  cordat <- cordat[-which(cordat$Sample %in% c("SARC","MESO","DLBC","UVM")),]
  return(cordat)
}
cordat <- corda(NEO)
A <- unique(cordat$Sample)
writexl::write_xlsx(NEO, "NEO_merge.xlsx")

table <- data.frame()
corMerge <- function(x){
  for (i in 1:length(A)) {
    B <- A[i]
    dat0 <- cordat[which(cordat$Sample %in% B), ]
    cor0 <- rcorr(as.matrix(dat0[, c("Neoantigen", x)]), type = "spearman")
    table0 <- data.frame(B, x,cor0$r[2], cor0$P[2])
    table <- rbind(table, table0)
  }
  colnames(table) <- c("Sample","Gene","R", "pvalue")
  return(table)
}
NEO_merge_cor <- lapply(colnames(cordat)[3:10], corMerge)
names(NEO_merge_cor) <- colnames(cordat)[3:10]
merge_cor <- as.data.frame(NEO_merge_cor)
write.xlsx(merge_cor, "NEO_merge_cor.xlsx")

dat <- data.frame()
for (i in 1:8) {
  dat <- rbind(dat, NEO_merge_cor[[i]])
}
dat$p.adj <- -log10((dat$pvalue)+0.001)
fac <- c("FDX1","LIAS", "LIPT1","DLD","DLAT","PDHA1","PDHB", "SLC31A1")
dat$Gene <- factor(dat$Gene, levels = fac, labels = fac)
p <- ggplot(dat, aes(Sample,R,fill=Gene,size =p.adj)) +  #划线
  geom_point(shape=21)+ #设置点
  scale_fill_lancet(alpha = 0.7)+
  theme_classic()+
  ylab("Correlation coefficient(spearman)")+
  labs(title = "NEO")+
  theme(
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=10,face="plain",color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1,  vjust = 1),
    legend.title=element_text(size=12,face="plain",color="black"),
    legend.background = element_blank(),legend.key.size =unit(0.5, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )+ guides(fill = guide_legend(override.aes = list(size=4)))
ggsave(filename =paste0(dir, "Graph/Corplot_","NEO",".svg"), plot = p,width = 8, height = 4)
B <- c(1,5,9,13,17,21,25,29)
readx <- function(x){
  dat <- readxl::read_xlsx(paste0(x,".xlsx"), sheet = 1)
  dat <- dat[,c(1,3,4,5,6)]
  dat <- dat[order(dat[,1]),]
  for (i in B) {
    dat0 <- readxl::read_xlsx(paste0(x,".xlsx"), sheet = i)
    dat0 <- dat0[,c(1,3,4,5,6)]
    dat0 <- dat0[order(dat0[,1]),]
    dat <- cbind(dat, dat0)
  }
  return(dat)
}
pro <- readx("prog")
pro <- pro[-which(pro$CancerCode %in% c("TARGET-ALL-R(N=99)","TARGET-LAML(N=142)","TARGET-NB(N=151)",
                                       "TCGA-DLBC(N=44)","TCGA-MESO(N=84)","TCGA-SKCM-M(N=347)",
                                       "TCGA-SKCM-P(N=97)","TCGA-THYM(N=117)","TCGA-UVM(N=74)","TCGA-SARC(N=254)")),]
write.xlsx(pro, "pro_OS_merge.xlsx")
tm <- forestploter::forest_theme(base_size = 6.5, ci_pch = 4, ci_col = "#D9534F", ci_lty = 1, ci_lwd = 2,
                                 ci_Theight = 0.2, refline_lwd = 1, refline_lty = "dashed", refline_col = "#337AB7",
                                 vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "#5CB85C")
# pro[, c(2:4)] <- sapply(pro[, c(2:4)], log2)
pro$PDHB_se <- (log(pro$PHDB_Upper)-log(pro$PHDB_Upper))/1.96
pro$SLC31A1 <- paste0(rep(" ", 34), collapse = " ")
pro <- pro[, c(1,53:60,2:52)]
p <- forestploter::forest(pro[, c(1:9)],
                          est = list(pro$FDX1_HR,pro$LIAS_HR,pro$LIPT1_HR,pro$DLD_HR, pro$DLAT_HR, pro$PHDA1_HR,
                                     pro$PHDB_HR, pro$SLC31A1_HR),
                          lower = list(pro$FDX1_Lower,pro$LIAS_Lower,pro$LIPT1_Lower, pro$DLD_Lower, pro$DLAT_Lower,
                                       pro$PHDA1_Lower, pro$PHDB_Lower, pro$SLC31A1_Lower),
                          upper = list(pro$FDX1_Upper,pro$LIAS_Upper,pro$LIPT1_Upper, pro$DLD_Upper, pro$DLAT_Upper,
                                       pro$PHDA1_Upper, pro$PHDB_Upper, pro$SLC31A1_Upper),
                          ref_line = 1,xlim = c(0,2.5), ticks_at = c(0,1,2),
                          ci_column = (2:9),theme = tm)
plot(p)
ggsave(filename =paste0(dir, "Graph/Forestplot_","OS",".svg"), plot = p,width = 16, height = 8)

stage <- readxl::read_xlsx("Stage.xlsx", sheet = 6)
stage <- separate(stage, Sample, into = c("Sample", "dis"), sep = "[(]")
stage1 <- stage[which(stage$Sample %in% c("GBMLGG","KIRP")), ]
fac <- c("Stage I","Stage II", "Stage III","Stage IV")
stage1$Stage <- factor(stage1$Stage, levels = fac)
p1<-ggboxplot(stage1,
              x="Stage",
              y="SLC31A1",
              color="Stage",
              fill="Stage",
              add = "jitter", facet.by = "Sample",
              bxp.errorbar.width = 1,
              width = 0.4,
              size=0.5,
              font.label = list(size=30), 
              palette = "Set1")+theme(panel.background =element_blank())+
  theme(legend.position = "right")+
  stat_compare_means(method="wilcox.test",hide.ns = F,
                     comparisons =list(c("Stage I","Stage II"), c("Stage I","Stage III"),c("Stage III", "Stage IV"),
                                       c("Stage II", "Stage III"),c("Stage I","Stage IV")),
                     label="p.signif")
ggsave(filename =paste0(dir, "Graph/Stageplot_SLC31A1",".svg"), plot = p1,width = 7, height = 4)

num <- seq(1,15,2)
esR <- function(x){
  es <- readxl::read_xlsx("Estimate.xlsx", sheet = x)
  es <- es[-1,]
  est1 <- es[, c(1,4,8,12)]
  colnames(est1) <- c("Sample","StromalScore","ImmuneScore","EstimateScore")
  est1$Index <- rep("R",nrow(est1))
  est2 <- es[, c(1,5,9,13)]
  colnames(est2) <- c("Sample","StromalScore","ImmuneScore","EstimateScore")
  est2$Index <- rep("pvalue",nrow(est2))
  est <- rbind(est1,est2)
}
esl <- lapply(num, esR)
names(esl) <- c("FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "SLC31A1")
est <- data.frame()
for (i in 1:8) {
  esl[[i]]$Gene <- rep(names(esl[i]), nrow(esl[[i]]))
  est <- rbind(est, esl[[i]])
}
UN <- c("TARGET-LAML(N=142)","TCGA-SARC(N=258)","TCGA-THYM(N=118)","TCGA-SKCM-P(N=101)","TCGA-SKCM-M(N=351)","TARGET-NB(N=153)",
  "TCGA-MESO(N=85)","TCGA-UVM(N=79)","TCGA-DLBC(N=46)","TARGET-ALL-R(N=99)")
est <- est[-which(est$Sample %in% UN),]
est <- separate(est, Sample, into = c("Dis", "Sample"), sep = "[-]")
est <- est[,-1]
est1 <- gather(est,key = ScoreCa, value = Score, -Sample, -Gene,-Index)
est1 <- spread(est1, key=Index, value = Score)
est1$Gene <- factor(est1$Gene, levels = fac)
est1[,c(4,5)] <- lapply(est1[,c(4,5)], as.numeric)
est1$p.adj <- -log10(est1$pvalue+0.05)
est1[which(est1$p.adj <1 & est1$p.adj>0.25), "p.adj"] <- 0.1
est1[which(est1$p.adj <=1.25 & est1$p.adj>=1), "p.adj"] <- 1
est1 <- unite(est1, "Gene", ScoreCa, Gene, sep = "_")
cols<-c('#523EFF',"#58D6A8",'#F5BC48','#ED2453')
pal<-colorRampPalette(cols)
p <- ggplot(est1, aes(Sample,Gene,fill=R,size =p.adj)) +  #划线
  geom_point(shape=21)+ #设置点
  scale_fill_gradientn(colours = pal(20))+
  theme_classic()+
  labs(title = "Score")+
  theme(
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=10,face="plain",color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1,  vjust = 1),
    legend.title=element_text(size=12,face="plain",color="black"),
    legend.background = element_blank(),legend.key.size =unit(0.5, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggsave(filename =paste0(dir, "Graph/EstimatePlot",".svg"), plot = p,width = 11, height = 8)


seq <- seq(4,268,4)
seq1 <- seq+1
seq2 <- seq(1,265,4)
name <- c()
for (i in seq2) {
  n <- unlist(strsplit(colnames(xcr),split = "[....]"))[i]
  name <- c(name,n)}
name <- c("Sample",name)
xcR <- function(x){
  xCell <- readxl::read_xlsx("xCell.xlsx", sheet = x)
  xCell <- xCell[-1,]
  xcr <- xCell[,c(1,seq)]
  xcp <- xCell[,c(1,seq1)]
  colnames(xcr) <- name
  colnames(xcp) <- name
  xcr$Index <- rep("R",nrow(xcr))
  xcp$Index <- rep("pvalue", nrow(xcp))
  xc <- rbind(xcr,xcp) %>% return()
  }
xc <- lapply(1:8, xcR)
names(xc) <- c("FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "SLC31A1")
xct <- data.frame()
for (i in 1:8) {
  xc[[i]]$Gene <- rep(names(xc[i]), nrow(xc[[i]]))
  xct <- rbind(xct, xc[[i]])
}
xct <- xct[-which(xct$Sample %in% UN),]
colnames(xct)
cell <- c("Sample","Gene","Index","CLP",
          "NKT","Th1_cells","Th2_cells")

xct1 <- xct[,which(colnames(xct) %in% cell)]
xct2 <- gather(xct1, key = Cell, value = values, -Sample, -Gene, -Index)
xct2 <- spread(xct2, key = Index, value = values)
xct2[,c(4,5)] <- lapply(xct2[,c(4,5)], as.numeric)
xct2$p.adj <- -log10(xct2$pvalue+0.05)
xct2[which(xct2$p.adj <1 & xct2$p.adj>0), "p.adj"] <- 0.1
xct2[which(xct2$p.adj <=1.25 & xct2$p.adj>=1), "p.adj"] <- 1
write.xlsx(unique(xct2$Gene), "fac1.xlsx", row.names = F)
fac1 <- readxl::read_xlsx("fac1.xlsx", col_names = F)
fac1 <- c(fac1[,1])
fac1 <- fac1[[1]]
xct2 <- unite(xct2, "Gene", Cell,Gene,  sep = "_")
xct2$Gene <- factor(xct2$Gene, levels = fac1, labels = fac1)

p <- ggplot(xct2, aes(Gene,Sample,fill=R,size =p.adj)) +  #划线
  geom_point(shape=21)+ #设置点
  scale_fill_gradientn(colours = pal(20))+
  theme_classic()+
  # labs(title = "Score")+
  theme(
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=10,face="plain",color="black"),
    axis.text.x = element_text(angle = 90, hjust = 1,  vjust = 1),
    legend.title=element_text(size=12,face="plain",color="black"),
    legend.background = element_blank(),legend.key.size =unit(0.5, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggsave(filename =paste0(dir, "Graph/xCellPlot",".svg"), plot = p, width = 10, height = 9)

gene <- c("FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "SLC31A1")
sm <- read.csv("GDC-PANCAN.mutect2_snv.tsv", header = T, sep = "\t")
sm1 <- sm[which(sm$gene %in% gene), ]
sm1 <- sm1[,c(1,2,9,8,3)]
colnames(sm1)=c("sample","gene", "variant_class","amino.acid.change","chr")
sm1 <- separate(sm1, variant_class, into = c("variant_class", "dis"),sep = "[;]")
sm1 <- sm1[,-4]
pdf(paste0(dir, "Graph/WaterfallPlot",".pdf"), width = 10, height = 8)
waterfall(sm1, fileType="Custom",variant_class_order = unique(sm1$variant_class))
dev.off()

unique(sm1$variant_class)

A <- unique(stage$Sample)
for (i in 1:37) {
  dat <- stage[which(stage$Sample == A[i]), ]
  t <- compare_means(SLC31A1~Stage, data = dat, method = "wilcox.test", 
                     paired = F) 
  print(t$p.signif)
}

drug <- readxl::read_xlsx("DTP_NCI60_ZSCORE.xlsx")
drug <- drug[-c(1:7),]
colnames(drug) <- drug[1,]
drug <- drug[-1, ]
table(drug$`FDA status`)
drug <- drug[drug$`FDA status` %in% c("FDA approved", "Clinical trial"), ]
colna <- as.data.frame(which(drug=="na", arr.ind = T))
drug <- drug[,-c(65,66)]
table(colna$col)
drug <- drug[,-c(40,42)]
exp <- readxl::read_xls("RNA__RNA_seq_composite_expression.xls")
exp <- exp[-c(1:9),]
colnames(exp) <- exp[1,]
exp <- exp[-1,]
exp <- exp[,-c(40,42)]
library(impute)
drug <- drug[,-1]
drug <- as.matrix(drug)
row.names(drug) <- drug[,1]
drugm <- drug[,-(1:5)]
dimnames <- list(rownames(drug), colnames(drugm))
data <- matrix(as.numeric(as.matrix(drugm)), nrow = nrow(drugm), dimnames = dimnames) 
data_1 <- impute.knn(data)
drug <- data_1$data
drug <- avereps(drug)
name <- c("FDX1","LIAS","LIPT1","DLD","DLAT","PDHA1","PDHB","SLC31A1")
exp <- exp[exp$`Gene name d` %in% name, ]
exp <- as.data.frame(exp)
rownames(exp) <- exp$`Gene name d`
exp <- exp[,-c(1:6)]

outTab <- data.frame()
for (i in rownames(exp)) {
  x <- as.numeric(exp[i,])
  for (j in row.names(drug)) {
    y <- as.numeric(drug[j, ])
    corT <- cor.test(x,y,method="spearman")
    cor <- corT$estimate
    pvalue <- corT$p.value
    if (pvalue<0.05&cor>0.4) {
      outVector <- cbind(i,j,cor,pvalue)
      outTab <- rbind(outTab, outVector)
      
    }
  }
}
outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))), ]
write.xlsx(outTab, "drug_gene_correlation_R0.4.xlsx", row.names = F)
outTab <- outTab[outTab$cor>0.44,]
corPlotNum <- 8 
plotList1 <- list()
for(i in 1:corPlotNum){
  Gene <- outTab[i, 1] 
  Drug <- outTab[i, 2] 
  x <- as.numeric(exp[Gene,]) 
  y <- as.numeric(drug[Drug,]) 
  cor <- sprintf("%.03f", as.numeric(outTab[i, 3])) 
  pvalue= 0 
  if(as.numeric(outTab[i, 4])< 0.001){ 
    pvalue= "p<0.001" 
  }else{ 
    pvalue=paste0( "p=",sprintf( "%.03f", as.numeric(outTab[i, 4])))
  } 
  df1 <- as.data.frame(cbind(x,y)) 
  p1=ggplot(data= df1, aes(x = x, y = y))+ 
    geom_point(size= 1.5, color="#ff4a12")+ 
    stat_smooth(method= "lm",se=FALSE, formula=y~x, color = "#0a008a",size=1.3)+ 
    labs(x= "Expression",y= "IC50",title = paste0(Gene, ", ",Drug),
         subtitle = paste0( "Cor=",cor, ", ",pvalue))+ 
    theme(axis.ticks = element_blank, 
          axis.text.y = element_blank,
          axis.text.x = element_blank)+ 
    theme_bw()
  plotList1[[i]]=p1
}

plotList2 <- list() 
corPlotNum<- 8 
for(i in 1:corPlotNum){ 
  Gene <- outTab[i, 1] 
  Drug <- outTab[i, 2] 
  x <- as.numeric(exp[Gene,]) 
  y <- as.numeric(drug[Drug,]) 
  df1 <- as.data.frame(cbind(x,y)) 
  colnames(df1)[ 2] <- "IC50" 
  df1$group <- ifelse(df1$x > median(df1$x), "high", "low") 
  compaired <- list(c( "low", "high")) 
  p1 <- ggboxplot(df1, x = "group", y = "IC50", fill = "group", 
                  palette = c("#FC270F","#522F8C"), add = "jitter", 
                  size = 0.5, xlab = paste0( "expression of ", Gene), 
                  ylab = paste0( "IC50 of ", Drug),width = 0.5) + 
    stat_compare_means(comparisons = compaired, method = "wilcox.test",
                       symnum.args= list(cutpoints = c( 0, 0.001, 0.01, 0.05, 1), 
                                         symbols = c( "***", "**", "*", "ns"))) 
  plotList2[[i]]=p1 
}
p <- ggarrange(plotlist= plotList1,nrow=2,ncol=4)
ggsave(filename =paste0(dir, "Graph/drugCorPlot",".svg"), plot = p, width = 8, height = 4)
p <- ggarrange(plotlist= plotList2,nrow=2,ncol=4)
ggsave(filename =paste0(dir, "Graph/drugPlot",".svg"), plot = p, width = 8.5, height = 8) 

# clin <- read.csv("Survival_SupplementalTable_S1_20171025_xena_sp.csv", header = T, sep = "\t")
# clin <- clin[clin$cancer.type.abbreviation == "KIRC", ]
# gene <- readxl::read_xlsx("Results.xlsx")
# ls <- clin$sample
# gene <- gene[gene$SampleName %in% ls, ]
# gene <- gene[gene$Sample == "KIRC(T=530,N=168)", ]
# colnames(gene)[1] <- "sample"
# int <- merge(clin, gene,by = "sample" ,all = T)
# int <- int[, -c(8,12,18:21,24:25, 34:35)]
# os <- int[, c(18:19, 27:34)]
# os <- na.omit(os)
# library(survival)
# library(survivalROC)
# library(survminer)
# os$OS <- as.numeric(os$OS)
# sfit <- coxph(Surv(OS.time,OS)~LIPT1, data = os)
# 
# summary(sfit)
# sfit <- coxph(Surv(OS.time,OS)~LIAS+FDX1, data = os)
# pre <- predict(sfit, type = "lp")
# os$pre <- predict(sfit, type = "lp")
# coxROC <- survivalROC(Stime = int1$OS.time, status = int1$OS, marker = pre, predict.time = 30, method = "KM")
# plot(coxROC$FP,coxROC$TP,type="l", xlab=round(coxROC$AUC,2))
# int1 <- int[, c(4,5,7,9,18,19,27:34)]
# int1 <- na.omit(int1)
# sfit <- coxph(Surv(OS.time,OS)~FDX1+LIPT1+LIAS+PDHA1, data = int1)
# int1$pre <- predict(sfit, type = "lp")
# int1[,c(2:4)] <- lapply(int1[,c(2:4)], as.factor)
# int2 <- int1[,-c(7:14)]
# sfit <- coxph(Surv(OS.time,OS)~., data = int2)
# summary(sfit)
# sfit <- survfit(Surv(OS.time,OS)~group, data = os)
# os$group <- ifelse(os$pre > median(os$pre), "high", "low") 
# ggsurvplot(sfit, conf.int=F, pval=TRUE)
# ggsurvplot(sfit,
#            # xlim=c(0,1100),
#            pval = T, ##添加P值
#            conf.int = TRUE,# 显示置信区间
#            conf.int.style='ribbon', ##设置置信区间风格
#            linetype = "strata", # 根据性别分组自动设置曲线类型
#            surv.median.line = "hv", # 设置中位生存期显示
#            ggtheme = theme_classic(), # 设置ggplot2主题
#            legend.labs=c('F','M'),##改变图例标签
#            legend.title='Sex',
#            legend = c(0.80,0.88), # 指定图例位置
#            palette = c("red", "navy"))


