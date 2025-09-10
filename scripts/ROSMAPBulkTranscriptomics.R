library(limma)
library(pheatmap)

setwd("~/ESRRA_Rerun/Bulk_ROSMAP/rawData")
phenotype=read.table("Phenotype.txt",header=T,row.names=1,sep="\t") #the data download from synapse database
trait=phenotype[,c(13,15,16,17)]
t=pheatmap(trait,clustering_method="ward.D2",show_rownames=F)
a=data.frame(cutree(t$tree_row,k=2))
anno=data.frame(Group=paste0("G",a[,1],sep=""))
rownames(anno)=rownames(a)
anno$Gender=phenotype$msex
anno$Status=ifelse(anno$Group=="G1","LOAD","ND")
anno$Group=NULL
anno$Age=phenotype$age_at_visit_max
ann_colors = list(
  Status = c(LOAD="Coral", ND="SlateBlue"),
  Gender=c(Female="RoyalBlue",Male="Violet")
)

##Figure 1A heatmap
pdf("1_HeatmapTraits.pdf")
pheatmap(trait,clustering_method="ward.D2", color = viridis(7),show_rownames=F,annotation_row=anno,annotation_colors = ann_colors)
dev.off()

pdf("1a_HeatmapTraits.pdf",height=5,width=4.5)
print(t)
dev.off()

#downloaded from Syn3388564/ROSMAP_RNAseq_FPKM_gene.tsv
GeneExprByFPKMFromROSMAP=read.table("geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
sampleList=intersect(rownames(phenotype),colnames(GeneExprByFPKMFromROSMAP))
length(sampleList)
#613
GeneExprByFPKMFromROSMAP=GeneExprByFPKMFromROSMAP[,sampleList]
all(rownames(phenotype)==colnames(GeneExprByFPKMFromROSMAP))
#TRUE

#### Significance Calculations ####
range(GeneExprByFPKMFromROSMAP)

data <- log2(GeneExprByFPKMFromROSMAP+1)
range(data)
phenotype$patientID = row.names(phenotype)
anno$patientID = row.names(anno)
phenotype=merge(phenotype, anno, by.x="patientID",by.y="patientID")
row.names(phenotype) = phenotype$patientID

phenotype.Male=phenotype[phenotype$msex=="Male", c("Status","pmi","age_at_visit_max")]
phenotype.Male[is.na(phenotype.Male)] <- 0
phenotype.Male$Status=factor(phenotype.Male$Status,levels=c("ND","LOAD"))
data.Male=data[,rownames(phenotype.Male)]
design.Male=model.matrix(~Status+pmi+age_at_visit_max,data =phenotype.Male)
fit.Male=lmFit(data.Male,design.Male)
fit.Male=eBayes(fit.Male)
top.Male=topTable(fit.Male, coef = "StatusLOAD",number=Inf,p.value=1, lfc=0)
#these soft cutoffs were accepted to idnetifed a number of DEGs in Male
top.Male$Pattern=ifelse(top.Male$P.Value < 0.01 & abs(top.Male$logFC) >= 0, ifelse(top.Male$logFC >= 0 ,'Up','Down'),'NoSig')
write.table(top.Male,file="Male_DEG_BulkROSMAP.txt",sep="\t",quote=F)
SigDEG.Male=top.Male[top.Male$Pattern%in%c("Up","Down"),]
table(top.Male$Pattern)
#Down NoSig    Up 
#264 49735  175
write.table(SigDEG.Male,file="Male_SigDEG_BulkROSMAP.txt",sep="\t",quote=F)

phenotype.Female=phenotype[phenotype$msex=="Female",c("Status","pmi","age_at_visit_max")]
phenotype.Female[is.na(phenotype.Female)] <- 0
phenotype.Female$Status=factor(phenotype.Female$Status,levels=c("ND","LOAD"))
data.Female=data[,rownames(phenotype.Female)]
design.Female=model.matrix(~Status+pmi+age_at_visit_max,data =phenotype.Female)
fit.Female=lmFit(data.Female,design.Female)
fit.Female=eBayes(fit.Female)
top.Female=topTable(fit.Female, coef = "StatusLOAD",number=Inf,p.value=1, lfc=0)
top.Female$Pattern=ifelse(top.Female$P.Value < 0.01 & abs(top.Female$logFC) >= 0, ifelse(top.Female$logFC >= 0,'Up','Down'),'NoSig')
write.table(top.Female,file="Female_DEG_BulkROSMAP.txt",sep="\t",quote=F)

SigDEG.Female=top.Female[top.Female$Pattern%in%c("Up","Down"),]
table(top.Female$Pattern)
#Down NoSig    Up 
#5215 43559  1490
write.table(SigDEG.Female,file="Female_SigDEG_BulkROSMAP.txt",sep="\t",quote=F)

phenotype.LOAD=phenotype[phenotype$Status=="LOAD",c("Gender","pmi","age_at_visit_max")]
phenotype.LOAD[is.na(phenotype.LOAD)] <- 0
phenotype.LOAD$Gender=factor(phenotype.LOAD$Gender,levels=c("Male","Female"))
data.LOAD=data[,rownames(phenotype.LOAD)]
design.LOAD=model.matrix(~Gender+pmi+age_at_visit_max,data =phenotype.LOAD)
fit.LOAD=lmFit(data.LOAD,design.LOAD)
fit.LOAD=eBayes(fit.LOAD)
top.LOAD=topTable(fit.LOAD, coef = "GenderFemale",number=Inf,p.value=1, lfc=0)
top.LOAD$Pattern=ifelse(top.LOAD$P.Value < 0.01 & abs(top.LOAD$logFC) >= 0, ifelse(top.LOAD$logFC >= 0,'Up','Down'),'NoSig')
write.table(top.LOAD,file="LOAD_DEG_BulkROSMAP.txt",sep="\t",quote=F)

SigDEG.LOAD=top.LOAD[top.LOAD$Pattern%in%c("Up","Down"),]
table(top.LOAD$Pattern)
#Down NoSig    Up 
#4632 44576 966
write.table(SigDEG.LOAD,file="LOAD_SigDEG_BulkROSMAP.txt",sep="\t",quote=F)

phenotype.ND=phenotype[phenotype$Status=="ND",c("Gender","pmi","age_at_visit_max")]
phenotype.ND[is.na(phenotype.ND)] <- 0
phenotype.ND$Gender=factor(phenotype.ND$Gender,levels=c("Male","Female"))
data.ND=data[,rownames(phenotype.ND)]
design.ND=model.matrix(~Gender+pmi+age_at_visit_max,data =phenotype.ND)
fit.ND=lmFit(data.ND,design.ND)
fit.ND=eBayes(fit.ND)
top.ND=topTable(fit.ND, coef = "GenderFemale",number=Inf,p.value=1, lfc=0)
top.ND$Pattern=ifelse(top.ND$P.Value < 0.01 & abs(top.ND$logFC) >= 0, ifelse(top.ND$logFC >= 0,'Up','Down'),'NoSig')
write.table(top.ND,file="ND_DEG_BulkROSMAP.txt",sep="\t",quote=F)

SigDEG.ND=top.ND[top.ND$Pattern%in%c("Up","Down"),]
table(top.ND$Pattern)
#Down NoSig    Up 
#406  349616  152
write.table(SigDEG.ND,file="ND_SigDEG_BulkROSMAP.txt",sep="\t",quote=F)

## Figure 1B
# ROSMAP: Male AD vs ND
hs_data=top.Male
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern)) +
  geom_point(alpha=0.4, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: Male, AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank())
pdf("2_MaleDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()

# ROSMAP: Female AD vs ND
hs_data=top.Female
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern)) +
  geom_point(alpha=0.4, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: Female, AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) 
pdf("3_FemaleDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()

hs_data=top.LOAD
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern)) +
  geom_point(alpha=0.6, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: LOAD, Female vs Male") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) 
pdf("4_LOADDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()

hs_data=top.ND
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern)) +
  geom_point(alpha=0.4, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: ND, Female vs Male") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank())
pdf("5_NDDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()

#### No Labels ####
##viocano plots for DEGs
hs_data=top.Male
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern, label =ID )) +
  geom_point(alpha=0.4, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: Male, AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
  ggrepel::geom_text_repel(
    data = subset(hs_data, hs_data$P.Value < 0.01),
    aes(label = ID),
    max.overlaps = 5,
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("2_MaleDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()

hs_data=top.Female
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern, label =ID )) +
  geom_point(alpha=0.4, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: Female, AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
  ggrepel::geom_text_repel(
    data = subset(hs_data, hs_data$P.Value < 0.01),
    aes(label = ID),
    max.overlaps = 5,
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("3_FemaleDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()

# Supplementary Figure 1A: Female ND vs Male ND
hs_data=top.ND
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern, label =ID )) +
  geom_point(alpha=0.4, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: ND, Female vs Male") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
  ggrepel::geom_text_repel(
    data = subset(hs_data, hs_data$P.Value < 0.01),
    aes(label = ID),
    max.overlaps = 5,
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("5_NDDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()

# Supplementary Figure 1B: Female AD vs Male AD
hs_data=top.LOAD
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(P.Value), colour=Pattern, label =ID )) +
  geom_point(alpha=0.4, aes(size=-log10(P.Value))) +
  theme_bw() + 
  scale_size(range = c(0,5))+
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene: LOAD, Female vs Male") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
  ggrepel::geom_text_repel(
    data = subset(hs_data, hs_data$P.Value < 0.01),
    aes(label = ID),
    max.overlaps = 5,
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("4_LOADDEG_Volcano.pdf",height=7,width=9)
print(t)
dev.off()



###Figure 4G, and Figure 6E gene expression visualiation, boxplots #############
GeneExpr=read.table("geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
SampleInfo=read.table("Phenotype.txt",header=TRUE,row.names=1,sep="\t")
Samples=intersect(colnames(GeneExpr),rownames(SampleInfo))
GeneExpr=GeneExpr[rowSums(GeneExpr)>0,Samples]
SampleInfo=SampleInfo[Samples,]
all(colnames(GeneExpr)==rownames(SampleInfo))
BigData=cbind(SampleInfo,t(GeneExpr))

gene=c("PDHA1","DHCR24","VLDLR","LRP1","SDHA","SDHD","CYP51A1")
BigDataCogdx=BigData[BigData$cogdx %in% c(1,2,4,5),]

tmp=BigDataCogdx[,c("msex","cogdx",gene)]
library(reshape2)
tmp1=melt(tmp,id=c(1,2))
colnames(tmp1)=c("Gender","cogdx","Gene","Expression")

g=ggplot(tmp1, aes(x= as.character(cogdx), y=log2(Expression+1), color=as.character(cogdx))) +  
  #geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.7)+
  #coord_flip() + 
  facet_grid( Gene~Gender, scales="free_y")+
  theme_bw() + geom_jitter(position=position_jitter(0.2),size=1)+
  scale_color_brewer(palette="Set2")+ggpubr::stat_compare_means()+
  #theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme (strip.background = element_blank(), strip.text.x = element_text(size=rel(1.0)),strip.text.y = element_text(size=rel(1.0),face="italic")) + 
  theme (axis.text.x =element_text(size=rel(1.0)), axis.ticks = element_blank()) + 
  guides(fill=FALSE) +
  scale_x_discrete(expand = c(0.2, 0)) + 
  scale_y_continuous(expand = c(0.1, 0.5))
pdf("Cogdx_PDHA1_Expr_pvalue.pdf",height=12,width=6)
print(g)
dev.off()

# PDHA1
PDHA_1 = tmp %>% 
  filter(cogdx == 1) %>%
  select(msex, cogdx, Expression = PDHA1)

pdha1 <- ggplot(PDHA_1, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(PDHA1 Expression + 1)",
    title = "PDHA1 Expression in Individuals with cogdx = 1"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("1_PDHA_1.pdf",height=6,width=6)
print(pdha1)
dev.off()

PDHA_2 = tmp %>% 
  filter(cogdx == 2) %>%
  select(msex, cogdx, Expression = PDHA1)

pdha2 <- ggplot(PDHA_2, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(PDHA1 Expression + 1)",
    title = "PDHA1 Expression in Individuals with cogdx = 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("2_PDHA_2.pdf",height=6,width=6)
print(pdha2)
dev.off()

PDHA_4 = tmp %>% 
  filter(cogdx == 4) %>%
  select(msex, cogdx, Expression = PDHA1)

pdha4 <- ggplot(PDHA_4, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(PDHA1 Expression + 1)",
    title = "PDHA1 Expression in Individuals with cogdx = 4"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("3_PDHA_4.pdf",height=6,width=6)
print(pdha4)
dev.off()

PDHA_5 = tmp %>% 
  filter(cogdx == 5) %>%
  select(msex, cogdx, Expression = PDHA1)

pdha5 <- ggplot(PDHA_5, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(PDHA1 Expression + 1)",
    title = "PDHA1 Expression in Individuals with cogdx = 5"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("4_PDHA_5.pdf",height=6,width=6)
print(pdha5)
dev.off()

# DHCR24
DHCR24_1 = tmp %>% 
  filter(cogdx == 1) %>%
  select(msex, cogdx, Expression = DHCR24)

dhcr241 <- ggplot(DHCR24_1, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(DHCR24 Expression + 1)",
    title = "DHCR24 Expression in Individuals with cogdx = 1"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("1_DHCR24_1.pdf",height=6,width=6)
print(dhcr241)
dev.off()

DHCR24_2 = tmp %>% 
  filter(cogdx == 2) %>%
  select(msex, cogdx, Expression = DHCR24)

dhcr242 <- ggplot(DHCR24_2, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(DHCR24 Expression + 1)",
    title = "DHCR24 Expression in Individuals with cogdx = 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("2_DHCR24_2.pdf",height=6,width=6)
print(dhcr242)
dev.off()

DHCR24_4 = tmp %>% 
  filter(cogdx == 4) %>%
  select(msex, cogdx, Expression = DHCR24)

dhcr244 <- ggplot(DHCR24_4, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(DHCR24 Expression + 1)",
    title = "DHCR24 Expression in Individuals with cogdx = 4"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("3_DHCR24_4.pdf",height=6,width=6)
print(dhcr244)
dev.off()

DHCR24_5 = tmp %>% 
  filter(cogdx == 5) %>%
  select(msex, cogdx, Expression = DHCR24)

dhcr245 <- ggplot(DHCR24_5, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(DHCR24 Expression + 1)",
    title = "DHCR24 Expression in Individuals with cogdx = 5"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("4_DHCR24_5.pdf",height=6,width=6)
print(dhcr245)
dev.off()

# VLDLR1
VLDLR_1 = tmp %>% 
  filter(cogdx == 1) %>%
  select(msex, cogdx, Expression = VLDLR)

vldlr1 <- ggplot(VLDLR_1, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(VLDLR1 Expression + 1)",
    title = "VLDLR1 Expression in Individuals with cogdx = 1"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("1_VLDLR_1.pdf",height=6,width=6)
print(vldlr1)
dev.off()

VLDLR_2 = tmp %>% 
  filter(cogdx == 2) %>%
  select(msex, cogdx, Expression = VLDLR)

vldlr2 <- ggplot(VLDLR_2, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(VLDLR1 Expression + 1)",
    title = "VLDLR1 Expression in Individuals with cogdx = 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("2_VLDLR_2.pdf",height=6,width=6)
print(vldlr2)
dev.off()

VLDLR_4 = tmp %>% 
  filter(cogdx == 4) %>%
  select(msex, cogdx, Expression = VLDLR)

vldlr4 <- ggplot(VLDLR_4, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(VLDLR1 Expression + 1)",
    title = "VLDLR1 Expression in Individuals with cogdx = 4"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("3_VLDLR_4.pdf",height=6,width=6)
print(vldlr4)
dev.off()

VLDLR_5 = tmp %>% 
  filter(cogdx == 5) %>%
  select(msex, cogdx, Expression = VLDLR)

vldlr5 <- ggplot(VLDLR_5, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(VLDLR1 Expression + 1)",
    title = "VLDLR1 Expression in Individuals with cogdx = 5"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("4_VLDLR_5.pdf",height=6,width=6)
print(vldlr5)
dev.off()

# LRP11
LRP1_1 = tmp %>% 
  filter(cogdx == 1) %>%
  select(msex, cogdx, Expression = LRP1)

lrp1 <- ggplot(LRP1_1, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(LRP11 Expression + 1)",
    title = "LRP11 Expression in Individuals with cogdx = 1"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("1_LRP1_1.pdf",height=6,width=6)
print(lrp1)
dev.off()

LRP1_2 = tmp %>% 
  filter(cogdx == 2) %>%
  select(msex, cogdx, Expression = LRP1)

lrp2 <- ggplot(LRP1_2, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(LRP11 Expression + 1)",
    title = "LRP11 Expression in Individuals with cogdx = 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("2_LRP1_2.pdf",height=6,width=6)
print(lrp2)
dev.off()

LRP1_4 = tmp %>% 
  filter(cogdx == 4) %>%
  select(msex, cogdx, Expression = LRP1)

lrp4 <- ggplot(LRP1_4, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(LRP11 Expression + 1)",
    title = "LRP11 Expression in Individuals with cogdx = 4"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("3_LRP1_4.pdf",height=6,width=6)
print(lrp4)
dev.off()

LRP1_5 = tmp %>% 
  filter(cogdx == 5) %>%
  select(msex, cogdx, Expression = LRP1)

lrp5 <- ggplot(LRP1_5, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(LRP11 Expression + 1)",
    title = "LRP11 Expression in Individuals with cogdx = 5"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("4_LRP1_5.pdf",height=6,width=6)
print(lrp5)
dev.off()

# SDHA
SDHA_1 = tmp %>% 
  filter(cogdx == 1) %>%
  select(msex, cogdx, Expression = SDHA)

shda1 <- ggplot(SDHA_1, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHA Expression + 1)",
    title = "SDHA Expression in Individuals with cogdx = 1"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("1_SDHA_1.pdf",height=6,width=6)
print(shda1)
dev.off()

SDHA_2 = tmp %>% 
  filter(cogdx == 2) %>%
  select(msex, cogdx, Expression = SDHA)

shda2 <- ggplot(SDHA_2, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHA Expression + 1)",
    title = "SDHA Expression in Individuals with cogdx = 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("2_SDHA_2.pdf",height=6,width=6)
print(shda2)
dev.off()

SDHA_4 = tmp %>% 
  filter(cogdx == 4) %>%
  select(msex, cogdx, Expression = SDHA)

shda4 <- ggplot(SDHA_4, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHA Expression + 1)",
    title = "SDHA Expression in Individuals with cogdx = 4"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("3_SDHA_4.pdf",height=6,width=6)
print(shda4)
dev.off()

SDHA_5 = tmp %>% 
  filter(cogdx == 5) %>%
  select(msex, cogdx, Expression = SDHA)

shda5 <- ggplot(SDHA_5, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHA Expression + 1)",
    title = "SDHA Expression in Individuals with cogdx = 5"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("4_SDHA_5.pdf",height=6,width=6)
print(shda5)
dev.off()

# SDHD
SDHD_1 = tmp %>% 
  filter(cogdx == 1) %>%
  select(msex, cogdx, Expression = SDHD)

sdhd1 <- ggplot(SDHD_1, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHD Expression + 1)",
    title = "SDHD Expression in Individuals with cogdx = 1"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("1_SDHD_1.pdf",height=6,width=6)
print(sdhd1)
dev.off()

SDHD_2 = tmp %>% 
  filter(cogdx == 2) %>%
  select(msex, cogdx, Expression = SDHD)

sdhd2 <- ggplot(SDHD_2, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHD Expression + 1)",
    title = "SDHD Expression in Individuals with cogdx = 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("2_SDHD_2.pdf",height=6,width=6)
print(sdhd2)
dev.off()

SDHD_4 = tmp %>% 
  filter(cogdx == 4) %>%
  select(msex, cogdx, Expression = SDHD)

sdhd4 <- ggplot(SDHD_4, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHD Expression + 1)",
    title = "SDHD Expression in Individuals with cogdx = 4"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("3_SDHD_4.pdf",height=6,width=6)
print(sdhd4)
dev.off()

SDHD_5 = tmp %>% 
  filter(cogdx == 5) %>%
  select(msex, cogdx, Expression = SDHD)

sdhd5 <- ggplot(SDHD_5, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(SDHD Expression + 1)",
    title = "SDHD Expression in Individuals with cogdx = 5"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("4_SDHD_5.pdf",height=6,width=6)
print(sdhd5)
dev.off()

# CYP51A1
CYP51A1_1 = tmp %>% 
  filter(cogdx == 1) %>%
  select(msex, cogdx, Expression = CYP51A1)

cyp511 <- ggplot(CYP51A1_1, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(CYP51A1 Expression + 1)",
    title = "CYP51A1 Expression in Individuals with cogdx = 1"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("1_CYP51A1_1.pdf",height=6,width=6)
print(cyp511)
dev.off()

CYP51A1_2 = tmp %>% 
  filter(cogdx == 2) %>%
  select(msex, cogdx, Expression = CYP51A1)

cyp512 <- ggplot(CYP51A1_2, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(CYP51A1 Expression + 1)",
    title = "CYP51A1 Expression in Individuals with cogdx = 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("2_CYP51A1_2.pdf",height=6,width=6)
print(cyp512)
dev.off()

CYP51A1_4 = tmp %>% 
  filter(cogdx == 4) %>%
  select(msex, cogdx, Expression = CYP51A1)

cyp514 <- ggplot(CYP51A1_4, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(CYP51A1 Expression + 1)",
    title = "CYP51A1 Expression in Individuals with cogdx = 4"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("3_CYP51A1_4.pdf",height=6,width=6)
print(cyp514)
dev.off()

CYP51A1_5 = tmp %>% 
  filter(cogdx == 5) %>%
  select(msex, cogdx, Expression = CYP51A1)

cyp515 <- ggplot(CYP51A1_5, aes(x = msex, y = log2(Expression + 1), fill = msex)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.5) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  ggpubr::stat_compare_means(
    method = "t.test",  # or "wilcox.test" for non-parametric
    comparisons = list(c("Female", "Male")),  # Compare Female vs Male
    label = "p.signif"  # Show significance stars
  ) +
  labs(
    x = "Gender",
    y = "log2(CYP51A1 Expression + 1)",
    title = "CYP51A1 Expression in Individuals with cogdx = 5"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("4_CYP51A1_5.pdf",height=6,width=6)
print(cyp515)
dev.off()

#### Significance of each gene comparing the cognitive scores between male and female ####
library(dplyr)

#   Female Male
#1    119   76
#2     98   53
#4    147   66
#5     18   15

cogdx_groups <- c(1, 2, 4, 5)  
results <- list()

for (cogdx_val in cogdx_groups) {
  df_sub <- tmp %>% filter(cogdx == cogdx_val)
  
  # Skip if no data for either gender
  if (!all(c("Female", "Male") %in% df_sub$msex)) {
    warning(paste("Skipping cogdx =", cogdx_val, "due to missing group data"))
    next
  }
  
  gene_columns <- setdiff(names(df_sub), c("msex", "cogdx"))
  
  for (gene in gene_columns) {
    # Extract expression values for females and males
    females <- df_sub %>% filter(msex == "Female") %>% pull(gene)
    males <- df_sub %>% filter(msex == "Male") %>% pull(gene)
    
    # Skip if any group has no observations
    if (length(females) == 0 || length(males) == 0) next
    
    # Perform Welch's t-test (var.equal = FALSE by default)
    t_test <- t.test(females, males)
    
    # Store results
    results[[paste0("cogdx_", cogdx_val, "_", gene)]] <- data.frame(
      cogdx = cogdx_val,
      gene = gene,
      female_mean = mean(females, na.rm = TRUE),
      male_mean = mean(males, na.rm = TRUE),
      t_statistic = t_test$statistic,
      p_value = t_test$p.value,
      df = t_test$parameter  # degrees of freedom
    )
  }
}

# Combine results and adjust p-values
results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")

# Save to file
write.table(results_df, file = "Cognitivescores_ERaTarget.txt", sep = "\t", quote = FALSE, row.names = TRUE)
