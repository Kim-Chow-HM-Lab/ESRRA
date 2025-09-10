#### ROSMAP Metabolomics Analysis ####

#### Library loading ####
library(plyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggplot2)

#### Preprocessing ####
setwd("~/Metabolomics/RawFiles")

FullClinicalMeta = read.csv("ROSMAP_clinical.csv", header = T, row.names = 1)
# 3584 patients included in this full clinical meta file for all the patients enrolled in ROSMAP

# Metabolite Processed Raw data is found in the .xlsx file downloaded on Synapse in the first sheet
Metab_Raw = read.table("ROSMAP_MetaboliteProcessed.txt", header=T, row.names=1, check.names=F)
PatientIDs = intersect(colnames(Metab_Raw), rownames(FullClinicalMeta))
length(PatientIDs)
# 500

MetaboliteInfo = read.table("MetaboliteInfo.txt",sep="\t",header=T)

#### Status Definitions ####
setwd("~/Metabolomics")
ClinicalMeta = FullClinicalMeta[PatientIDs,]
# 500 patients remaining from selecting the information for the 500 patients included in the metabolite raw file
ClinicalMeta_ADScores = ClinicalMeta[,c("braaksc","ceradsc","cogdx","dcfdx_lv")]
ClinicalMeta_ADScores$ceradsc_mod = 5-ClinicalMeta_ADScores$ceradsc
# cerad scoring was adjusted because the smaller numbers mean AD pathology, while the bigger numbers mean no plaques (1 for definite, 2 for probable, 3 for possible, 4 for noAD). This way the heatmap doesn't look nice hence was adjusted.
ClinicalMeta_ADScores$ceradsc = NULL

t=pheatmap(ClinicalMeta_ADScores,clustering_method="ward.D2",show_rownames=F)
pdf("1_heatmap.pdf")
print(t)
dev.off()

sampleGroup=cutree(t$tree_row,k=2)
table(sampleGroup)
#sampleGroup
#1   2 
#277 223 

ClinicalMeta_Covariates = ClinicalMeta[,c("msex","educ","apoe_genotype","age_at_visit_max","pmi")]
ClinicalMeta_Covariates$Sex = ifelse(ClinicalMeta_Covariates$msex ==0, "Female","Male")
ClinicalMeta_Covariates$Age = as.numeric(ifelse(ClinicalMeta_Covariates$age_at_visit_max=="90+","90",ClinicalMeta_Covariates$age_at_visit_max))
ClinicalMeta_Covariates$Apoe = as.character(ClinicalMeta_Covariates$apoe_genotype)
ClinicalMeta_Covariates$Status = ifelse(sampleGroup==1, "ND","AD")
ClinicalMeta_Covariates$msex=NULL
ClinicalMeta_Covariates$age_at_visit_max=NULL
ClinicalMeta_Covariates$apoe_genotype=NULL

ApoeColor = brewer.pal(6,"GnBu")
names(ApoeColor) = sort(unique(ClinicalMeta_Covariates$Apoe))
HeatmapColors = list(
  Sex = c(Female="deeppink1", Male = "cadetblue3"),
  Age = c("white","tomato3"),
  educ = c("white","seagreen4"),
  pmi = c("white","darkmagenta"),
  Apoe = ApoeColor,
  Status = c(AD = "firebrick1",ND="blue1"))

# Supplementary Figure 14a
pdf("2_heatmap_fullcovariates.pdf")
pheatmap(ClinicalMeta_ADScores, clustering_method = "ward.D2",color = inferno(10), show_rownames = F, annotation_row=ClinicalMeta_Covariates, annotation_colors = HeatmapColors)
dev.off()

all(rownames(ClinicalMeta_Covariates)==rownames(ClinicalMeta_ADScores))
write.table(cbind(ClinicalMeta_Covariates,ClinicalMeta_ADScores), file = "ClinicalStatusInfo.txt",sep="\t",quote=F)

#### DEM identification: Function ####
DEM_Calculation <- function(metab_id, model_data) {
  readings <- cbind(model_data[,1:ncol(ClinicalMeta_Covariates)], model_data[,metab_id]) #ncol(clinical.female)==ncol(clinical.male)
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),] #remove the sample with missing value
  model <- glm(reading ~ Status+educ+pmi+Age+Sex,data = readings, family = gaussian)
  data.frame(
    metab_id = metab_id,
    effect = as.matrix(summary(model)$coefficients)["StatusAD","Estimate"],
    pval = as.matrix(summary(model)$coefficients)["StatusAD","Pr(>|t|)"]
  )
}

DEM_Calculation_Sex <- function(metab_id, model_data) {
  readings <- cbind(model_data[,1:ncol(ClinicalMeta_Covariates)], model_data[,metab_id]) #ncol(clinical.female)==ncol(clinical.male)
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),] #remove the sample with missing value
  model <- glm(reading ~ Status+educ+pmi+Age,data = readings, family = gaussian)
  data.frame(
    metab_id = metab_id,
    effect = as.matrix(summary(model)$coefficients)["StatusAD","Estimate"],
    pval = as.matrix(summary(model)$coefficients)["StatusAD","Pr(>|t|)"]
  )
}

all(colnames(Metab_Raw) == rownames(ClinicalMeta_ADScores))
# True
ClinicalMeta_ADScores[is.na(ClinicalMeta_ADScores)] = 0

#### All Samples - AD vs ND ####
MetabClinicalIntegrated =cbind(ClinicalMeta_Covariates, t(Metab_Raw))
MetabClinicalIntegrated$Status=factor(MetabClinicalIntegrated$Status,levels=c("ND","AD"))
# Create data frame with effects, pvals, and padj for each metabolite
DEM_All <- ldply(colnames(MetabClinicalIntegrated[,-1:-ncol(ClinicalMeta_Covariates)]), 
                 DEM_Calculation,
                 model_data = MetabClinicalIntegrated)
DEM_All$padj <- p.adjust(DEM_All$pval, method = "BH")
write.table(DEM_All,file="DEM_All.txt",sep="\t",row.names=F,quote=F)

# Volcano Plots
DEM_All$Pattern = ifelse(DEM_All$padj>0.05, "NotSig", ifelse(DEM_All$effect>0, "UpinAD", "DowninAD"))
DEM_All = merge(DEM_All, MetaboliteInfo, by.x = "metab_id", by.y = "Metabolite")
write.table(DEM_All, file = "DEM_All_Info.txt", sep="\t", quote=F, row.names=F)

a=ggplot(data = DEM_All, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.01),-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_All, padj<0.05),
    aes(label = Name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

a2=ggplot(data = DEM_All, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_All, padj<0.05),
    aes(label = Name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

a3=ggplot(data = DEM_All, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_All, padj<0.05),
    aes(label = NA),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

pdf("3_DEM_All_ADvND_Volcano.pdf",height=7, width=8)
print(a)
dev.off()

pdf("3a_DEM_All_ADvND_Volcano.pdf",height=7, width=8)
print(a2)
dev.off()

pdf("3b_DEM_All_ADvND_Volcano.pdf",height=7, width=8)
print(a3)
dev.off()

#### Female Samples - AD vs ND ####
MetabClinicalIntegrated_Female=MetabClinicalIntegrated[MetabClinicalIntegrated$Sex=="Female",] #352
Metab_RawFemales=Metab_Raw[,rownames(MetabClinicalIntegrated_Female)]
all(colnames(Metab_RawFemales)==rownames(MetabClinicalIntegrated_Female)) #True
MetabClinicalIntegrated_Female$Status=factor(MetabClinicalIntegrated_Female$Status,levels=c("ND","AD"))
summary(MetabClinicalIntegrated_Female[,1:ncol(ClinicalMeta_Covariates)])
#educ            pmi            Sex               Age            Apoe           Status  
#Min.   : 5.00   Min.   : 2.500   Length:352         Min.   :71.23   Length:352         ND:191  
#1st Qu.:13.00   1st Qu.: 5.308   Class :character   1st Qu.:87.08   Class :character   AD:161  
#Median :16.00   Median : 6.583   Mode  :character   Median :90.00   Mode  :character           
#Mean   :15.52   Mean   : 7.967                      Mean   :87.75                              
#3rd Qu.:18.00   3rd Qu.: 8.646                      3rd Qu.:90.00                              
#Max.   :25.00   Max.   :32.583                      Max.   :90.00                              
#               NA's   :2                                                                      
# Create data frame with effects, pvals, and padj for each metabolite
DEM_female <- ldply(colnames(MetabClinicalIntegrated_Female[,-1:-ncol(ClinicalMeta_Covariates)]), 
                    DEM_Calculation_Sex,
                    model_data = MetabClinicalIntegrated_Female)
DEM_female$padj <- p.adjust(DEM_female$pval, method = "BH")
write.table(DEM_female,file="DEM_Female.txt",sep="\t",row.names=F,quote=F)

# Volcano Plots
DEM_female$Pattern = ifelse(DEM_female$padj>0.05, "NotSig", ifelse(DEM_female$effect>0, "UpinAD_female", "DowninAD_female"))
DEM_female = merge(DEM_female, MetaboliteInfo, by.x = "metab_id", by.y = "Metabolite")
write.table(DEM_female, file = "DEM_female_Info.txt",sep="\t",quote=F,row.names = F)

f=ggplot(data = DEM_female, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.01),-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In female: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_female, padj<0.05),
    aes(label = Name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

f2=ggplot(data = DEM_female, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In female: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_female, padj<0.05),
    aes(label = Name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

f3=ggplot(data = DEM_female, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In female: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_female, padj<0.05),
    aes(label = NA),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

pdf("4_DEM_female_Volcano.pdf",height=7, width=8)
print(f)
dev.off()

pdf("4a_DEM_female_Volcano.pdf",height=7, width=8)
print(f2)
dev.off()

pdf("4b_DEM_female_Volcano.pdf",height=7, width=8)
print(f3)
dev.off()

#### Male Samples - AD vs ND ####
MetabClinicalIntegrated_Male=MetabClinicalIntegrated[MetabClinicalIntegrated$Sex=="Male",] #148
Metab_RawMales=Metab_Raw[,rownames(MetabClinicalIntegrated_Male)]
all(colnames(Metab_RawMales)==rownames(MetabClinicalIntegrated_Male)) #True
MetabClinicalIntegrated_Male$Status=factor(MetabClinicalIntegrated_Male$Status,levels=c("ND","AD"))
summary(MetabClinicalIntegrated_Male[,1:ncol(ClinicalMeta_Covariates)])
#educ            pmi             Sex               Age            Apoe           Status 
#Min.   : 8.00   Min.   : 0.8667   Length:148         Min.   :72.30   Length:148         ND:86  
#1st Qu.:14.00   1st Qu.: 5.0667   Class :character   1st Qu.:84.46   Class :character   AD:62  
#Median :16.00   Median : 6.8500   Mode  :character   Median :88.93   Mode  :character          
#Mean   :16.77   Mean   : 8.2136                      Mean   :86.72                             
#3rd Qu.:19.00   3rd Qu.: 8.8125                      3rd Qu.:90.00                             
#Max.   :30.00   Max.   :43.8833                      Max.   :90.00                                                                                                 
# Create data frame with effects, pvals, and padj for each metabolite
DEM_Male <- ldply(colnames(MetabClinicalIntegrated_Male[,-1:-ncol(ClinicalMeta_Covariates)]), 
                    DEM_Calculation_Sex,
                    model_data = MetabClinicalIntegrated_Male)
DEM_Male$padj <- p.adjust(DEM_Male$pval, method = "BH")
write.table(DEM_Male,file="DEM_Male.txt",sep="\t",row.names=F,quote=F)

# Volcano Plots
DEM_Male$Pattern = ifelse(DEM_Male$padj>0.05, "NotSig", ifelse(DEM_Male$effect>0, "UpinAD_Male", "DowninAD_Male"))
DEM_Male = merge(DEM_Male, MetaboliteInfo, by.x = "metab_id", by.y = "Metabolite")
write.table(DEM_Male, file = "DEM_Male_Info.txt",sep="\t",quote=F,row.names = F)

m=ggplot(data = DEM_Male, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.01),-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In Male: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_Male, padj<0.05),
    aes(label = Name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

m2=ggplot(data = DEM_Male, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In Male: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_Male, padj<0.05),
    aes(label = Name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

m3=ggplot(data = DEM_Male, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =Name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In Male: AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEM_Male, padj<0.05),
    aes(label = NA),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

pdf("5_DEM_Male_Volcano.pdf",height=7, width=8)
print(m)
dev.off()

pdf("5a_DEM_Male_Volcano.pdf",height=7, width=8)
print(m2)
dev.off()

pdf("5b_DEM_Male_Volcano.pdf",height=7, width=8)
print(m3)
dev.off()

#### Permutation Analysis ####
# To account for the different numbers in the male vs female samples
table(MetabClinicalIntegrated$Sex, MetabClinicalIntegrated$Status)
#       ND  AD
#Female 191 161
#Male    86  62

PermutationNumber = 50
# Initialize a data frame to store ALL results (not just counts)
perm_results <- data.frame(
  Permutation = integer(),
  Metabolite = character(),
  effect = numeric(),
  pval = numeric(),
  padj = numeric(),
  Pattern = character(),
  stringsAsFactors = FALSE
)

# Initialize DEGNumber (now a data frame instead of matrix)
DEGNumber <- data.frame(
  Permutation = 1:PermutationNumber,
  UpDEMNumber = rep(NA, PermutationNumber),
  DownDEMNumber = rep(NA, PermutationNumber)
)

for (i in 1:PermutationNumber) {
  # Sample and subset data
  ADDownsample <- sample(rownames(MetabClinicalIntegrated[MetabClinicalIntegrated$Sex == "Female" & MetabClinicalIntegrated$Status == "AD", ]), 62, replace = FALSE)
  NDDownsample <- sample(rownames(MetabClinicalIntegrated[MetabClinicalIntegrated$Sex == "Female" & MetabClinicalIntegrated$Status == "ND", ]), 85, replace = FALSE)
  DownsampledFemales <- MetabClinicalIntegrated[c(ADDownsample, NDDownsample), ]
  
  # Run differential analysis
  DEMFemaleDownsampled <- ldply(colnames(metabolite_data), 
                                DEM_Calculation_Sex,
                                model_data = DownsampledFemales)
  DEMFemaleDownsampled$padj <- p.adjust(DEMFemaleDownsampled$pval, method = "BH")
  DEMFemaleDownsampled$Pattern <- ifelse(
    DEMFemaleDownsampled$padj < 0.05 & abs(DEMFemaleDownsampled$effect) >= 0,
    ifelse(DEMFemaleDownsampled$effect >= 0, "Up", "Down"),
    "NonSig"
  )
  
  # Store counts
  DEGNumber$UpDEMNumber[i] <- sum(DEMFemaleDownsampled$Pattern == "Up")
  DEGNumber$DownDEMNumber[i] <- sum(DEMFemaleDownsampled$Pattern == "Down")
  
  # Append full results for this permutation
  DEMFemaleDownsampled$Permutation <- i
  perm_results <- rbind(perm_results, DEMFemaleDownsampled)
}

# Save results
write.table(DEGNumber, "Permutation_DEM_Counts.txt",sep="\t",quote=F, row.names = FALSE)
write.table(perm_results, "DEM_PermutationFemale_Results.txt", sep="\t", quote=F, row.names=F)

# Volcano Plot for downsampled female samples
f4=ggplot(data = DEMFemaleDownsampled, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =metab_id)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,4)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.01),-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In Female (Downsampled): AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEMFemaleDownsampled, padj<0.05),
    aes(label = metab_id),
    size = 3,
    max.overlaps=18,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

f5=ggplot(data = DEMFemaleDownsampled, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =metab_id)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,4)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In Female (Downsampled): AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEMFemaleDownsampled, padj<0.05),
    aes(label = metab_id),
    size = 3,
    max.overlaps=18,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

f6=ggplot(data = DEMFemaleDownsampled, aes(x = effect, y = -log10(padj), size=-log10(padj),colour=Pattern, label =metab_id)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,4)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-adj)",title="DEM In Female (Downsampled): AD vs ND") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(DEMFemaleDownsampled, padj<0.05),
    aes(label = NA),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
  )

pdf("6_DEM_FemaleDownsampled_Volcano.pdf",height=7, width=8)
print(f4)
dev.off()

pdf("6a_DEM_FemaleDownsampled_Volcano.pdf",height=7, width=8)
print(f5)
dev.off()

pdf("6b_DEM_FemaleDownsampled_Volcano.pdf",height=7, width=8)
print(f6)
dev.off()

observed_up_DEMs <- sum(DEM_female$Pattern=="UpinAD_female")
observed_down_DEMs <- sum(DEM_female$Pattern=="DowninAD_female")

pdf("7_Histogram_PermutatedCounts.pdf")
ggplot(DEGNumber, aes(x = UpDEMNumber)) +
  geom_histogram(bins = 15, fill = "white", color = "black") +
  geom_vline(xintercept = observed_up_DEMs, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Permuted Up-DEMs Distribution", x = "Number of Up-DEMs", y = "Frequency")
dev.off()

pdf("7a_Histogram_PermutatedCounts.pdf")
ggplot(DEGNumber, aes(x = DownDEMNumber)) +
  geom_histogram(bins = 15, fill = "white", color = "black") +
  geom_vline(xintercept = observed_down_DEMs, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Permuted Down-DEMs Distribution", x = "Number of Down-DEMs", y = "Frequency")
dev.off()

# Supplementary Figure 14d
pdf("8_ECDFPlotforSignificance.pdf")
ggplot(perm_results, aes(x = -log10(padj), color = as.factor(Permutation))) +
  stat_ecdf(geom = "step", alpha = 0.3) +
  labs(title = "ECDF of Permuted P-values", x = "-Log10(Adjusted P-value)", y = "Cumulative Probability")
dev.off()

# Count occurrences of each metab_id as significant (Up/Down) across permutations
dem_frequency <- perm_results %>%
  filter(Pattern %in% c("Up", "Down")) %>%  
  group_by(metab_id) %>%
  summarize(
    Total_Significant = n(),  
    Frequency_Pct = (n() / PermutationNumber) * 100,  
    Avg_Effect = mean(effect),  
    Avg_padj = mean(padj)     
  ) %>%
  arrange(desc(Total_Significant))  
write.table(dem_frequency, "FrequencyofDEM.txt",sep="\t", quote=F, row.names=F)

threshold_pct <- 50
robust_DEMs <- dem_frequency %>% filter(Frequency_Pct >= threshold_pct)

# Supplementary Figure 14c
pdf("9_BarPlot_FrequentDEMs.pdf")
top_n <- 50
dem_frequency %>%
  slice_max(Frequency_Pct, n = top_n) %>%
  ggplot(aes(x = reorder(metab_id, Frequency_Pct), y = Frequency_Pct)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top DEMs by Significance Frequency Across Permutations",
    x = "Metabolite",
    y = "Percentage of Permutations Where Significant (%)"
  ) +
  theme_minimal()
dev.off()

perm_results_annotated <- perm_results %>%
  left_join(dem_frequency, by = "metab_id")

library(ggplot2)
library(ggrepel) 

perm_results_annotated <- perm_results_annotated %>%
  group_by(metab_id) %>%
  mutate(
    Label = ifelse(row_number() == 1 & Pattern != "NonSig" & Frequency_Pct >= 75, metab_id, NA)
  ) %>%
  ungroup()

CustomLabel <- c("glycerophosphoethanolamine", "glycerophosphorylcholine..GPC.","myo.inositol","N.acetyl.aspartyl.glutamate..NAAG.")

pdf("10_VolcanoPlot_DeduplicatedLabels.pdf",height=8, width=10)
ggplot(perm_results_annotated, aes(x = effect, y = -log10(padj), color = Frequency_Pct)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(
    aes(label = Label), 
    max.overlaps = 18,
    size = 3,
    nudge_x = ifelse(perm_results_annotated$effect > 0, 0.5, -0.5), # Push labels left/right based on effect
    nudge_y = 2, # Push labels slightly upward
    direction = "y", # Prioritize vertical adjustment
    hjust = "outward", # Align labels outward from points
    segment.color = "gray", # Color of the connecting line
    segment.size = 0.5, # Thickness of the connecting line
    box.padding = 0.5, # Padding around labels
    force = 1 # Adjust repulsion force
  ) +
  scale_color_gradient(low = "gray", high = "red",
                       name = "Significance\nFrequency (%)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Permutation Results: Effect Size vs. Significance",
    subtitle = "Color = Frequency of significance across permutations",
    x = "Effect Size",
    y = "-log10(Adjusted P-value)"
  ) +  
  theme_bw() +
  coord_cartesian(clip = "off") # Allow labels to extend beyond plot area
dev.off()

pdf("10a_VolcanoPlot.pdf",height=8, width=10)
ggplot(perm_results_annotated, aes(x = effect, y = -log10(padj), color = Frequency_Pct)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(
    aes(label = NA), 
    max.overlaps = 18,
    size = 3,
    nudge_x = ifelse(perm_results_annotated$effect > 0, 0.5, -0.5), # Push labels left/right based on effect
    nudge_y = 2, # Push labels slightly upward
    direction = "y", # Prioritize vertical adjustment
    hjust = "outward", # Align labels outward from points
    segment.color = "gray", # Color of the connecting line
    segment.size = 0.5, # Thickness of the connecting line
    box.padding = 0.5, # Padding around labels
    force = 1 # Adjust repulsion force
  ) +
  scale_color_gradient(low = "gray", high = "red",
                       name = "Significance\nFrequency (%)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Permutation Results: Effect Size vs. Significance",
    subtitle = "Color = Frequency of significance across permutations",
    x = "Effect Size",
    y = "-log10(Adjusted P-value)"
  ) +  
  theme_bw() +
  coord_cartesian(clip = "off") # Allow labels to extend beyond plot area
dev.off()

library(pheatmap)
library(tibble)
library(tidyverse)

top_50_metabs <- perm_results %>%
  filter(Pattern %in% c("Up", "Down")) %>%  
  count(metab_id, name = "Frequency") %>%   
  arrange(desc(Frequency)) %>%              
  slice_head(n = 50) %>%                    
  pull(metab_id)     

heatmap_data <- perm_results %>%
  mutate(Significant = ifelse(Pattern %in% c("Up", "Down"), 1, 0)) %>%
  select(metab_id, Permutation, Significant) %>%
  pivot_wider(names_from = Permutation, values_from = Significant, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("metab_id")

pdf("11_Heatmap_DEMFrequency.pdf", height=4,width=8)
pheatmap(
  heatmap_data[rownames(heatmap_data) %in% top_50_metabs, ],
  color = c("white", "red"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,          # Show metabolite IDs
  fontsize_row = 7,              # Adjust label size
  main = "Top 50 Most Frequently Significant Metabolites"
)
dev.off()


#### Box Plots ####
ClinicalMeta$age = ifelse(ClinicalMeta$age_at_visit_max=="90+", 90, ClinicalMeta$age_at_visit_max)
ClinicalMeta$Sex = ifelse(ClinicalMeta$msex==0, "Female","Male")
all(rownames(ClinicalMeta)==colnames(Metab_Raw))
Clinical_Metab = cbind(ClinicalMeta, t(Metab_Raw))
Clinical_Metab[1:6,1:20]
Clinical_Metab_cogdx = Clinical_Metab[Clinical_Metab$cogdx %in% c(1,2,4,5),]

NAAG=ggplot(Clinical_Metab_cogdx, aes(x=as.character(cogdx), y=`N.acetyl.aspartyl.glutamate..NAAG.`,color=as.character(cogdx)))+
  geom_boxplot()+labs(title="",x="Cogdx score",y="Expression")+ #ylim(c(7,9))+
  #ggpubr::stat_compare_means()+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  theme(legend.position="none")+ geom_jitter(shape=1, position=position_jitter(0.2))+
  labs(title="NAAG")+
  theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),color="black"),axis.text.y = element_text(size=rel(1.0)))+
  facet_grid(~Sex)

pdf("12_NAAG.pdf",width=6,height=3)
print(NAAG)
dev.off()

Aspartate=ggplot(Clinical_Metab_cogdx, aes(x=as.character(cogdx), y=`aspartate`,color=as.character(cogdx)))+
  geom_boxplot()+labs(title="",x="Cogdx score",y="Expression")+ #ylim(c(7,9))+
  #ggpubr::stat_compare_means()+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  theme(legend.position="none")+ geom_jitter(shape=1, position=position_jitter(0.2))+
  labs(title="Aspartate")+
  theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),color="black"),axis.text.y = element_text(size=rel(1.0)))+
  facet_grid(~Sex)

pdf("13_Aspartate.pdf",width=6,height=3)
print(Aspartate)
dev.off()

Glutamate=ggplot(Clinical_Metab_cogdx, aes(x=as.character(cogdx), y=`glutamate`,color=as.character(cogdx)))+
  geom_boxplot()+labs(title="",x="Cogdx score",y="Expression")+ #ylim(c(7,9))+
  #ggpubr::stat_compare_means()+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  theme(legend.position="none")+ geom_jitter(shape=1, position=position_jitter(0.2))+
  labs(title="Aspartate")+
  theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),color="black"),axis.text.y = element_text(size=rel(1.0)))+
  facet_grid(~Sex)

pdf("14_Glutamate.pdf",width=6,height=3)
print(Glutamate)
dev.off()



