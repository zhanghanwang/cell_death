rm(list=ls())
gc()
options(stringsAsFactors = F)

library(survival)
library(survminer)
library(edgeR)
library(dplyr)
library(survivalROC)
library(timeROC)
library(glmnet)
library(pheatmap)
library(gplots)
library(colorRamps)
library(RColorBrewer)
library(tidyverse)
library(future.apply)
library(ggrisk)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(ggthemes)
plan(multisession)
options(future.globals.maxSize= 891289600)
setwd("E:/proj/LGG/death_pathway3/")


########################################################
## download                ############################
## TCGA SNV                #####################################
proj="TCGA-LGG"
proj1=tolower(gsub("-","_",proj))
##您提供的代码似乎与使用 TCGAbiolinks 软件包查询和准备屏蔽体细胞突变数据有关。它主要是查询体细胞突变数据、下载数据，然后处理并保存相关信息。
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(
  project = proj, 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open")
GDCdownload(query)
GDCprepare(query, save = T,save.filename = paste0(proj1,"_SNP.Rdata"))
load(paste0(proj1,"_SNP.Rdata"))
mut=assayData(data)
mut$t_vaf=mut$t_alt_count/mut$t_depth
save(mut,file=paste0(proj1,"_mut.RData"))


##这部分代码包括使用 maftools 软件包以及 dplyr 和 tidyr 进行计算，并创建肿瘤突变负荷 (TMB) 的汇总得分。
library(maftools)
library(dplyr)
library(tidyr)
maf= read.maf(maf= mut)
##您可以使用 maftools 软件包中的 read.maf() 函数，将突变数据从 mut 对象读取到名为 maf 的 MAF 对象中。
math_score=math.score(maf= maf, vafCutOff= 0.05)
##使用 maftools 软件包中的 math.score() 函数计算每个样本的数学分数。vafCutOff 参数设置为 0.05，以指定变异等位基因频率截止值

TMB_score <- mut %>% group_by(Tumor_Sample_Barcode) %>% summarise(TCGA_sum = n())##使用 dplyr 中的 group_by()，按肿瘤样本条形码对突变数据进行分组。汇总每个样本的突变数 (TCGA_sum)。

TMB_score$TMB=round(TMB_score$TCGA_sum/38,3)##用突变总数除以 38（归一化因子）计算 TMB 分数。结果存储在 TMB_score 数据框中。

maf_tmb_score = merge(math_score,TMB_score,by="Tumor_Sample_Barcode")
##使用 maftools 软件包中的 mergeL() 根据共同列 Tumor_Sample_Barcode 合并 math_score 和 TMB_score 数据帧。

maf_tmb_score=as.data.frame(maf_tmb_score)
##使用 as.data.frame() 将合并后的数据帧 maf_tmb_score 转换为常规数据帧。
table(substr(maf_tmb_score$Tumor_Sample_Barcode,14,16))##检查 Tumor_Sample_Barcode 第 14 至 16 位字符的分布情况。
save(maf_tmb_score,file=paste0(proj1,"_tmb.RData"))
##这些步骤包括计算突变分数和肿瘤突变负荷、合并结果，以及准备数据供进一步分析。

## TCGA CNV###############

TCGA.LGG.gistic <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/TCGA-LGG.gistic.tsv")

`gencode.v22.annotation.gene.(1)` <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/gencode.v22.annotation.gene (1).probeMap")

# 查看每列的数据类型
sapply(TCGA.LGG.gistic, class)
gtf23 = `gencode.v22.annotation.gene.(1)`

# 删除包含缺失值的行
TCGA.LGG.gistic <- na.omit(TCGA.LGG.gistic)

TCGA.LGG.gistic =as.data.frame(lapply(TCGA.LGG.gistic,as.numeric))
TCGA.LGG.gistic$gene_name = gtf23[match(rownames(TCGA.LGG.gistic),rownames(gtf23)),"gene"] ## ID转换，将ENSG00000223972.5 转换为 DDX11L1
TCGA.LGG.gistic = limma::avereps(TCGA.LGG.gistic[,!colnames(TCGA.LGG.gistic) %in% "gene_name"],ID = TCGA.LGG.gistic$gene_name) ## 去除重复gene_name，重复的取平均值
TCGA.LGG.gistic = as.data.frame(TCGA.LGG.gistic)
save(TCGA.LGG.gistic,file="TCGA_LGG_gistic.RData")


## TCGA micRNA             #####################################
query <- GDCquery(project = proj,data.category = "Transcriptome Profiling",data.type = "miRNA Expression Quantification")
GDCdownload(query)
GDCprepare(query, save = T,save.filename = paste0(proj1,"_miRNA.Rdat"))

load(paste0(proj1,"_miRNA.Rdat"))
mic_count=data[c("miRNA_ID",colnames(data)[grepl("read_count",colnames(data))])]
rownames(mic_count)=mic_count$miRNA_ID
mic_count$miRNA_ID=NULL
colnames(mic_count)=gsub("read_count_","",colnames(mic_count))
mic_rpm=data[c("miRNA_ID",colnames(data)[grepl("reads_per_million",colnames(data))])]
rownames(mic_rpm)=mic_rpm$miRNA_ID
mic_rpm$miRNA_ID=NULL
colnames(mic_rpm)=gsub("reads_per_million_miRNA_mapped_","",colnames(mic_rpm))
save(mic_count,file = paste0(proj1,"_mic_count.RData"))
save(mic_rpm,file = paste0(proj1,"_mic_rpm.RData"))

## TCGA pheno              ########################################################
phen <- read.delim("G:/fangshe/diverse cell-death patterns/TCGA-LGG.GDC_phenotype.tsv")
colnames(phen)

# 定义列名
column_names <- colnames(phen)

# 检查列名中是否包含"_cm"
matches <- grepl("codeletion", column_names)

# 打印出包含"_cm"的列名
print(column_names[matches])

phen$ldh1_mutation=ifelse(phen$ldh1_mutation_tested %in% "YES" & phen$ldh1_mutation_found %in% "YES",1,ifelse(phen$ldh1_mutation_tested %in% "YES" & phen$ldh1_mutation_found %in% "NO",0,NA))##基于条件逻辑，将 ldh1_mutation 列中的值进行处理，将满足条件的数据设置为1，不满足条件的设置为0，未知情况为NA。

phen$os=ifelse(phen$vital_status.demographic %in% "Alive",0,ifelse(phen$vital_status.demographic %in% "Dead",1,NA))
phen$days_to_last_follow_up.diagnoses=as.numeric(phen$days_to_last_follow_up.diagnoses)

phen$days_to_death.demographic=as.numeric(phen$days_to_death.demographic)

phen$os_time=ifelse(is.na(phen$days_to_death.demographic),phen$days_to_last_follow_up.diagnoses,phen$days_to_death.demographic)

phen$rfs=ifelse(phen$new_tumor_event_after_initial_treatment %in% "NO",0,ifelse(phen$new_tumor_event_after_initial_treatment %in% "YES",1,NA))

phen$days_to_new_tumor_event_after_initial_treatment=as.numeric(phen$days_to_new_tumor_event_after_initial_treatment)

phen$rfs_time=ifelse(phen$rfs %in% 1,phen$days_to_new_tumor_event_after_initial_treatment,ifelse(phen$rfs %in% 0,phen$os_time,NA))

phen=unique(phen[,c(1,88,6,64,54,63,33,133,59,134:137)])
colnames(phen) = c("Sample","Gender","Age","radiation_therapy","Grade","primary_therapy_outcome","followup_treatment_success","ldh1_mutation","person_neoplasm_cancer_status","os","os_time","rfs","rfs_time")

unique(phen$new_tumor_event_after_initial_treatment)
which(colnames(phen) %in% "person_neoplasm_cancer_status")

rownames(phen)=phen$Sample
table(phen$radiation_therapy)
phen=phen[phen$radiation_therapy %in% "NO",]

summary(phen$rfs_time)
phen$rfs_time=phen$rfs_time/365
phen$os_time=phen$os_time/365
phen_NO = phen
save(phen,file = "tcga_lgg_phen.RData")
surv=phen[,c(10:13)]
save(surv,file = "tcga_lgg_surv.RData")

## TCGA mRNA             #####################################
rm(list=ls())
gc()

load("gencode_annotation_gtf23.RData")
load("tcga_gtex_phenotype.RData")
load("tcga_lgg_phen.RData")

tcga_gtex_phenotype = tcga_gtex_phenotype[substr(tcga_gtex_phenotype$sample,1,15) %in% substr(phen$Sample,1,15) |                                    tcga_gtex_phenotype$X_study %in% "GTEX",]

samp = tcga_gtex_phenotype[tcga_gtex_phenotype$tcga_type %in% "LGG","sample"] %>% as.character()

#可以过滤表型数据,只保留TCGA_LGG 队列中的相关样本（"LGG "类型）。从筛选出的表型数据中提取与 "LGG "类型相匹配的样本名称

library(data.table)
tpm_sva=fread("E:/database/TCGA/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz",data.table = F)
tpm_sva[1:3,1:3]
rownames(tpm_sva) = tpm_sva$sample
tpm_sva$sample = NULL

table(samp %in% colnames(tpm_sva))
tpm_sva = tpm_sva[,samp]
tpm_sva = 2**tpm_sva-2**min(tpm_sva)
tpm_sva = log(tpm_sva+1)

#使用 data.table 软件包中的 fread 函数从 gzip 文件中读取基因表达数据。将行名设置为 "sample "列，并删除 "sample "列。检查样本名称是否与 tpm_sva 中的列名匹配。过滤 tpm_sva 中的列，只保留 samp 向量中的样本。对 tpm_sva 进行数据转换，将 TPM 值转换为 log2(TPM + 1) 值。检查数据中的基因 ID是否与注释数据中的基因 ID 一致。

table(rownames(tpm_sva) %in% gtf23$gene_id)
tpm_sva$gene_id=rownames(tpm_sva)
tpm_sva=merge(gtf23[,-2],tpm_sva,by="gene_id")
tpm_sva=limma::avereps(tpm_sva[,3:ncol(tpm_sva)],ID=tpm_sva$gene_name)
save(tpm_sva,file="tcga_gtex_lgg_tpm_all.RData")
boxplot(tpm_sva[,1:20])

tpm_sva=tpm_sva[,substr(colnames(tpm_sva),1,4)=="TCGA" & substr(colnames(tpm_sva),14,15)<10]
tpm_sva=as.data.frame(tpm_sva)
save(tpm_sva,file="tcga_gtex_lgg_tpm.RData")##TCGA和GTEX的表达谱



## CGGC tcga mrna tpm_sva  #######################################
rm(list=ls())

CGGA.mRNAseq_325_clinical <- read.delim("G:/fangshe/cuproptosis/mRNAseq_325/CGGA.mRNAseq_325_clinical.20200506.txt", row.names=1)
CGGA.mRNAseq_693_clinical <- read.delim("G:/fangshe/diverse cell-death patterns/CGGA.mRNAseq__Counts-genes/CGGA.mRNAseq_693_clinical.20200506.txt", row.names=1)

pheno=rbind(CGGA.mRNAseq_325_clinical,CGGA.mRNAseq_693_clinical)
table(pheno$IDH_mutation_status)
table(pheno$Radio_status..treated.1.un.treated.0.)
pheno=pheno[pheno$Radio_status..treated.1.un.treated.0. %in% 1,]
class(pheno$IDH_mutation_status)
pheno$IDH_mutation_status=factor(pheno$IDH_mutation_status,levels = c("Wildtype","Mutant")) %>% as.numeric(.)-1 ##对 IDH 突变状态进行了转换，将 "Wildtype" 编码为0，"Mutant" 编码为1
table(pheno$Grade)
pheno$Grade=factor(pheno$Grade,levels = c("WHO II","WHO III"),labels = c("G2","G3"))

pheno=pheno[pheno$Grade %in% c("G2","G3"),c(1,5,6,7,8,11,4)] ##只保留 "Grade" 为 "G2" 和 "G3" 的数据

colnames(pheno)=c("Sample","Gender","Age","os_time","os","ldh1_mutation","Grade")
pheno$os_time=pheno$os_time/365
save(pheno,file="phen_cgga_lgg.RData")

load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_phen.RData")
colnames(phen)

phen$Sample=substr(phen$Sample,1,15)

phen_sva=rbind(pheno,phen[,c("Sample","Gender","Age","os_time","os","ldh1_mutation","Grade")])
# phen_sva$Sample=colnames(tpm_all)[match(substr(phen_sva$Sample,1,12),substr(colnames(tpm_all),1,12))]
save(phen_sva,file="phen_sva_lgg.RData")

CGGA.mRNAseq_325 <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/CGGA.mRNAseq_325.RSEM-genes.20200506.txt/CGGA.mRNAseq_325.RSEM-genes.20200506.txt")
CGGA.mRNAseq_693 <- read.delim("G:/fangshe/diverse cell-death patterns/CGGA.mRNAseq_693.RSEM-genes.20200506.txt/CGGA.mRNAseq_693.RSEM-genes.20200506.txt")

fpkm = merge(CGGA.mRNAseq_325,CGGA.mRNAseq_693,by="Gene_Name")
rownames(fpkm) = fpkm$Gene_Name
fpkm$Gene_Name = NULL
fpkmToTpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
tpm1 = apply(fpkm,2,fpkmToTpm)
tpm1 = as.data.frame(tpm1)
min(tpm1)
boxplot(tpm1[,1:10])
tpm1=log(tpm1+1)
table(colnames(tpm1) %in% pheno$Sample)
tpm1=tpm1[,colnames(tpm1) %in% pheno$Sample]

load("tcga_gtex_lgg_tpm.RData")##这个是什么文件
boxplot(tpm_sva[,1:10])
min(tpm_sva)

com_gene = intersect(rownames(tpm_sva),rownames(tpm1))
tpm_all = cbind(tpm_sva[com_gene,],tpm1[com_gene,])

min(tpm_all)
anyNA(tpm_all)
coldata=rbind(data.frame(sample=colnames(tpm_sva),group="TCGA"),
              data.frame(sample=colnames(CGGA.mRNAseq_325)[-1],group="cgga325"),
              data.frame(sample=colnames(CGGA.mRNAseq_693)[-1],group="cgga693"))
table(coldata$sample %in% colnames(tpm_all))
table(coldata$sample %in% phen_sva$Sample)
coldata=coldata[coldata$sample %in% colnames(tpm_all),]
coldata=coldata[coldata$sample %in% phen_sva$Sample,]
save(coldata,file="tpm_sva_lgg_coldata.RData")
library(sva)
tpm_sva=ComBat(dat = tpm_all[,coldata$sample],batch = coldata$group)
tpm_sva=as.data.frame(tpm_sva)
save(tpm_sva,file="tpm_sva_lgg_tpm.RData")##放疗病人矩阵
## Rembrandt ########################################
rm(list=ls())
gse_pheno <- read.delim("E:/database/CGGA/Rembrandt_mRNA_array_475_clinical.txt")
unique(gse_pheno$Grade)
gse_pheno=gse_pheno[gse_pheno$Grade %in% c("WHO I","WHO II","WHO III") & !gse_pheno$Histology %in% "GBM",]
rownames(gse_pheno)=gse_pheno$gene
gse_pheno$Grade=factor(gse_pheno$Grade,levels = c("WHO II","WHO III"),labels = c("G2","G3"))
colnames(gse_pheno)
colnames(gse_pheno)=c("Sample","Histology","Grade","Gender","Age","os_time","os","X1p19q_Codeletion_status")
summary(gse_pheno$os_time)
gse_pheno$os_time=gse_pheno$os_time/365
save(gse_pheno,file="Rembrandt_lgg_gse_pheno.RData")

gse_exp <- read.delim("E:/database/CGGA/Rembrandt_mRNA_array_475.txt")
gse_exp=limma::avereps(gse_exp[,2:ncol(gse_exp)],ID=gse_exp$gene)
gse_exp=as.data.frame(gse_exp[,colnames(gse_exp) %in% gse_pheno$Sample])
save(gse_exp,file="Rembrandt_lgg_gse_exp.RData")

####################################
rm(list=ls())

load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_phen.RData")
# 获取原始行名
old_row_names <- rownames(tcga_risk_exp)
# 在每个行名前面加上字母"A"
new_row_names <- paste0( old_row_names,"A")
# 将新的行名赋值给数据框
rownames(tcga_risk_exp) <- new_row_names
table(rownames(phen) %in% rownames(tcga_risk_exp))
#phen=phen[(phen$Sample) %in% rownames(tcga_risk_exp),]
phen=cbind(phen[rownames(tcga_risk_exp),],tcga_risk_exp[,c("Risk_Level","Risk_Score")])

a= c("NA","NA.1","NA.2","NA.3","NA.4","NA.5","NA.6","NA.7")
phen <- phen[!(rownames(phen) %in% a),]

unique(phen$Gender)
phen$Gender=toupper(gsub(" ","",phen$Gender))
phen$Gender=factor(phen$Gender,levels = c("FEMALE","MALE"))
save(phen,file = "phen_lgg.RData")

##########################################################
################# 数据处理 ###############################
##########################################################
## gene list   #######################################
rm(list=ls())
load("E:/database/database/genelist/death_gene1.RData")
colnames(death_gene)=c("Gene","Type" )
death_gene[death_gene$Type %in% "Necrosis","Type"]="Necroptosis"
table(death_gene$Type)

load("E:/database/database/GTF_annotation/gencode_annotation_gtf23.RData")
gene_list <- read.delim("E:/proj/LGG/death_pathway/gene_list.txt")
table(gene_list$Type)
gene_list=rbind(gene_list[,1:2],death_gene)
gene_list=unique(gene_list)

table(gene_list$Gene %in% gtf23$gene_name)
unique(gene_list[!gene_list$Gene %in% gtf23$gene_name,])
table(gene_list$Type )
gene_list=gene_list[gene_list$Gene %in% gtf23$gene_name,]
gene_list$Type=gsub(" ","",gene_list$Type)
save(gene_list,file="gene_list.RData")
## PCD基因的放疗表达谱#####

###整理CGGA的放疗临床信息
rm(list=ls())

CGGA.mRNAseq_325_clinical <- read.delim("G:/fangshe/cuproptosis/mRNAseq_325/CGGA.mRNAseq_325_clinical.20200506.txt", row.names=1)
CGGA.mRNAseq_693_clinical <- read.delim("G:/fangshe/diverse cell-death patterns/CGGA.mRNAseq__Counts-genes/CGGA.mRNAseq_693_clinical.20200506.txt", row.names=1)

pheno_no=rbind(CGGA.mRNAseq_325_clinical,CGGA.mRNAseq_693_clinical)
table(pheno_no$IDH_mutation_status)
table(pheno_no$Radio_status..treated.1.un.treated.0.)
pheno_no=pheno_no[pheno_no$Radio_status..treated.1.un.treated.0. %in% 0,]
class(pheno_no$IDH_mutation_status)
pheno_no$IDH_mutation_status=factor(pheno_no$IDH_mutation_status,levels = c("Wildtype","Mutant")) %>% as.numeric(.)-1 ##对 IDH 突变状态进行了转换，将 "Wildtype" 编码为0，"Mutant" 编码为1

table(pheno_no$Grade)
pheno_no$Grade=factor(pheno_no$Grade,levels = c("WHO II","WHO III"),labels = c("G2","G3"))

pheno_no=pheno_no[pheno_no$Grade %in% c("G2","G3"),c(1,5,6,7,8,11,4)] ##只保留 "Grade" 为 "G2" 和 "G3" 的数据

colnames(pheno_no)=c("Sample","Gender","Age","os_time","os","ldh1_mutation","Grade")
pheno_no$os_time=pheno_no$os_time/365

save(pheno_no,file="phenno_cgga_lgg.RData")

###整理CGGA的表达矩阵###
CGGA.mRNAseq_325 <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/CGGA.mRNAseq_325.RSEM-genes.20200506.txt/CGGA.mRNAseq_325.RSEM-genes.20200506.txt")
CGGA.mRNAseq_693 <- read.delim("G:/fangshe/diverse cell-death patterns/CGGA.mRNAseq_693.RSEM-genes.20200506.txt/CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
load("phenno_cgga_lgg.RData")
load("phen_cgga_lgg.RData")
load("gene_list.RData")
fpkm = merge(CGGA.mRNAseq_325,CGGA.mRNAseq_693,by="Gene_Name")
rownames(fpkm) = fpkm$Gene_Name
fpkm$Gene_Name = NULL
fpkmToTpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
tpm = apply(fpkm,2,fpkmToTpm)
tpm = as.data.frame(tpm)
min(tpm1)
boxplot(tpm[,1:10])
table(colnames(tpm) %in% rownames(pheno_no))
tpm1=tpm[,colnames(tpm) %in% rownames(pheno_no)]
tpm2=tpm[,colnames(tpm) %in% rownames(pheno)]
tpm1 = tpm1[rownames(tpm1) %in% gene_list$Gene,]
tpm2 = tpm2[rownames(tpm2) %in% gene_list$Gene,]
# 输出文件
write.table(tpm1, "CGGA_matrix_no.txt",sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(tpm2, "CGGA_matrix_yes.txt",sep = "\t",
            row.names = TRUE, col.names = NA)

## TIDE                    ############################################
load("G:/fangshe/cell_death_qj2/death_pathway3/gencode_annotation_gtf23.RData")
table(gtf23$gene_type)
mrna=gtf23[gtf23$gene_type %in% "protein_coding",]
load("tpm_sva_lgg_tpm.RData")
tpm_sva=as.data.frame(tpm_sva[rownames(tpm_sva) %in% mrna$gene_name,])
tpm_sva=scale(tpm_sva,scale = T,center = T)
boxplot(tpm_sva[,1:10])
tpm_sva=as.data.frame(tpm_sva)
tpm_sva1=tpm_sva[,sample(736,350)]
tpm_sva2=tpm_sva[,!colnames(tpm_sva) %in% colnames(tpm_sva1)]
write.table(tpm_sva1,file = "tcga_lgg_tide_tpm1.tsv",quote = F,sep = "\t",row.names = T)
write.table(tpm_sva2,file = "tcga_lgg_tide_tpm2.tsv",quote = F,sep = "\t",row.names = T)
lgg_tide1 <- read.csv("C:/Users/WIN 10/Downloads/tcga_lgg_tide_tpm1.csv")
lgg_tide2 <- read.csv("C:/Users/WIN 10/Downloads/tcga_lgg_tide_tpm2.csv")
wilcox.test(lgg_tide2$TIDE,lgg_tide1$TIDE)
lgg_tide=rbind(lgg_tide2,lgg_tide1)
save(lgg_tide,file="lgg_tide.RData")

## ssgsea                  #####################################
rm(list=ls())
proj1="lgg_sva"
load("tpm_sva_lgg_tpm.RData")

library(GSVA)
gene_set<-read.delim("G:/fangshe/cell_death_qj2/death_pathway3/mm3_Metagene_Immunity.txt")
ssgsea_list<- split(as.matrix(gene_set)[,1], gene_set[,2])

gsva_matrix<- gsva(as.matrix(tpm_sva), ssgsea_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix=as.data.frame(t(gsva_matrix))
save(gsva_matrix,file=paste0(proj1,"_ssgsea.RData"))

library(clusterProfiler)
kegg_geneset=read.gmt("G:/fangshe/cell_death_qj2/death_pathway3/c2.cp.kegg.v2022.1.Hs.symbols.gmt") 
kegg_geneset <- split(as.matrix(kegg_geneset)[,2], kegg_geneset[,1])
gsva_matrix<- gsva(as.matrix(tpm_sva), kegg_geneset,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
save(gsva_matrix,file=paste0(proj1,"_ssgsea_kegg.RData"))


load("G:/fangshe/cell_death_qj2/death_pathway3/gene_list.RData")
gene_list <- split(as.matrix(gene_list)[,2], gene_list[,1])
gsva_matrix <- gsva(as.matrix(tpm_sva), gene_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
save(gsva_matrix,file=paste0(proj1,"_ssgsea_death_pathway.RData"))

## CIBERSORT               ############################################
source("D:/R_code/CIBERSORT.R")
write.table(tpm_sva,"cib_data.tsv",sep = "\t",row.names = T ,col.names = T,quote = F)
cib_result=CIBERSORT("./cib_data.tsv","D:/R_code/cibersort/LM22.txt",perm = 500,absolute = T  ) %>% as.data.frame()
save(cib_result,file=paste0(proj1,"_cib_raw_result.RData"))

## xcell                   #####################################
source("D:/R_code/xCell-master/xCell-master/R/xCell.R")
load("D:/R_code/xCell-master/xCell-master/data/xCell.data.rda")
xcell_result=xCellAnalysis(expr = tpm_sva) %>% as.data.frame() 
save(xcell_result,file=paste0(proj1,"_xcell.RData"))

## EPIC                    ################################
library(EPIC)
boxplot(tpm_sva[,1:10])
min(2**tpm_sva)
epic_score <- EPIC(bulk = 2**tpm_sva-min(2**tpm_sva))
save(epic_score,file=paste0(proj1,"_epic_score.RData"))

## MCPcounter              ################################
library(MCPcounter)
probesets <- read.table("E:/database/pan_cancer/data/mcp_count_data/probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character")
genes <- read.table("E:/database/pan_cancer/data/mcp_count_data/genes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
MCPcounterScore <- MCPcounter.estimate(tpm_sva, featuresType = "HUGO_symbols",probesets = probesets,genes = genes)
save(MCPcounterScore,file=paste0(proj1,"_MCPcounterScore.RData"))

## quantiseq               ################################
## remotes::install_github("icbi-lab/immunedeconv")
library(immunedeconv)
quantiseq_score <- immunedeconv::deconvolute(2**tpm_sva-min(2**tpm_sva), "quantiseq")
save(quantiseq_score,file=paste0(proj1,"_quantiseq_score.RData"))

library(immunedeconv)
timer_score <- immunedeconv::deconvolute(2**tpm_sva-min(2**tpm_sva), "timer",indications=rep("lgg",ncol(tpm_sva)))
timer_score=as.data.frame(timer_score)
rownames(timer_score)=timer_score$cell_type
timer_score$cell_type=NULL
save(timer_score,file=paste0(proj1,"_timer_score.RData"))

## IPS                     #####################################
library(ggplot2)
library(grid)
library(gridExtra)
gene_expression<- as.data.frame(tpm_sva)
sample_names<- colnames(gene_expression)

ipsmap<- function (x) {
  if(x<=0){ips<-0}else{if(x>=3){ips<-10}else{ips<-round(x*10/3, digits=0)}}
  return(ips)
}
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function(x){
  za<-NULL
  if (x>=3){za=1000}else{if(x<=-3){za=1}else{za=round(166.5*x+500.5,digits=0)}}
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if(x>=2){za2=1000}else{if(x<=-2){za2=1}else{za2=round(249.75*x+500.5,digits=0)}}
  return(my_palette2[za2])
}
IPSG<-read.table("E:/database/pan_cancer/data/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes<-as.vector(unique(IPSG$NAME))
IPS<-NULL
MHC<-NULL
CP<-NULL
EC<-NULL
SC<-NULL
AZ<-NULL
GVEC<-row.names(gene_expression)
VEC<-as.vector(IPSG$GENE)
ind<-which(is.na(match(VEC,GVEC)))
MISSING_GENES<-VEC[ind]
dat<-IPSG[ind,]
if(length(MISSING_GENES)>0){cat("differently named or missing genes: ",MISSING_GENES,"\n")}
for(x in 1:length(ind)){print(IPSG[ind,])}

for (i in 1:length(sample_names)) {	
  GE<-gene_expression[[i]]
  mGE<-mean(GE)
  sGE<-sd(GE)
  Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1<-IPSG$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k<-k+1
  }
  WG<-MIG*WEIGHT
  MHC[i]<-mean(WG[1:10],na.rm=T)
  CP[i]<-mean(WG[11:20],na.rm=T)
  EC[i]<-mean(WG[21:24],na.rm=T)
  SC[i]<-mean(WG[25:26],na.rm=T)
  AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i],na.rm=T)
  IPS[i]<-ipsmap(AZ[i])
}
ips_score<-data.frame(sample=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
save(ips_score,file=paste0(proj1,"_ips_score.RData"))

## oncoPredict             #####################################
##您提供的代码似乎是 "oncoPredict "分析流程的一部分，这是一个用于药物反应预测和生物标记发现的软件包。下面是代码的详细说明
proj1 = 'TCGA'
library(oncoPredict)
dir="G:/fangshe/cell_death_qj2/death_pathway3/osfstorage-archive/DataFiles/Training_Data"
GDSC2_Expr = readRDS(file=file.path(dir,"GDSC2_Expr_RMA_Log_Transformed.rds"))##这是个啥文件
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))##这是个啥文件
load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
 
tpm_TCGA =tpm_sva[,substr(colnames(tpm_sva),1,4) %in% "TCGA"]
##您可以从 dir变量指定的目录中加载基因表达数据和药物反应数据。GDSC2_Expr 从 "GDSC2_Expr_RMA_Log_Transformed.rds "文件加载。GDSC2_Res 从 "GDSC2_Res.rds "文件中加载。该数据被指数化，可能是为了反向进行对数变换。
GDSC2_Res <- exp(GDSC2_Res)

calcPhenotype(
  trainingExprData = GDSC2_Expr,
  trainingPtype = GDSC2_Res,
  testExprData = as.matrix(tpm_TCGA),
  batchCorrect = 'eb', 
  powerTransformPhenotype = TRUE,       removeLowVaryingGenes = 0.2,
  minNumSamples = 10, 
  printOutput = TRUE, 
  removeLowVaringGenesFrom = 'rawData' )
##您可以使用 "oncoPredict "软件包中的 calcPhenotype 函数，根据训练基因表达数据 (GDSC2_Expr) 和训练药物反应数据 (GDSC2_Res) 计算药物敏感性预测值。
DrugPredictions_GDSC2 <- data.table::fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
save(DrugPredictions_GDSC2,file=paste0(proj1,"_oncoPredict_gdsc2.RData"))
##该函数需要几个参数来控制计算过程，如批次校正、表型的幂次转换、最小样本数等。计算结果保存在 "calcPhenotype_Output "目录下名为 "DrugPredictions.csv "的 CSV 文件中。


## com_gene                ######################################
load("tpm_sva_lgg_tpm.RData")
load("Rembrandt_lgg_gse_exp.RData")
com_gene = intersect(rownames(tpm_sva),rownames(gse_exp))
save(com_gene,file = "com_gene.RData")

load("tpm_sva_lgg_tpm.RData")
load("phen_sva_lgg.RData")

exp=as.data.frame(t(tpm_sva[com_gene,colnames(tpm_sva) %in% phen_sva$Sample]))
exp$os_time=phen_sva[match(rownames(exp),phen_sva$Sample),"os_time"]
exp$os=phen_sva[match(rownames(exp),phen_sva$Sample),"os"]
exp=exp[exp$os_time>0,]
exp=na.omit(exp)
save(exp,file = "train_data.RData")

load("Rembrandt_lgg_gse_exp.RData")
load("Rembrandt_lgg_gse_pheno.RData")
table(gse_pheno$Sample %in% colnames(gse_exp)) #检查 gse_pheno 中的样本名是否与 gse_exp 中的列名一致。
gse_exp=as.data.frame(t(gse_exp[com_gene,]))
gse_exp$os_time=gse_pheno[rownames(gse_exp),6]
gse_exp$os=gse_pheno[rownames(gse_exp),7]
gse_exp=gse_exp[gse_exp$os_time>0,]
gse_exp=na.omit(gse_exp) #从gse_pheno 数据中提取相关生存信息（os_time 和os）。过滤掉 os_time 小于或等于 0 的样本。删除缺失值（na.omit）的行。
save(gse_exp,file = "test_data.RData")

## km_cox                  ####################################
rm(list=ls())
load("com_gene.RData")

load("train_data.RData")
km_cox=do.call(rbind,future_lapply(com_gene,function(x){
  exp$group=ifelse(exp[,x]>median(exp[,x]),"H","L")
  if(length(unique(exp$group))==2){
    surv_diff <- survdiff(Surv(os_time, os) ~ group, data = exp)
    km_p=1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
    cox_res=summary(coxph(Surv(os_time, os)~exp[,x],data = exp))
    data.frame(gene=x,km_p=km_p,cox_p=cox_res[["coefficients"]][5],
               cox_hr=cox_res[["conf.int"]][1],low_CI=cox_res[["conf.int"]][3],up_CI=cox_res[["conf.int"]][4])
  }
}))
km_cox$gene=as.character(km_cox$gene)
save(km_cox,file = "km_cox.RData")

## tumor_normal_de         #################################################
rm(list=ls())
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_gtex_lgg_tpm_all.RData")
library(limma)
condition=data.frame(sample=colnames(tpm_sva),group=ifelse(substr(colnames(tpm_sva),1,4)=="TCGA" & substr(colnames(tpm_sva),14,15)<10,"T","N"))
table(condition$group)
mm=model.matrix(~0+as.character(condition$group))
colnames(mm)=c("Nor","Tum")
res=lmFit(tpm_sva[,condition$sample],mm) %>% contrasts.fit(., makeContrasts(Tum-Nor, levels = mm)) %>% eBayes(.,trend = T) %>% topTable(., number =Inf) %>% .[order(.$adj.P.Val),] %>% na.omit()
table(res$adj.P.Val<0.01 & abs(res$logFC)>1)
save(res,file="tcga_gtex_lgg_norm_tumor_de.RData")



load("multi_cox.Rdata")
multi_cox_sum=summary(multi_cox)
multi_cox_coefficients=as.data.frame(na.omit(multi_cox_sum$coefficients))
rownames(multi_cox_coefficients)=gsub("\`","",rownames(multi_cox_coefficients))
multi_cox_gene=rownames(multi_cox_coefficients)

dat=as.data.frame(t(tpm_sva[multi_cox_gene,]))
dat$Group=ifelse(substr(colnames(tpm_sva),1,4)=="TCGA" & substr(colnames(tpm_sva),14,15)<10,"Tumor","Normal")
table(dat$Group)
dat$samole=rownames(dat)
dat1=reshape2::melt(dat, c("Group","samole") )

pdf("lgg_norm_tumor_de_boxplot.pdf",5,5)
for (i in multi_cox_gene) {
  p = ggplot(dat1[dat1$variable %in% i, ], aes(Group, value)) +
    geom_boxplot(aes(fill = Group)) + 
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black') +
    stat_compare_means(aes(group = Group), hide.ns = TRUE) +
    scale_fill_manual(values = c("#88c4e8", "#c74546")) +
    labs(title = "", x = "", y = paste0(i, " Expression")) +
    theme(text = element_text(size = 16))  # 调整字体大小为 14（可以根据需要修改大小）
  
  print(p)
}
dev.off()
##########################################################################
##################### 主线 #########################################
##########################################################################
## RR_RS_de                #################################################
rm(list=ls())
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_gtex_lgg_tpm.RData")
tpm_sva=tpm_sva[rowSums(tpm_sva>0.5)>ncol(tpm_sva)*0.2,] #过滤 tpm_sva 矩阵中的行，只保留 20% 以上的值大于 0.5 的行。
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_phen.RData")
table(phen$radiation_therapy)
unique(phen$primary_therapy_outcome)
table(phen$followup_treatment_success) #使用表格函数和 unique 函数检查和分析表型数据中的各列。
phen$group=ifelse(phen$primary_therapy_outcome %in% c("Complete Remission/Response"),"RR",ifelse(phen$primary_therapy_outcome %in% c("Progressive Disease"),"RS",NA))
phen=phen[!is.na(phen$group),] #您可以删除组列中有缺失值的行，并调整样本列以匹配基因表达数据的格式。

phen$Sample=substr(phen$Sample,1,15)
phen=phen[phen$Sample %in% colnames(tpm_sva),]

library(limma)
mm=model.matrix(~0+as.character(phen$group))
colnames(mm)=c("RR","RS")
#res=lmFit(tpm_sva[,phen$Sample],mm) %>% contrasts.fit(., makeContrasts(RS-RR, levels = mm)) %>% eBayes(.,trend = T) %>% topTable(., number =Inf) %>% .[order(.$adj.P.Val),] %>% na.omit()

# Fit linear model
fit <- lmFit(tpm_sva[, phen$Sample], mm)

# Apply contrasts
contrast_fit <- contrasts.fit(fit, makeContrasts(RS - RR, levels = mm))

# Perform empirical Bayes moderation
eb_fit <- eBayes(contrast_fit, trend = TRUE)

# Extract top differentially expressed genes
top_genes <- topTable(eb_fit, number = Inf)

# Order results by adjusted p-value and remove rows with NA
ordered_genes <- top_genes[order(top_genes$adj.P.Val),]
final_result <- ordered_genes[!is.na(ordered_genes$adj.P.Val),]

table(res$P.Value<0.01 & abs(res$logFC)>0.5)
table(res$adj.P.Val<0.05 & abs(res$logFC)>0.5)
save(res,file="tcga_lgg_rr_rs_de.RData")
write.csv(res,file = "tcga_lgg_rr_rs_de.csv")
## model                   #######################################
rm(list=ls())

load("G:/fangshe/cell_death_qj2/death_pathway3/km_cox.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/train_data.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/test_data.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_rr_rs_de.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/gene_list.RData")

p_v=0.05   
km_p=0.01   
cor_v=0.5 
res$logFC=-res$logFC

res=res[rownames(res) %in% gene_list$Gene,]##取交集
de_gene=rownames(res[res$adj.P.Val<p_v & abs(res$logFC)> cor_v,])##筛选P值为0.05logfc大于0.5的
de_gene_filt = km_cox[km_cox$gene %in% de_gene & km_cox$cox_p < 0.05,"gene"]##筛选啥的

tcga_exp =exp[substr(rownames(exp),1,4) %in% "TCGA",]
lasso_x=as.matrix(tcga_exp[,de_gene_filt])
lasso_y=tcga_exp[,c("os_time","os")]
colnames(lasso_y)=c("time","status")
lasso_y=as.matrix(lasso_y)
#默认的统计图，设定label = TRUE可以给曲线加上注释，
lasso <- glmnet(lasso_x,lasso_y, family = "cox", alpha = 1 , nlambda = 1000)
plot(lasso,xvar = "lambda", label = F)

library(glmnet)
set.seed(12345678)
lasso_cv_fit = cv.glmnet(lasso_x, 
                         lasso_y, 
                         family = "cox",
                         type.measure = "auc",
                         maxit = 10000)
plot(lasso_cv_fit)##

####构建多因素cox模型##
#1）获取lasso筛选出的基因
#λ值重新建模,选择lambda.min
lasso_cv_fit$lambda.min

model_lasso_min <- glmnet(x=lasso_x, y=lasso_y,family = "cox",lambda=lasso_cv_fit$lambda.min)
model_lasso_min$beta

coef.min = coef(lasso_cv_fit, s = lasso_cv_fit$lambda.min)
#multiCox <- coxph(Surv(OS.time, OS) ~ ., data =  train_surv_expr_vennlasso)

lasso_gene = coef.min@Dimnames[[1]][coef.min@i+1]##这句代码干什么用的

save(lasso_cv_fit,file = "lasso_cv_fit.RData")
save(lasso,file = "lasso.Rdata")

multi_cox_data = tcga_exp[,c("os_time","os",lasso_gene)]
multi_cox= coxph(Surv(os_time, os) ~ . , data =multi_cox_data)

multi_cox_sum = summary(multi_cox)

multi_cox_coefficients=as.data.frame(na.omit(multi_cox_sum$coefficients))
rownames(multi_cox_coefficients)=gsub("\`","",rownames(multi_cox_coefficients))

# Create a risk score plot using ggrisk
ggrisk(multi_cox,code.highrisk = 'High rtCDI',#高风险标签，默认为 ’High’
       code.lowrisk = 'Low rtCDI', #低风险标签，默认为 ’Low’
       title.A.ylab='Risk Score', #A图 y轴名称
       title.B.ylab='Survival Time(year)', #B图 y轴名称，注意区分year month day
        title.A.legend='Risk Group', #A图图例名称
       title.B.legend='Status', #B图图例名称
        title.C.legend='Expression', #C图图例名称
       relative_heights=c(0.1,0.1,0.01,0.15), #A、B、热图注释和热图C的相对高度
       color.A=c(low='#88c4e8',high='#c74546'),#A图中点的颜色
       color.B=c(code.0='#88c4e8',code.1='#c74546'), #B图中点的颜色
       color.C=c(low='#88c4e8',median='#f6f5ee',high='#c74546'), #C图中热图颜色
       vjust.A.ylab=1, #A图中y轴标签到y坐标轴的距离,默认是1
       vjust.B.ylab=2 #B图中y轴标签到y坐标轴的距离,默认是2
       )

dev.off()

# multi_cox_coefficients=multi_cox_coefficients[multi_cox_coefficients$`Pr(>|z|)`<1,]####筛掉多因素中不显著的

save(multi_cox_coefficients,file = "multi_cox_coefficients.RData")

multi_cox_gene=rownames(multi_cox_coefficients)
tcga_wight=as.numeric(multi_cox_coefficients[,1])
multi_cox_gene  #"MMP9"     "HSPB1"    "RAB34"    "SFRP2"    "DIRAS3"   "RRM2"    "NOG"      "ATP6V1G2"

tcga_wight  #  -0.3825510  0.1216995  0.1365023 -0.2482744  0.3181137  0.4182687 -0.3167153 -0.2313615
save(multi_cox_coefficients,file = "multi_cox_coefficients.Rdata")
save(multi_cox,file = "multi_cox.Rdata")

multi_cox_data = tcga_exp[,c("os_time","os",lasso_gene)]
multi_cox= coxph(Surv(os_time, os) ~ . , data =multi_cox_data)

#######验证#######
risk_exp=exp[,c("os_time","os",as.character(multi_cox_gene))]
tcga_risk_exp =risk_exp[substr(rownames(exp),1,4) %in% "TCGA",]
risk=as.character()
for(i in 1:nrow(tcga_risk_exp)){risk=append(risk,sum(tcga_risk_exp[i,multi_cox_gene]*tcga_wight))}

tcga_risk_exp$Risk_Score=as.numeric(risk)
tcga_risk_exp=na.omit(tcga_risk_exp)
tcga_risk_exp$Risk_Level=ifelse(tcga_risk_exp$Risk_Score> median(tcga_risk_exp$Risk_Score),"H","L")
#res.cut <- surv_cutpoint(risk_exp, time = "os_time", event = "os", variables = c("Risk_Score"))
#risk_exp$Risk_Level=ifelse(risk_exp$Risk_Score> res.cut$cutpoint[[1]],"H","L")
table(tcga_risk_exp$Risk_Level)
save(tcga_risk_exp,file = "lgg_tcga_risk_exp.RData")


surv_pvalue(survfit(Surv(os_time, os)~ Risk_Level ,data = tcga_risk_exp))$pval

with(tcga_risk_exp,
     ROC_riskscore <<- timeROC(T=tcga_risk_exp$os_time, 
        delta=tcga_risk_exp$os,
        marker=tcga_risk_exp$Risk_Score,
        cause=1,
        weighting="cox",
        times=c(1,2,3,5),ROC=TRUE))
# t=1       t=2       t=3       t=5 
# 0.9370280 0.9049877 0.9104910 0.7843480 


pdf("TCGA_lgg_ggsurvplot_roc_ggrisk.pdf")
# 创建一个基础的ROC曲线图
p <- plot(ROC_riskscore, time = 1, col = "#c74546", add = F,title = "")
# 添加2年和3年时间点的ROC曲线，并手动指定颜色
p <- p + plot(ROC_riskscore, time = 2, add = T, col = "#88c4e8")+ theme(text=element_text(size=20))
p <- p + plot(ROC_riskscore, time = 3, add = T, col = "#93cc82")+ theme(text=element_text(size=30)) #change font size of all text

# 添加文本标签显示AUC值
p <- p + text(0.5, 0.2, paste("1-Year AUC = ", round(ROC_riskscore$AUC[1], 3)), cex = 1.3)
p <- p + text(0.5, 0.15, paste("2-Year AUC = ", round(ROC_riskscore$AUC[2], 3)), cex = 1.3)
p <- p + text(0.5, 0.1, paste("3-Year AUC = ", round(ROC_riskscore$AUC[3], 3)), cex = 1.3)

##tcga km曲线
pdf("TCGA_lgg_ggsurvplot_km_ggrisk.pdf")
colnames(tcga_risk_exp)[12] = "rtCDI"
ggsurvplot(survfit(Surv(os_time,os) ~ rtCDI,
           data = tcga_risk_exp),
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           legend.title ="rtCDI",
           palette = c("#c74546","#88c4e8"),
           xlab="Time (years)",
           ylab = "Overall survival probability",
           pval.method = T,legend.labs = c("High", "Low"))
dev.off()


###GEO验证
gse_risk_exp=gse_exp[,c(multi_cox_gene,"os_time","os")]
risk=as.character()
for(i in 1:nrow(gse_risk_exp)){risk=append(risk,sum(gse_risk_exp[i,multi_cox_gene]*tcga_wight))}

gse_risk_exp$Risk_Score=as.numeric(risk)
gse_risk_exp$Risk_Level=ifelse(gse_risk_exp$Risk_Score> median(gse_risk_exp$Risk_Score),"H","L")

table(gse_risk_exp$Risk_Level)
save(gse_risk_exp,file = "geo_lgg_risk_exp.RData")
surv_pvalue(survfit(Surv(os_time,os)~Risk_Level ,data = gse_risk_exp))$pval


# Create a risk score plot using ggrisk
multi_cox_data = gse_risk_exp[,c("os_time","os",multi_cox_gene)]
multi_cox= coxph(Surv(os_time, os) ~ . , data =multi_cox_data)
# 创建数据分布对象
ddist <- datadist(multi_cox_data)
pdf("Rembrandt_ggrisk.pdf")
ggrisk(multi_cox,code.highrisk = 'High rtCDI',#高风险标签，默认为 ’High’
       code.lowrisk = 'Low rtCDI', #低风险标签，默认为 ’Low’
       title.A.ylab='Risk Score', #A图 y轴名称
       title.B.ylab='Survival Time(year)', #B图 y轴名称，注意区分year month day
       title.A.legend='Risk Group', #A图图例名称
       title.B.legend='Status', #B图图例名称
       title.C.legend='Expression', #C图图例名称
       relative_heights=c(0.1,0.1,0.01,0.15), #A、B、热图注释和热图C的相对高度
       color.A=c(low='#88c4e8',high='#c74546'),#A图中点的颜色
       color.B=c(code.0='#88c4e8',code.1='#c74546'), #B图中点的颜色
       color.C=c(low='#88c4e8',median='#f6f5ee',high='#c74546'), #C图中热图颜色
       vjust.A.ylab=1, #A图中y轴标签到y坐标轴的距离,默认是1
       vjust.B.ylab=2 #B图中y轴标签到y坐标轴的距离,默认是2
       )
dev.off()

pdf("Rembrandt_lgg_ggsurvplot_roc_ggrisk.pdf")
with(gse_risk_exp,
     ROC_riskscore <<-
timeROC(T=gse_risk_exp$os_time,
        delta=gse_risk_exp$os,
        marker=gse_risk_exp$Risk_Score,
        cause=1,weighting="cox",
        times=c(1,2,3,5),
        ROC=TRUE))
#t=1       t=2       t=3       t=5 
#0.6971649 0.7945122 0.7919501 0.7412095
plot(ROC_riskscore, time = 1, col = "#c74546", add = F,title = "")
plot(ROC_riskscore, time = 2, col = "#88c4e8", add = T)
plot(ROC_riskscore, time = 3, col = "#93cc82", add = T)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)),cex=1.3)
text(0.5,0.15,paste("2-Year AUC = ",round(ROC_riskscore$AUC[2],3)),cex=1.3)
text(0.5,0.1,paste("3-Year AUC = ",round(ROC_riskscore$AUC[3],3)),cex=1.3)
##geo km曲线
ggsurvplot(survfit(Surv(os_time,os) ~ Risk_Level,
           data = gse_risk_exp),
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           legend.title ="rtCDI",
           palette = c("#c74546","#88c4e8"),
           xlab="Time (years)",
           ylab = "Overall survival probability",
           pval.method = T,legend.labs = c("High", "Low"))
dev.off()

##CGGA验证

cgga_risk_exp =exp[substr(rownames(exp),1,4) %in% "CGGA",c(multi_cox_gene,"os_time","os")]
# Create a risk score plot using ggrisk
multi_cox_data = cgga_risk_exp[,c("os_time","os",multi_cox_gene)]
multi_cox= coxph(Surv(os_time, os) ~ . , data =cgga_risk_exp)
# 创建数据分布对象
ddist <- datadist(cgga_risk_exp)
pdf("CGGA_ggrisk.pdf")
ggrisk(multi_cox,code.highrisk = 'High rtCDI',#高风险标签，默认为 ’High’
       code.lowrisk = 'Low rtCDI', #低风险标签，默认为 ’Low’
       title.A.ylab='Risk Score', #A图 y轴名称
       title.B.ylab='Survival Time(year)', #B图 y轴名称，注意区分year month day
       title.A.legend='Risk Group', #A图图例名称
       title.B.legend='Status', #B图图例名称
       title.C.legend='Expression', #C图图例名称
       relative_heights=c(0.1,0.1,0.01,0.15), #A、B、热图注释和热图C的相对高度
       color.A=c(low='#88c4e8',high='#c74546'),#A图中点的颜色
       color.B=c(code.0='#88c4e8',code.1='#c74546'), #B图中点的颜色
       color.C=c(low='#88c4e8',median='#f6f5ee',high='#c74546'), #C图中热图颜色
       vjust.A.ylab=1, #A图中y轴标签到y坐标轴的距离,默认是1
       vjust.B.ylab=2 #B图中y轴标签到y坐标轴的距离,默认是2
)
dev.off()

#计算风险评分
risk=as.character()
for(i in 1:nrow(cgga_risk_exp)){risk=append(risk,sum(cgga_risk_exp[i,multi_cox_gene]*tcga_wight))}
cgga_risk_exp$Risk_Score=as.numeric(risk)
cgga_risk_exp$Risk_Level=ifelse(cgga_risk_exp$Risk_Score> median(cgga_risk_exp$Risk_Score),"H","L")

table(cgga_risk_exp$Risk_Level)
save(cgga_risk_exp,file = "cgga_lgg_risk_exp.RData")
surv_pvalue(survfit(Surv(os_time,os)~Risk_Level ,data = cgga_risk_exp))$pval


pdf("CGGA_lgg_ggsurvplot_roc_ggrisk.pdf")
with(cgga_risk_exp,
     ROC_riskscore <<-
timeROC(T=cgga_risk_exp$os_time, 
        delta=cgga_risk_exp$os,
        marker=cgga_risk_exp$Risk_Score,
        cause=1,weighting="cox",
        times=c(1,2,3,5),ROC=TRUE))

#t=1       t=2       t=3       t=5 
##0.7917526 0.8446370 0.8282861 0.7826257
plot(ROC_riskscore, time = 1, col = "#c74546", add = F,title = "")
plot(ROC_riskscore, time = 2, col = "#88c4e8", add = T)
plot(ROC_riskscore, time = 3, col = "#93cc82", add = T)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)),cex=1.2)
text(0.5,0.15,paste("2-Year AUC = ",round(ROC_riskscore$AUC[2],3)),cex=1.2)
text(0.5,0.1,paste("3-Year AUC = ",round(ROC_riskscore$AUC[3],3)),cex=1.2)
##cgga km曲线
ggsurvplot(survfit(Surv(os_time,os) ~ Risk_Level,
                   data = cgga_risk_exp),
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           legend.title ="rtCDI",
           palette = c("#CA0020","#0571B0"),
           xlab="Time (years)",
           ylab = "Overall survival probability",
           pval.method = T,legend.labs = c("High", "Low"))
dev.off()

########其他肿瘤验证##################
TCGA.BRCA.GDC_phenotype <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/TCGA-BRCA.GDC_phenotype.tsv", row.names=1)

TCGA.BRCA.htseq_fpkm <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/TCGA-BRCA.htseq_fpkm.tsv", row.names=1)

TCGA.BRCA.survival <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/TCGA-BRCA.survival.tsv", row.names=1)
`gencode.v22.annotation.gene.(1)` <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/gencode.v22.annotation.gene (1).probeMap", row.names=1)
load("multi_cox_coefficients.RData")

gtf23 = `gencode.v22.annotation.gene.(1)`
pheno= TCGA.BRCA.GDC_phenotype
pheno$sample =  rownames(pheno)
pheno$sample <- gsub("-",".",pheno$sample)

# pheno$days_to_last_follow_up.diagnoses=as.numeric(pheno$days_to_last_follow_up.diagnoses)
# pheno$days_to_death.demographic=as.numeric(pheno$days_to_death.demographic)
pheno=pheno[pheno$radiation_therapy %in% "YES",]
pheno=unique(pheno[,c(140,77)])
colnames(pheno) = c("sample","radiation_therapy")
rownames(pheno) = pheno$sample
#table(pheno$radiation_therapy)
rownames(TCGA.BRCA.survival) <- gsub("-",".",rownames(TCGA.BRCA.survival))
comgene <- intersect(pheno$sample,rownames(TCGA.BRCA.survival))
pheno = pheno[comgene,]
TCGA.BRCA.survival = TCGA.BRCA.survival[comgene,]
#save(pheno,file = "tcga_GBM_phen.RData")
# surv=phen[,c(10:13)]
# save(surv,file = "tcga_lgg_surv.RData")

#fpkm转tpm
fpkm = TCGA.BRCA.htseq_fpkm
fpkmToTpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
tpm =  apply(fpkm,2,fpkmToTpm)
tpm1 = as.data.frame(tpm)
min(tpm1)
boxplot(tpm1[,1:10])
tpm1=log(tpm1+1)
table(colnames(tpm1) %in% pheno$sample)
tpm1=tpm1[,colnames(tpm1) %in% pheno$sample]
table(rownames(tpm1) %in% rownames(gtf23))
tpm1$gene_id=rownames(tpm1)
tpm1$gene_name = gtf23[match(tpm1$gene_id,rownames(gtf23)),"gene"] ## ID转换

tpm2 = limma::avereps(tpm1[,!colnames(tpm1) %in% "gene"],ID = tpm1$gene_name) ## 去除重复gene_name，重复的取平均值
tpm2 = as.data.frame(tpm2) 
#tpm2=merge(gtf23[,-2],tpm1,by="gene_id")
rownames(tpm2) = tpm2$gene_name
tpm2 = tpm2[,-593]

exp = as.data.frame(t(tpm2[rownames(multi_cox_coefficients),colnames(tpm2) %in% pheno$sample]))
# comgene <- intersect(colnames(exp),rownames(TCGA.BRCA.survival))
# exp = exp[,comgene]
# TCGA.BRCA.survival = TCGA.BRCA.survival[comgene,]
exp$os_time=TCGA.BRCA.survival[match(rownames(exp),rownames(TCGA.BRCA.survival)),"OS.time"]
exp$os=TCGA.BRCA.survival[match(rownames(exp),rownames(TCGA.BRCA.survival)),"OS"]
exp=exp[exp$os_time > 0,]
exp=na.omit(exp)
a = rownames(exp)
exp <- as.data.frame(sapply(exp, as.numeric))
rownames(exp) = a
save(exp,file = "BRCA_train_data.RData")
multi_cox_gene = rownames(multi_cox_coefficients)
#risk_BRCA=exp[,c("os_time","os",as.character(multi_cox_gene))]

risk=as.character()
for(i in 1:nrow(exp)){risk=append(risk,sum(exp[i,multi_cox_gene]*tcga_wight))}

exp$Risk_Score=as.numeric(risk)
exp=na.omit(exp)
exp$Risk_Level=ifelse(exp$Risk_Score> median(exp$Risk_Score),"H","L")
#res.cut <- surv_cutpoint(risk_GBM, time = "os_time", event = "os", variables = c("Risk_Score"))
#risk_exp$Risk_Level=ifelse(risk_exp$Risk_Score> res.cut$cutpoint[[1]],"H","L")
table(exp$Risk_Level)
risk_exp = exp
save(risk_exp,file = "BRCA_risk_exp.RData")

surv_pvalue(survfit(Surv(os_time, os)~ Risk_Level ,data = risk_exp))$pval
pdf("BRCA_ggsurvplot_roc_ggrisk.pdf")
with(risk_exp,
     ROC_riskscore <<- timeROC(T=risk_exp$os_time, 
                               delta=risk_exp$os,
                               marker=risk_exp$Risk_Score,
                               cause=1,
                               weighting="cox",
                               times=c(365,730,1095,1825),ROC=TRUE))
#t=1       t=2       t=3       t=5 
#0.9175524 0.9067425 0.8921965 0.8051401 
plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 730, col = "blue", add = T)
plot(ROC_riskscore, time = 1095, col = "purple", add = T)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)),cex=1.2)
text(0.5,0.15,paste("2-Year AUC = ",round(ROC_riskscore$AUC[2],3)),cex=1.2)
text(0.5,0.1,paste("3-Year AUC = ",round(ROC_riskscore$AUC[3],3)),cex=1.2)
##tcga km曲线
ggsurvplot(survfit(Surv(os_time,os) ~ Risk_Level,
                   data = risk_exp),
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           legend.title ="Risk_Level",
           palette = c("#CA0020","#0571B0"),
           xlab="Time (years)",
           ylab = "Overall survival probability",
           pval.method = T,legend.labs = c("High", "Low"))
dev.off()

## nomo                    ###############################################
rm(list=ls())

options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_lgg.RData")
phen=phen[,c(2:3,5,8,10,11,14,15)]
colnames(phen)[8] = 'rtCDI'

table(phen$Gender)
phen=phen[order(phen$rtCDI),]
phen$Age=as.numeric(phen$Age)
phen$Gender=as.numeric(phen$Gender)-1
phen$Grade <- ifelse(phen$Grade == "G2", 2,
                                ifelse(phen$Grade == "G3", 3,1))
save(phen,file="lgg_phen_rtCDI.RData") 
load("lgg_phen_rtCDI.RData")

pdf("lgg_ggforest_phen_risk_score.pdf")
ggforest(coxph(Surv(os_time, os) ~. ,data = phen[,c(1:3,5,6,8)]),main = "TCGA Cohort",cpositions = c(0.05,0.15,0.35),refLabel = "reference")
dev.off()

#multiCox <- coxph(Surv(OS.time, OS) ~ ., data =  train_surv_expr_vennlasso)

library(rms)
library(foreign)
library(survival)
colnames(phen)
pheno1=phen[,c("os_time","os","rtCDI","Age","Gender")]
##对数据进行打包
ddist <- datadist(pheno1)
options(datadist='ddist')
units(pheno1$os_time) <- "year"

f <- coxph(Surv(os_time,os)~ rtCDI+Age+Gender, data=pheno1) ##为什么这里有两个构建模型的代码，
library(regplot)
pdf("lgg_risk_nomo2.pdf",10,8)
regplot(f,plots = c("density", "no plot"),observation=pheno1[12,],failtime = c(1, 2, 3), prfail = TRUE,points = TRUE )
dev.off()
## 评价COX回归的预测效果
## 计算c-index
rcorrcens(Surv(os_time,os) ~ predict(f), data = pheno1)
# nom <- nomogram(f, fun=list(function(x) surv(1, x),function(x) surv(2, x),function(x) surv(3, x)),
#                 fun.at=c(0.05,seq(0.05,0.95,by=0.1),0.95),           
#                 funlabel=c("1-year Survival Probability","2-year Survival Probability", "3-year Survival Probability"))
# pdf("lgg_risk_nomo_score.pdf",10,8)
# plot(nom, xfrac=0.4)
# dev.off()

# validate(f, method="boot", B=1000, dxy=T)
# rcorrcens(Surv(os_time, os) ~ predict(f), data = pheno1)
#             C     Dxy   aDxy  SD    Z   P   n
# predict(f) 0.756 0.513 0.513 0.029 17.6 0 482

f1 <- cph(Surv(os_time, os)~ rtCDI+Age+Gender,x=T, y=T, surv=T, time.inc = 1, data=pheno1) 


DynNom(f,pheno1)

#cal1 <- calibrate(f1, cmethod="KM", method="boot", u= 5, m= 150, B= 1000)
# 使用calibrate函数创建一个校准对象cal1
cal1 <- calibrate(f1, 
                  cmethod = 'KM',   #表示使用Kaplan-Meier（KM）估计方法进行校准
                  method = "boot",  # 表示使用自助法（Bootstrap）进行校准，Bootstrap 是一种统计方法，它通过从原始数据中有放回地进行重采样来估计参数的不确定性和分布。在这里，Bootstrap 用于生成多个随机样本来估计校准曲线的分布，以便获得更可靠的校准结果。
                  u = 1,          # 设置时间间隔，需要与之前模型中定义的time.inc一致
                  m = 90,           # 每次抽样的样本量，根据样本量来确定，标准曲线一般将所有样本分为3组（在图中显示3个点）
                  B = 1000)         # 抽样次数
f2 <- cph(Surv(os_time, os)~ rtCDI+Age+Gender,x=T, y=T, surv=T, time.inc = 2, data=pheno1) 
cal2 <- calibrate(f2, cmethod="KM", method="boot", u= 5, m= 90, B= 1000)
f3 <- cph(Surv(os_time, os)~ rtCDI+Age+Gender,x=T, y=T, surv=T, time.inc = 3, data=pheno1) 
cal3 <- calibrate(f3, cmethod="KM", method="boot", u= 5, m= 90, B= 1000)
 


pdf("lgg_nomo_calibration_score.pdf")
plot(cal1,add=F,conf.int=T,subtitles = T,cex.subtitles=0.8,lwd=2,lty=1,errbar.col="#0571B0",
     xlim=c(0.0,1),ylim=c(0.0,1),xlab="", ylab="",col="red")
plot(cal2,add=T,conf.int=T,subtitles = F,cex.subtitles=0.8,  lwd=2, lty=1, errbar.col="#92C5DE", 
     xlim=c(0.0,1),ylim=c(0.0,1),xlab="", ylab="", col="#407600")
plot(cal3,add=T,conf.int=T,subtitles = F,cex.subtitles=0.8,  lwd=2, lty=1, errbar.col="#CA0020", 
     xlim=c(0.0,1),ylim=c(0.0,1),xlab="", ylab="", col="#407600")
legend("topleft", legend=c("1 year","2 year", "3 year"), col=c("#0571B0","#92C5DE","#CA0020"), lwd=2)
abline(0,1,lty=3,lwd=1,col="grey")
dev.off()

source("G:/fangshe/cell_death_qj2/death_pathway3/STDCA/stdca.R")
ff <- cph(Surv(os_time, os) ~ rtCDI, data=pheno1,x=T, y=T, surv=T)
pheno1$one.years.Survival.Probabilitynew = c(1- (summary(survfit(ff, newdata = pheno1),times = 1)$surv))
pheno1$two.years.Survival.Probabilitynew = c(1- (summary(survfit(ff, newdata = pheno1),times = 2)$surv))
pheno1$three.years.Survival.Probabilitynew = c(1- (summary(survfit(ff, newdata = pheno1),times = 3)$surv))

pdf("lgg_nomo_stdca_score.pdf")
stdca(data=pheno1,outcome="os", ttoutcome="os_time", timepoint=1,predictors="one.years.Survival.Probabilitynew",xstop=0.2, smooth=TRUE)
stdca(data=pheno1,outcome="os", ttoutcome="os_time", timepoint=2,predictors="two.years.Survival.Probabilitynew",xstop=0.6, smooth=TRUE)
stdca(data=pheno1,outcome="os", ttoutcome="os_time", timepoint=3,predictors="three.years.Survival.Probabilitynew",xstop=0.8,smooth=TRUE)
dev.off()


##单因素cox##################
rm(list=ls())

load("G:/fangshe/cell_death_qj2/death_pathway3/km_cox.RData")
#load("G:/fangshe/cell_death_qj2/death_pathway3/train_data.RData")
#load("G:/fangshe/cell_death_qj2/death_pathway3/test_data.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_rr_rs_de.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/gene_list.RData")

p_v=0.05   
km_p=0.01   
cor_v=0.5 
res$logFC=-res$logFC

res=res[rownames(res) %in% gene_list$Gene,]##取交集
de_gene=rownames(res[res$adj.P.Val<p_v & abs(res$logFC)> cor_v,])##筛选P值为0.05logfc大于0.5的
de_gene_filt = km_cox[km_cox$gene %in% de_gene & km_cox$cox_p < 0.05,]##筛选啥的

# tabletext[c("V2", "V3", "V4")] <- lapply(tabletext[c("V2", "V3", "V4")], function(x) as.numeric(as.character(x)))

# 绘制森林图
# 输入表格的制作
library(forestplot)
# tabletext <- cbind(c("Gene",de_gene_filt$gene),
#                    c("HR",format(round(as.numeric(de_gene_filt$cox_hr),3),nsmall = 3)),
#                    c("lower 95%CI",format(round(as.numeric(de_gene_filt$low_CI),3),nsmall = 3)),
#                    c("upper 95%CI",format(round(as.numeric(de_gene_filt$up_CI),3),nsmall = 3)),
#                    c("pvalue",format(round(as.numeric(de_gene_filt$cox_p),3),nsmall = 3)))
# tabletext = as.data.frame(tabletext)
hz <- paste(round(de_gene_filt$cox_hr,3),
            "(",round(de_gene_filt$low_CI,3),
            "-",round(de_gene_filt$up_CI,3),")",sep = "")


tabletext <- cbind(c(NA,"Gene",de_gene_filt$gene),
                   c(NA,"Coefficient",round(de_gene_filt$cox_hr,3)),
                   c(NA,"P value",ifelse(de_gene_filt$cox_p<0.05,"P < 0.05",round(de_gene_filt$cox_p,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
pdf("forestplot.pdf",10,10)
forestplot(labeltext=tabletext, 
           graph.pos=3,  #为Pvalue箱线图所在的位置
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,de_gene_filt$cox_hr),
           lower=c(NA,NA,de_gene_filt$low_CI ), #95%置信区间下限
           upper=c(NA,NA,de_gene_filt$up_CI), #95%置信区间上限
           boxsize=0.3,lwd.ci=2,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=1,      #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0.5, 1,1.5), #横坐标刻度
           lwd.xaxis=1,            #X轴线宽
           lineheight = unit(1.5,"cm"), #固定行高
           graphwidth = unit(.3,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                           "39" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.5),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")
dev.off()
## TCGA                    ###############################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/TCGA-LGG_clinical_all_list.Rdata")
clinical_drug=clinical$clinical_drug
table(clinical_drug$therapy_types)
samp=unique(as.character(clinical_drug[clinical_drug$therapy_types %in% "Immunotherapy","bcr_patient_barcode"]))

load("lgg_sva_risk_exp.RData")
risk_exp=risk_exp[substr(rownames(risk_exp),1,12) %in% samp,]
risk_exp$Risk_Level=ifelse(risk_exp$Risk_Score>median(risk_exp$Risk_Score),"H","L")
surv_pvalue(survfit(Surv(os_time, os)~ Risk_Level ,data = risk_exp))$pval
timeROC(T=risk_exp$os_time, delta=risk_exp$os,marker=risk_exp$Risk_Score,cause=1,weighting="cox",times=c(1,2,3,5),ROC=TRUE)$AUC


library(survivalROC)
auc1= survivalROC(Stime=risk_exp$os_time,status=risk_exp$os,marker = risk_exp$Risk_Score,predict.time = 1, span = 0.25*nrow(risk_exp)^(-0.20),method = "KM")
auc2= survivalROC(Stime=risk_exp$os_time,status=risk_exp$os,marker = risk_exp$Risk_Score,predict.time = 2, span = 0.25*nrow(risk_exp)^(-0.20),method = "KM")
auc3= survivalROC(Stime=risk_exp$os_time,status=risk_exp$os,marker = risk_exp$Risk_Score,predict.time = 3, span = 0.25*nrow(risk_exp)^(-0.20),method = "KM")
auc5= survivalROC(Stime=risk_exp$os_time,status=risk_exp$os,marker = risk_exp$Risk_Score,predict.time = 5, span = 0.25*nrow(risk_exp)^(-0.20),method = "KM")
message(paste("auc1",round(auc1$AUC,3),"auc2",round(auc2$AUC,3),"auc3",round(auc3$AUC,3),"auc5",round(auc5$AUC,3),sep = "\t"))

pdf("lgg_tcga_immunotherapy_ggsurvplot_roc_ggrisk.pdf")
ggsurvplot(survfit(Surv(os_time,os) ~ Risk_Level, data = risk_exp),pval = TRUE,risk.table = TRUE,risk.table.col = "strata",linetype = "strata",ggtheme = theme_bw(),legend.title ="Risk_Level",palette = c("#CA0020","#0571B0"),xlab="Time (years)",ylab = "Overall survival probability",pval.method = T,legend.labs = c("High", "Low"))
plot(auc5$FP, auc5$TP, type="l",col="#0571B0", xlim=c(0,1), ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",main="TCGA immunotherapy Cohort",lwd=3)
lines(auc2$FP, auc2$TP, type="l",col="#F4A582", xlim=c(0,1), ylim=c(0,1),lwd=3)
lines(auc3$FP, auc3$TP, type="l",col="#CA0020", xlim=c(0,1), ylim=c(0,1),lwd=3)
abline(0,1,col="gray",lty=4,lwd=3)
legend(0.4,0.4,c(paste("Year = 2,  AUC =",round(auc2$AUC,4)), 
                 paste("Year = 3,  AUC =",round(auc3$AUC,4)),
                 paste("Year = 5,  AUC =",round(auc5$AUC,4))),lty= 1 ,lwd= 2,
       col=c("#F4A582","#CA0020","#0571B0"),bty = "n",seg.len=1,cex=1)
dev.off()



risk_exp$measure_of_response=clinical_drug[match(substr(rownames(risk_exp),1,12),clinical_drug$bcr_patient_barcode),"measure_of_response"]
risk_exp$measure_of_response=as.character(risk_exp$measure_of_response)
chisq.test(table(risk_exp$measure_of_response,risk_exp$Risk_Level))

ggplot(risk_exp[risk_exp$measure_of_response != "",], aes(measure_of_response, Risk_Score)) + 
  geom_boxplot(aes(fill = measure_of_response))+xlab(label = "")+
  stat_compare_means(aes(group=measure_of_response))+
  geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  scale_fill_brewer(palette = "Set1")+ylab("Risk_Score")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))


##########################################################################
##################### 作图 ###############################################
##########################################################################
## gene enrichment         ##################################
rm(list=ls())

load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_rr_rs_de.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/gene_list.RData")

p_v=0.01   
km_p=0.01   
cor_v=0.5 
res=res[rownames(res) %in% gene_list$Gene,]
de_gene=res[res$adj.P.Val<0.05 & abs(res$logFC)>0.5,]
write.table(de_gene,sep = "\t",quote = F,row.names = F)
de_gene = rownames(de_gene)
library(clusterProfiler)
library(org.Hs.eg.db)
ensem=bitr(de_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
go=enrichGO(ensem$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,keyType = "ENTREZID")
go=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')##enterID转化为基因名
kegg=enrichKEGG(ensem$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,qvalueCutoff = 0.8,pAdjustMethod = "none")
kegg = setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")##enterID转化为基因名

save(go,kegg,file = "genelist_rrrs_de_gene_enrich.RData")
table(go[,1]) #查看BP,CC,MF的统计数目

##可视化
#1、go柱状图
go=go[order(go@result$qvalue),]
go_MF<-go[go$ONTOLOGY=="MF",][1:10,]
go_CC<-go[go$ONTOLOGY=="CC",][1:10,]
go_BP<-go[go$ONTOLOGY=="BP",][1:10,]
go_enrich_df<-data.frame(ID=c(go_BP$ID, go_CC$ID, go_MF$ID),
                         Description=c(go_BP$Description, go_CC$Description, go_MF$Description),
                         GeneNumber=c(go_BP$Count, go_CC$Count, go_MF$Count),
                         type=factor(c(rep("biological process", 10), rep("cellular component", 10),rep("molecular function",10)),levels=c("molecular function", "cellular component", "biological process")))
## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=5, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
# Assuming go_enrich_df$Description is a character vector
labels <- sapply(go_enrich_df$Description, function(desc) {
  if (!is.na(desc)) {
    # Shorten the names or perform other modifications here
    shortened_name <- shorten_names(desc)
    return(shortened_name)
  } else {
    # Handle NA values if necessary
    return("NA_Label")  # Replace with your desired label for NAs
  }
})
names(labels) = rev(1:nrow(go_enrich_df))
## colors for bar // green, blue, orange

CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
pdf("de_PCD_goenrih.pdf",8,6)
ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_test() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
dev.off()
#coord_flip(...)横向转换坐标：把x轴和y轴互换，没有特殊参数

#1、kegg柱状图
###KEGG可视化###

#数据导入#

#数据处理#

display_number = 25#显示数量设置

kk_result = as.data.frame(kegg)[1:display_number[1], ]

kk = as.data.frame(kk_result)

rownames(kk) = 1:nrow(kk)

kk$order=factor(rev(as.integer(rownames(kk))),labels = rev(kk$Description))

library(ggplot2)
pdf("de_PCD_keggenrih.pdf",8,6)
#柱状图#

ggplot(kk,aes(y=order,x=Count,fill=pvalue))+
  
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "#8DA1CB",high ="#66C3A5" )+
  labs(title = "KEGG Pathways Enrichment",  #设置标题、x轴和Y轴名称
       x = "Gene number",
       y = "Pathway")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
dev.off()
#coord_flip(...)横向转换坐标：把x轴和y轴互换，没有特殊参数



#2、基因-通路关联网络图
#GO通路-基因网络图
y <- c("Fatty acid biosynthesis","Fatty acid metabolism","PPAR signaling pathway","Insulin signaling pathway","Fatty acid degradation","Glucagon signaling pathway","Insulin resistance")
go_resule = go@result
foldchange = de_gene$logFC
names(foldchange) <- rownames(de_gene)
head(foldchange)
# 把想要展示的通路名称赋值给y#
# cnetplot(go, 
#          foldChange =foldchange, 
#          showCategory = 5, ## 展示前5个通路
#          node_label = "all", # category | gene | all | none
#          colorEdge = TRUE,
#          circular=TRUE,)
p1 <- cnetplot(
  go,
  foldChange = foldchange,
  showCategory = 7, 
  node_label = "all",
  colorEdge = TRUE,
  circular = TRUE
  # max.overlaps = 100, # Increase max.overlaps to avoid unlabeled data points
  )
p1+scale_color_gradientn(colours = c("#407600","#0571B0", "#CA0020","#F4A582","#88c4e8"))
dev.off()


##gene maftools ###################
load("G:/fangshe/cell_death_qj2/death_pathway3/gene_list.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_mut.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_rr_rs_de.RData")
mut$Tumor_Sample_Barcode=substr(mut$Tumor_Sample_Barcode,1,15)
res=res[rownames(res) %in% gene_list$Gene,]
de_gene=res[res$adj.P.Val<0.05 & abs(res$logFC)>0.5,]

mut$Tumor_Sample_Barcode=substr(mut$Tumor_Sample_Barcode,1,15)

library(maftools)
vc_cols = c("#33A02C","#FF7F00","#6A3D9A","#B15928","#1F78B4","#E31A1C","black")
names(vc_cols) = c('Missense_Mutation','Nonsense_Mutation','Splice_Site','Frame_Shift_Del',"CNV_loss","CNV_gain",'Multi_Hit')
var_maf4 = read.maf(maf= mut[mut$Hugo_Symbol %in% rownames(de_gene),])

pdf("gene_maftools.pdf",12,8)
oncoplot(maf = var_maf4,genes = rownames(de_gene),top=10)

plotmafSummary(maf = var_maf4,rmOutlier = TRUE, addStat = "median",dashboard = TRUE,titvRaw = FALSE,titleSize = c(1.2, 1))
dev.off()


## risk PCA                #############################################################
rm(list=ls())
load("G:/fangshe/cell_death_qj2/death_pathway3/gene_list.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_tcga_risk_exp.RData")

tcga_risk_exp=tcga_risk_exp[rownames(tcga_risk_exp) %in% colnames(tpm_sva),]

library(Rtsne)
library(ggplot2)
tsne_out=Rtsne(t(tpm_sva[rownames(tpm_sva) %in% gene_list$Gene,rownames(tcga_risk_exp)]),pca=FALSE) #dim=3
dat_plt=as.data.frame(tsne_out$Y)
rownames(dat_plt)=rownames(tcga_risk_exp)
dat_plt$txt=rownames(dat_plt)
dat_plt$Risk_Level=tcga_risk_exp[match(dat_plt$txt,rownames(tcga_risk_exp)),"Risk_Level"]
dat_plt$Risk_Level=factor(dat_plt$Risk_Level,labels=c("Low","High"),levels= c("L","H"))

pdf("lgg_Risk_Level_pca.pdf")
ggplot(dat_plt,aes(x=dat_plt[,1],y=dat_plt[,2],colour=Risk_Level))+geom_point()+
  scale_color_manual(values=c("#CA0020","#0571B0","#4DAF4A"))+
  theme(plot.title = element_text(hjust = 0.5))+labs(x = "",y = "",title ="tSNE Plot" )

ggplot(data = dat_plt, aes(x=V1, y=V2,color=Risk_Level)) + geom_point(size=2,aes(color =  dat_plt$Risk_Level)) +
  stat_ellipse(aes(fill = dat_plt$Risk_Level), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  scale_color_manual(values =c("#CA0020","#0571B0","#4DAF4A")) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'),  plot.title = element_text(hjust = 0.5),legend.title = element_blank()) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5)+
  labs(x = "",y = "",title ="" )
dev.off()

## risk HL go_kegg_gsea    ###################################
rm(list=ls())
library(forcats)

load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
library(limma)
load("G:/fangshe/cell_death_qj2/death_pathway3/phen_sva_lgg.RData")

#load("G:/fangshe/cell_death_qj2/death_pathway3/kegg_pathway_category")
phen_sva=phen_sva[phen_sva$Sample %in% colnames(tpm_sva),]
table(phen_sva$Sample %in% colnames(tpm_sva))
#1、差异分析
mm=model.matrix(~0+as.character(phen_sva$Risk_Level))
colnames(mm)=c("H","L")
res=lmFit(tpm_sva[,phen_sva$Sample],mm) %>%
  contrasts.fit(., makeContrasts(H-L, levels = mm))%>%
  eBayes(.,trend = T) %>% 
  topTable(., number =Inf)%>% 
  .[order(.$adj.P.Val),] %>%
  na.omit()
table(res$adj.P.Val<0.05 )
table(res$adj.P.Val<0.05 & abs(res$logFC)>0.5)
save(res,file="lgg_sva_HL_de_ssgsea_kegg_res.RData")


#2、GO、KEGG、GSEA分析
load("lgg_sva_HL_de_ssgsea_kegg_res.RData")
de_gene=rownames(res[res$adj.P.Val<0.05 & abs(res$logFC)>0.5,])
de_gene_df = res[res$adj.P.Val<0.01 & abs(res$logFC)>0.5,]
library(clusterProfiler)
library(org.Hs.eg.db)
#2.1 GO 和KEGG

ensem=bitr(gene,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")
#GO、KEGG富集
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。

go=enrichGO(ensem$ENTREZID,
            OrgDb = org.Hs.eg.db,
            ont = "ALL",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.02,
            keyType = "ENTREZID")
kegg=enrichKEGG(ensem$ENTREZID, 
                organism = 'hsa', 
                keyType = 'kegg', 
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

kegg = setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
save(go,kegg,file = "risk_HL_de_gene_enrich.RData")

##可视化GO、KEGG
load("risk_HL_de_gene_enrich.RData")
simple_go <- simplify(go) #去除冗余的GO
kegg_result <- kegg@result
go_result <- simple_go@result

##可视化
#1、go柱状图
go =  go_result[order(go_result$p.adjust),]
go_MF<-go[go$ONTOLOGY=="MF",][1:10,]
go_CC<-go[go$ONTOLOGY=="CC",][1:10,]
go_BP<-go[go$ONTOLOGY=="BP",][1:10,]
go_enrich_df<-data.frame(ID=c(go_BP$ID, go_CC$ID, go_MF$ID),
                         Description=c(go_BP$Description, go_CC$Description, go_MF$Description),
                         GeneNumber=c(go_BP$Count, go_CC$Count, go_MF$Count),
                         type=factor(c(rep("biological process", 10), rep("cellular component", 10),rep("molecular function",10)),levels=c("molecular function", "cellular component", "biological process")))
## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=6, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
# Assuming go_enrich_df$Description is a character vector
labels <- sapply(go_enrich_df$Description, function(desc) {
  if (!is.na(desc)) {
    # Shorten the names or perform other modifications here
    shortened_name <- shorten_names(desc)
    return(shortened_name)
  } else {
    # Handle NA values if necessary
    return("NA_Label")  # Replace with your desired label for NAs
  }
})
names(labels) = rev(1:nrow(go_enrich_df))
## colors for bar // green, blue, orange

CPCOLS <- c("#88c4e8", "#c74546", "#93cc82")
library(ggplot2)
pdf("de_risk_HL_goenrih.pdf",8,6)
ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_test() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
dev.off()
#coord_flip(...)横向转换坐标：把x轴和y轴互换，没有特殊参数

#1、kegg柱状图
###KEGG可视化###

#数据导入#

#数据处理#

display_number = 30#显示数量设置

kk_result = as.data.frame(kegg)[1:display_number[1], ]

kk = as.data.frame(kk_result)

rownames(kk) = 1:nrow(kk)

kk$order=factor(rev(as.integer(rownames(kk))),labels = rev(kk$Description))

library(ggplot2)
pdf("de_risk_HL_keggenrih.pdf",8,6)
#柱状图#

ggplot(kk,aes(y=order,x=Count,fill=pvalue))+
  
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "#88c4e8",high ="#c74546" )+
  labs(title = "KEGG Pathways Enrichment",  #设置标题、x轴和Y轴名称
       x = "Gene number",
       y = "Pathway")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
dev.off()
#coord_flip(...)横向转换坐标：把x轴和y轴互换，没有特殊参数


#2.2 GSEA
rm(list=ls())
library(msigdbr)
load("lgg_sva_HL_de_ssgsea_kegg_res.RData")
res$genesymbol = rownames(res)
# 多个gene symbol的直接删除，方便演示
res<- res[!grepl("/",res$genesymbol),]
gene = rownames(res)
ensem=bitr(gene,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")
gene_entrezid <- merge(ensem,res,by.x = "SYMBOL", by.y = "genesymbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist,decreasing = T)

head(genelist)
#使用msigdbr包从msigdb数据库下载人类的C5注释集，大家常用的GO、KEGG的数据其实都是包括在msigdb数据库中的
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

set.seed(456)
gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH"
)
save(gsea_res,file ="de_risk_HL_gsea_result.RData",row.names = T)

#可视化
load("G:/fangshe/cell_death_qj3/cell_death_qj3/de_risk_HL_gsea_result.RData")
library(enrichplot)
library(ggplot2)
gsea = gsea_res@result
##桥图
ridgeplot(gsea_res,
          showCategory = 20,
          fill = "p.adjust", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE
)+
  theme(axis.text.y = element_text(size=8))

##一张图展示多条通路
geneSetID = c('GOBP_MITOTIC_SISTER_CHROMATID_SEGREGATION',
              'GOBP_SISTER_CHROMATID_SEGREGATION',
              'GOBP_B_CELL_MEDIATED_IMMUNITY',
              'GOBP_T_CELL_ACTIVATION',
              'GOCC_DENDRITIC_SHAFT',
              'GOBP_NEUROTRANSMITTER_TRANSPORT')

pdf("de_risk_HL_GSEA.pdf",11,8)
# change line color
gseaNb(object = gsea_res,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.05,pvalY = 0.05,
       curveCol = jjAnno::useMyCol('paired',6))
dev.off()


####GSVA
rm(list=ls())
library(TCGAbiolinks)
load("lgg_sva_HL_de_ssgsea_kegg_res.RData")
load()







## 模型基因结果富集分析     #####################################    
rm(list=ls())
load("multi_cox_coefficients.RData")
#1、加载包
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包

#2、基因id转换，kegg和go富集用的ID类型是ENTREZID）

a <- rownames(multi_cox_coefficients)
a <- as.data.frame(a)

gene.df <- bitr(a$a, fromType="SYMBOL",
                toType="ENTREZID", OrgDb='org.Hs.eg.db')
a <- inner_join(a,gene.df,by=c("a"="SYMBOL"))

#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     
gene <- gene.df$ENTREZID

#3、GO富集
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.05,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
#4、将结果保存到当前路径
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- ego_CC@result
ego_result_MF <- ego_MF@result
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)#或者这样也能得到ego_ALL一样的结果
write.csv(ego_ALL,file = "ego_ALL.csv",row.names = T)
write.csv(ego_result_BP,file = "ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "ego_result_MF.csv",row.names = T)
write.csv(ego,file = "ego.csv",row.names = T)

#5、但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
display_number = c(10, 10, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- ego_result_MF[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#8DA1CB","#FD8D62","#66C3A5")#设定颜色

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置

#1、KEGG富集
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)

#2、可视化
###柱状图
hh <- kk@result#自己记得保存结果哈！
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels=rev(hh$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱


##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(hh)){
  description_splite=strsplit(hh$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  hh$Description[i]=description_collapse
  hh$Description=gsub(pattern = "NA","",hh$Description)
}
  
ggplot(hh,aes(y=order,x=Count,fill=pvalue))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#CC0000",high ="#2f5688" )+#颜色自己可以换
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

###气泡图
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()




## de ssgsea               ###################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_ssgsea.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_lgg.RData")

gsva_matrix=apply(gsva_matrix,1,function(x){(x-mean(x)) / sd(x)})
gsva_matrix=as.data.frame(t(gsva_matrix))

gsva_matrix$Risk_Level=phen[rownames(gsva_matrix),"Risk_Level"]
gsva_matrix=na.omit(gsva_matrix)
names(gsva_matrix)[29] = 'rtCDI'
gsva_matrix$sample=rownames(gsva_matrix)
exp1=reshape2::melt(gsva_matrix,id.vars =c("sample","rtCDI"))


gen_sort=do.call(rbind,lapply(unique(exp1$variable), function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))
gen_sort=na.omit(gen_sort)
gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)

# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("H" = "#c74546", "L" = "#88c4e8")
pdf("risk_ssgsea_boxplot.pdf",10,6)
# 设置 "High" 和 "Low" 的颜色
ggplot(exp1, aes(x = Gene, y = value)) +
  geom_boxplot(aes(fill = rtCDI)) +
  xlab(label = "") +
  stat_compare_means(aes(group = rtCDI), label = "p.signif", na.rm = TRUE, hide.ns = TRUE) +
  ylab("Fraction") +
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()
## de ssgsea celldeath              ###################################
rm(list=ls())
options(stringsAsFactors = F)
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_ssgsea_death_pathway.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_lgg.RData")
# gsva_matrix=apply(gsva_matrix,1,function(x){(x-mean(x)) / sd(x)})
# load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_phen.RData")

gsva_matrix=as.data.frame(t(gsva_matrix))
rownames(gsva_matrix) <- c("Alkaliptosis", "Apoptosis", "Autophagy","Cuproptosis","Disulfidptosis","Entosis","Ferroptosis","ICD","LCD","Necroptosis","NETosis","Oxeiptosis","PANoptosis","Parthanatos","Pyroptosis")

gsva_matrix$Risk_Level=phen[rownames(gsva_matrix),"Risk_Level"]
gsva_matrix=na.omit(gsva_matrix)

gsva_matrix$sample=rownames(gsva_matrix)
exp1=reshape2::melt(gsva_matrix,id.vars =c("sample","Risk_Level"))
names(exp1)[2] <- "rtCDI"

gen_sort=do.call(rbind,lapply(unique(exp1$variable),function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))##排序

gen_sort=na.omit(gen_sort)
gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)

##
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("H" = "#c74546", "L" = "#88c4e8")
pdf("risk_ssgsea_celldeath_boxplot.pdf",10,6)
ggplot(exp1, aes(x = Gene, y = value)) +
  geom_boxplot(aes(fill = rtCDI)) +
  xlab(label = "") +
  stat_compare_means(aes(group = rtCDI), label = "p.signif", na.rm = TRUE, hide.ns = TRUE) +
  ylab("Score") +
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(text=element_text(size=17))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

##km曲线
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_ssgsea_death_pathway.RData")
# load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_lgg.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
gsva_matrix=as.data.frame(t(gsva_matrix))
rownames(gsva_matrix) <- c("Alkaliptosis", "Apoptosis", "Autophagy","Cuproptosis","Disulfidptosis","Entosis","Ferroptosis","ICD","LCD","Necroptosis","NETosis","Oxeiptosis","PANoptosis","Parthanatos","Pyroptosis")

gsva_matrix$os = tcga_risk_exp[rownames(gsva_matrix),"os"]
gsva_matrix$os_time = tcga_risk_exp[rownames(gsva_matrix),"os_time"]
gsva_matrix=na.omit(gsva_matrix)
cl3 = gsva_matrix
gene <- c("Apoptosis")
cl3$gene <- ifelse(cl3[,gene] > median(cl3[,gene]), "high", "low")

#绘制生存曲线：
fit2 <- survfit(Surv(os_time, os) ~ gene, data = cl3)

ggsurvplot(survfit(Surv(os_time,os) ~ gene,
                   data = cl3),
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           legend.title ="Apoptosis",
           palette = c("#c74546","#88c4e8"),
           xlab="Time (years)",
           ylab = "Overall survival probability",
           pval.method = T,legend.labs = c("High", "Low"))

## de ssgsea RR_RS              ###################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_ssgsea_death_pathway.RData")
# load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_lgg.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_phen.RData")
# gsva_matrix=apply(gsva_matrix,1,function(x){(x-mean(x)) / sd(x)})
table(phen$radiation_therapy)
unique(phen$primary_therapy_outcome)
table(phen$followup_treatment_success) #使用表格函数和 unique函数检查和分析表型数据中的各列。
phen$group=ifelse(phen$primary_therapy_outcome %in% c("Complete Remission/Response"),"RS",ifelse(phen$primary_therapy_outcome %in% c("Progressive Disease"),"RR",NA))
phen=phen[!is.na(phen$group),] #您可以删除组列中有缺失值的行，并调整样本列以匹配基因表达数据的格式。
table(phen$group)
gsva_matrix=as.data.frame(t(gsva_matrix))
rownames(gsva_matrix) <- c("Alkaliptosis", "Apoptosis", "Autophagy","Cuproptosis","Disulfidptosis","Entosis","Ferroptosis","ICD","LCD","Necroptosis","NETosis","Oxeiptosis","PANoptosis","Parthanatos","Pyroptosis")

gsva_matrix$group=phen[rownames(gsva_matrix),"group"]
gsva_matrix=na.omit(gsva_matrix)
table(gsva_matrix$group)
gsva_matrix$sample=rownames(gsva_matrix)
exp1=reshape2::melt(gsva_matrix,id.vars =c("sample","group"))
# exp1$value = as.data.frame(lapply(exp1$value,as.numeric))
# colnames(exp1)[4] <- "value" 
gen_sort=do.call(rbind,lapply(unique(exp1$variable),function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))##排序

gen_sort=na.omit(gen_sort)
gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)

##
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("RS" = "#c74546", "RR" = "#88c4e8")
pdf("risk_ssgsea_RRRS_boxplot.pdf",10,6)
ggplot(exp1, aes(x = Gene, y = value)) +
  geom_boxplot(aes(fill = group)) +
  xlab(label = "") +
  stat_compare_means(aes(group = group), label = "p.signif", na.rm = TRUE, hide.ns = TRUE) +
  ylab("Score") +
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(text=element_text(size=17))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

##km曲线
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_ssgsea_death_pathway.RData")
# load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_lgg.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
gsva_matrix=as.data.frame(t(gsva_matrix))
rownames(gsva_matrix) <- c("Alkaliptosis", "Apoptosis", "Autophagy","Cuproptosis","Disulfidptosis","Entosis","Ferroptosis","ICD","LCD","Necroptosis","NETosis","Oxeiptosis","PANoptosis","Parthanatos","Pyroptosis")

gsva_matrix$os = tcga_risk_exp[rownames(gsva_matrix),"os"]
gsva_matrix$os_time = tcga_risk_exp[rownames(gsva_matrix),"os_time"]
gsva_matrix=na.omit(gsva_matrix)
cl3 = gsva_matrix
gene <- c("Apoptosis")
cl3$gene <- ifelse(cl3[,gene] > median(cl3[,gene]), "high", "low")

#绘制生存曲线：
fit2 <- survfit(Surv(os_time, os) ~ gene, data = cl3)
# ggsurvplot(
#   fit2,
#   data = cl3,
#   censor.shape="|", censor.size = 4,
#   conf.int = TRUE,
#   conf.int.style = "ribbon",
#   conf.int.alpha = 0.2,
#   pval = TRUE,
#   palette = "lancet",
#   #surv.median.line = "hv",
#   ggtheme =theme_bw(),
#   legend = "top",
#   legend.labs = c("High","Low"),
#   xlab = "OS_time(days)",
#   ylab = "Survival probablity",
#   title = "Survival curves",
#   break.x.by = 1000,
#   break.y.by = 0.2,
#   risk.table = TRUE,
#   risk.table.col = "strata",
#   risk.table.height = 0.2,
#   risk.table.y.text = FALSE
# )

ggsurvplot(survfit(Surv(os_time,os) ~ gene,
                   data = cl3),
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           legend.title ="Apoptosis",
           palette = c("#c74546","#88c4e8"),
           xlab="Time (years)",
           ylab = "Overall survival probability",
           pval.method = T,legend.labs = c("High", "Low"))

## de EPIC                 ###################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_epic_score.RData")
epic_score=as.data.frame(epic_score$cellFractions)

load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_lgg.RData")
epic_score$Risk_Level=phen[rownames(epic_score),"Risk_Level"]
epic_score=na.omit(epic_score)

epic_score$sample=rownames(epic_score)
exp1=reshape2::melt(epic_score,id.vars =c("sample","Risk_Level"))

gen_sort=do.call(rbind,lapply(unique(exp1$variable), function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))
gen_sort=na.omit(gen_sort)
gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)
exp1=exp1[exp1$variable != "otherCells",]
names(exp1)[2] = 'rtCDI'
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("H" = "#c74546", "L" = "#88c4e8")

pdf("Risk_epic_boxplot.pdf",8,5)
ggplot(exp1, aes(Gene, value)) + geom_boxplot(aes(fill = rtCDI))+xlab(label = "")+
  stat_compare_means(aes(group=rtCDI),label = "p.signif",na.rm = T,hide.ns = T)+
  scale_fill_brewer(palette = "Set1")+ylab("Fraction")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()


## de MCPcounter           ###################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_MCPcounterScore.RData")
MCPcounterScore=as.data.frame(t(MCPcounterScore))

load("phen_lgg.RData")
MCPcounterScore$Risk_Level=phen[rownames(MCPcounterScore),"Risk_Level"]
MCPcounterScore=na.omit(MCPcounterScore)

MCPcounterScore$sample=rownames(MCPcounterScore)
exp1=reshape2::melt(MCPcounterScore,id.vars =c("sample","Risk_Level"))

gen_sort=do.call(rbind,lapply(unique(exp1$variable), function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))
gen_sort=na.omit(gen_sort)
gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)
names(exp1)[2]='rtCDI'
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("H" = "#c74546", "L" = "#88c4e8")
pdf("Risk_MCPcounterScore_boxplot.pdf",7,5)
ggplot(exp1, aes(Gene, value)) + geom_boxplot(aes(fill = rtCDI))+xlab(label = "")+
  stat_compare_means(aes(group=rtCDI),label = "p.signif",na.rm = T,hide.ns = T)+
  scale_fill_brewer(palette = "Set1")+ylab("Fraction")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()
## de quantiseq            ###################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_quantiseq_score.RData")
quantiseq_score=as.data.frame(quantiseq_score)
rownames(quantiseq_score)=quantiseq_score$cell_type
quantiseq_score$cell_type=NULL
quantiseq_score=as.data.frame(t(quantiseq_score))

load("phen_lgg.RData")
quantiseq_score$Risk_Level=phen[rownames(quantiseq_score),"Risk_Level"]
quantiseq_score=na.omit(quantiseq_score)

quantiseq_score$sample=rownames(quantiseq_score)
exp1=reshape2::melt(quantiseq_score,id.vars =c("sample","Risk_Level"))

gen_sort=do.call(rbind,lapply(unique(exp1$variable), function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))
gen_sort=na.omit(gen_sort)
gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)
exp1=exp1[exp1$variable != "otherCells",]
names(exp1)[2] = 'rtCDI'
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("H" = "#c74546", "L" = "#88c4e8")
pdf("Risk_quantiseq_boxplot.pdf",7,5)
ggplot(exp1, aes(Gene, value)) + geom_boxplot(aes(fill = rtCDI))+xlab(label = "")+
  stat_compare_means(aes(group=rtCDI),label = "p.signif",na.rm = T,hide.ns = T)+
  scale_fill_brewer(palette = "Set1")+ylab("Fraction")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

## de TIMER            ###################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_timer_score.RData")


timer_score =as.data.frame(t(timer_score))

load("phen_lgg.RData")
timer_score$Risk_Level=phen[rownames(timer_score),"Risk_Level"]
timer_score=na.omit(timer_score)

timer_score$sample=rownames(timer_score)
exp1=reshape2::melt(timer_score,id.vars =c("sample","Risk_Level"))

gen_sort=do.call(rbind,lapply(unique(exp1$variable), function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))
gen_sort=na.omit(gen_sort)
gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)
exp1=exp1[exp1$variable != "otherCells",]
names(exp1)[2] = 'rtCDI'
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("H" = "#c74546", "L" = "#88c4e8")
pdf("Risk_timer_boxplot.pdf",6,4)
ggplot(exp1, aes(Gene, value)) + geom_boxplot(aes(fill = rtCDI))+xlab(label = "")+
  stat_compare_means(aes(group=rtCDI),label = "p.signif",na.rm = T,hide.ns = T)+
  scale_fill_brewer(palette = "Set1")+ylab("Fraction")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

## de CIBERSORT            ###################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_cib_raw_result.RData")
cib_result =as.data.frame(t(cib_result))
cib_result <- cib_result[rowSums(cib_result != 0) > 0, ]
cib_result <- cib_result[, colSums(cib_result != 0) > 0]
load("phen_lgg.RData")
cib_result$Risk_Level=phen[rownames(cib_result),"Risk_Level"]
cib_result=na.omit(cib_result)

cib_result$sample=rownames(cib_result)
exp1=reshape2::melt(cib_result,id.vars =c("sample","Risk_Level"))

# 使用逻辑条件筛选出某一列不包含0的行
column_name <- "value"  # 选择你想要检查的列
exp1 <- exp1[exp1[, column_name] != 0, ]

gen_sort=do.call(rbind,lapply(unique(exp1$variable), function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))
gen_sort=na.omit(gen_sort)


gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)
exp1=exp1[exp1$variable != "otherCells",]
names(exp1)[2] = 'rtCDI'
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("H" = "#c74546", "L" = "#88c4e8")
pdf("Risk_timer_boxplot.pdf",6,4)
ggplot(exp1, aes(Gene, value)) + geom_boxplot(aes(fill = rtCDI))+xlab(label = "")+
  stat_compare_means(aes(group=rtCDI),label = "p.signif",na.rm = T,hide.ns = T)+
  scale_fill_brewer(palette = "Set1")+ylab("Fraction")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

## risk maftools risk      ###################################
rm(list=ls())

##这些步骤准备突变数据和相关临床数据，以便使用 maftools 实现可视化，让您可以根据风险评分探索和分析突变模式。请记住，您需要安装并加载 maftools 软件包才能成功执行此代码。
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_mut.RData")
mut$Tumor_Sample_Barcode=substr(mut$Tumor_Sample_Barcode,1,15)

#使用变量名 mut 从 "tcga_lgg_mut.RData "中加载突变数据。从 "lgg_sva_risk_exp.RData "中加载风险表达数据，并将风险得分和样本名称一起存储在 cluster_clin 数据框中。

load("lgg_tcga_risk_exp.RData")##
cluster_clin = data.frame(Tumor_Sample_Barcode=rownames(tcga_risk_exp),Risk_Score=tcga_risk_exp$Risk_Level)

table(cluster_clin$Tumor_Sample_Barcode %in% mut$Tumor_Sample_Barcode)

#使用 table() 将 cluster_clin$Tumor_Sample_Barcode 中的样本名称与 mut$Tumor_Sample_Barcode 中的样本名称进行比较，检查是否匹配。

cluster_clin$Risk_Score=factor(cluster_clin$Risk_Score,levels = c("H","L"),labels = c("High","Low"))

#将 cluster_clin$Risk_Score 中的风险分数转换为因子变量，其中 "H "表示 "高"，"L "表示 "低"

mut=mut[mut$Tumor_Sample_Barcode %in% cluster_clin$Tumor_Sample_Barcode,]

cluster_clin=cluster_clin[cluster_clin$Tumor_Sample_Barcode %in% mut$Tumor_Sample_Barcode,]
#过滤突变数据，只保留肿瘤_样本_条形码与 cluster_clin 中样本名称相匹配的行，同样，过滤 cluster_clin 数据，只保留与过滤突变数据中的样本名称相对应的行。

library(maftools)
vc_cols = c("#33A02C","#FF7F00","#6A3D9A","#B15928","#1F78B4","#E31A1C","black")
names(vc_cols) = c('Missense_Mutation','Nonsense_Mutation','Splice_Site','Frame_Shift_Del',"CNV_loss","CNV_gain",'Multi_Hit')
#为不同的变异分类（如错义突变、无义突变等）定义名为 vc_cols 的颜色向量。使用 names() 为每个颜色类别指定有意义的名称。这些颜色将用于在可视化中区分不同类型的变异。

#读取和筛选突变数据：
var_maf1 = read.maf(maf = mut[mut$Tumor_Sample_Barcode %in% cluster_clin[cluster_clin$Risk_Score %in% "High","Tumor_Sample_Barcode"],],clinicalData = cluster_clin)

getClinicalData(x = var_maf1)

var_maf2 = read.maf(maf= mut[mut$Tumor_Sample_Barcode %in% cluster_clin[cluster_clin$Risk_Score %in% "Low","Tumor_Sample_Barcode"],],clinicalData = cluster_clin)
#使用read.maf函数来读取一个MAF(Mutation Annotation Format) 文件并且你根据一些条件从两个数据框mut和 cluster_clin筛选数据 

var_maf4 = read.maf(maf= mut,clinicalData = cluster_clin)
#用 maftools 软件包中的 read.maf() 函数创建了三个不同的 MAF 对象（var_maf1、var_maf2 和 var_maf4）。对于 var_maf1，您要根据 cluster_clin 数据过滤突变数据，只保留风险分数为 "高 "的样本。对于 var_maf2，过滤突变数据，只保留风险分数为 "低 "的样本。对于 var_maf4，则使用所有突变数据，无需过滤。

#设置配色方案
fabcolors = c("#c74546","#4d97cd")
names(fabcolors) = c("High","Low")
fabcolors = list(Risk_Score = fabcolors)

pdf("tcga_lgg_risk_HL_maftools.pdf",16,12)
oncoplot(maf = var_maf1,top = 15 ,showTumorSampleBarcodes = F ,clinicalFeatures = "Risk_Score", annotationColor = fabcolors)
oncoplot(maf = var_maf2,top = 15 ,showTumorSampleBarcodes = F ,clinicalFeatures = "Risk_Score", annotationColor = fabcolors)
oncoplot(maf = var_maf4,top = 15 ,showTumorSampleBarcodes = F ,clinicalFeatures = "Risk_Score", annotationColor = fabcolors,sortByAnnotation = TRUE)
dev.off()




load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_tmb.RData")
maf_tmb_score$Tumor_Sample_Barcode=substr(maf_tmb_score$Tumor_Sample_Barcode,1,15)
table(maf_tmb_score$Tumor_Sample_Barcode %in% cluster_clin$Tumor_Sample_Barcode)

maf_tmb_score=maf_tmb_score[maf_tmb_score$Tumor_Sample_Barcode %in% cluster_clin$Tumor_Sample_Barcode,]
maf_tmb_score$Risk_Score=as.character(cluster_clin[match(maf_tmb_score$Tumor_Sample_Barcode,cluster_clin$Tumor_Sample_Barcode),"Risk_Score"])
maf_tmb_score=maf_tmb_score[!is.na(maf_tmb_score$Risk_Score),]

pdf("tcga_lgg_risk_tmb_boxplot.pdf")
# 定义 Risk_Score 值的颜色映射
score_colors <- c("Low" = "#88c4e8", "High" = "#c74546")
ggplot(maf_tmb_score, aes(Risk_Score, TMB)) + 
  geom_boxplot(aes(fill = Risk_Score)) +
  geom_jitter(width = 0.3, alpha = 0.7, color = 'black') +
  stat_compare_means(aes(group = Risk_Score), hide.ns = TRUE) +
  # 使用手动指定的颜色映射
  scale_fill_manual(values = score_colors) +
  labs(title = "", x = "", y = "TMB") +
  ylim(0,6)

ggplot(maf_tmb_score, aes(Risk_Score, MATH)) +
  geom_boxplot(aes(fill=Risk_Score))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=Risk_Score),hide.ns = T)+# 使用手动指定的颜色映射
  scale_fill_manual(values = score_colors) +
  labs(title = "")+xlab(label = "")+ylab("MATH")
dev.off()



## CNV                     ###################################
rm(list=ls())
load("lgg_tcga_risk_exp.RData")
coldata=data.frame(sample=rownames(tcga_risk_exp),Risk_Level=tcga_risk_exp$Risk_Level)

cnv <- read.delim("G:/fangshe/cell_death_qj3/cell_death_qj3/TCGA-LGG.gistic.tsv.gz",row.names = 1)
table(coldata$sample %in% colnames(cnv))
coldata$sample2=colnames(cnv)[match(coldata$sample,gsub("[.]","-",substr(colnames(cnv),1,15)))]
rownames(coldata)=coldata$sample2
coldata=coldata[coldata$sample2 %in% colnames(cnv),]
cnv=cnv[,coldata$sample2]
coldata=coldata[colnames(cnv),]

a=do.call(rbind,lapply(rownames(cnv),function(x){
  data.frame(gene=x,num=sum(abs(as.numeric(cnv[x,]))),p=chisq.test(table(as.character(cnv[x,]),coldata$Risk_Level))$p.value)
}))
a$fdr=p.adjust(a$p,method = "fdr",n=nrow(a))
table(a$p<0.01)
save(a,file = "cnv_de.RData")

load("G:/fangshe/cell_death_qj3/cell_death_qj3/gencode_annotation_gtf22.RData")
#load("lgg_tcga_risk_exp.RData")
table(gtf22$gene_id %in% a$gene)
b=merge(gtf22,a,by.x="gene_id",by.y="gene",all.y=T)
b=b[order(b$num,decreasing = T),]
b=b[!duplicated(b$gene_name),]
rownames(b)=b$gene_name

table(b$gene_type)
gene=b[1:15,"gene_name"]

cnv$gene_name=gtf22[match(rownames(cnv),gtf22$gene_id),"gene_name"]
cnv=cnv[cnv$gene_name %in% gene,]
rownames(cnv)=cnv$gene_name
cnv$gene_name=NULL

library(pheatmap)
library(gplots)
library(colorRamps)
library(RColorBrewer)
coldata$Risk_Level=factor(coldata$Risk_Level,levels = c("L","H"),labels = c("Low","High"))
coldata=coldata[order(coldata$Risk_Level),]
coldata$sample=NULL
coldata$sample2=NULL

ann_colors = list(rtCDI=c(High="#c74546", Low="#88c4e8"))
pdf("cnv_top15_heatmap.pdf",6,4)
pheatmap(cnv[gene,rownames(coldata)],color= c("#88c4e8","grey95","#c74546"),border_color = NA, 
         annotation_col = coldata,annotation_colors = ann_colors,
         scale = "none", cluster_row = T,cluster_cols = F,show_rownames=T,
         show_colnames=F,clustering_method = "ward.D2",clustering_distance_cols = "euclidean",legend = TRUE)
dev.off()

## risk PD1                #########################################################
rm(list=ls())
options(stringsAsFactors = F)


load("phen_sva_lgg.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
tpm_tcga =tpm_sva[,substr(colnames(tpm_sva),1,4) %in% "TCGA"]
exp=as.data.frame(t(tpm_tcga[c("PDCD1","CD274"),]))

load("lgg_tcga_risk_exp.RData")
exp$Risk_Score=tcga_risk_exp[rownames(exp),"Risk_Score"]
exp$Risk_Level=tcga_risk_exp[rownames(exp),"Risk_Level"]

names(exp)[3] <- "rtCDI"
pdf("risk_pd1_pdl1.pdf")
ggplot(data=exp[!is.na(exp$rtCDI),], aes(y=rtCDI, x=CD274))+
  geom_point()+ theme(legend.position = "none")+  
  stat_smooth(method="glm",se=T)+ylab("rtCDI")+xlab("CD274 Expression")+
  stat_cor(method = "pearson",aes(color="#c74546"),label.sep = "\n", p.digits = 2)+ 
  theme_bw()+theme(panel.spacing.y=unit(0, "lines"),legend.position = "none")+  theme(text=element_text(size=20)) #change font size of all text)

ggplot(data=exp[!is.na(exp$rtCDI),], aes(y=rtCDI, x=PDCD1))+
  geom_point()+ theme(legend.position = "none")+  
  stat_smooth(method="glm",se=T)+ylab("rtCDI")+xlab("PD-1 Expression")+
  stat_cor(method = "pearson",aes(color="#c74546"),label.sep = "\n", p.digits = 2)+ 
  theme_bw()+theme(panel.spacing.y=unit(0, "lines"),legend.position = "none")+  theme(text=element_text(size=20))

ggplot(exp[!is.na(exp$Risk_Level),], aes(Risk_Level, CD274)) +
  geom_boxplot(aes(fill=Risk_Level),color = 'black', fill = c("#88c4e8", "#c74546"))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=Risk_Level),hide.ns = T)+
  labs(title = "")+xlab(label = "")+ylab("CD274 Expression")+ theme(text=element_text(size=20))

ggplot(exp[!is.na(exp$Risk_Level),], aes(Risk_Level, PDCD1)) +
  geom_boxplot(aes(fill = Risk_Level), color = 'black', fill = c("#88c4e8", "#c74546")) + # 手动设置填充颜色
  geom_jitter(width = 0.3, alpha = 0.7, color = 'black') +
  stat_compare_means(aes(group = Risk_Level), hide.ns = TRUE) +
  labs(title = "")+xlab(label = "")+ylab("PDCD1 Expression") +
  theme(text = element_text(size = 20))
dev.off()
  

##基因相关性
load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
load("lgg_tcga_risk_exp.RData")
tpm_tcga =tpm_sva[,substr(colnames(tpm_sva),1,4) %in% "TCGA"]
tpm_tcga = as.data.frame(t(tpm_tcga))
tpm_tcga$Risk_Score = tcga_risk_exp[match(rownames(tpm_tcga),rownames(tcga_risk_exp)),"Risk_Score"]
tpm_tcga = na.omit(tpm_tcga)
tpm_tcga = as.data.frame(t(tpm_tcga))
singleGene_cor <- function(gene){
  y <- as.numeric(tpm_tcga[gene,])
  rownames <- rownames(tpm_tcga)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(tpm_tcga[x,]), y, type='pearson')
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value)
  }))
}

dd <- singleGene_cor('PDCD1') #如果更换基因，只用更新基因名即可。
save(dd,file = "PDCD1Immune_Checkpoint.Rdata")
library(dplyr)
library(tidyr)
corSig <- dd %>% 
  filter(p.value < 0.05) %>% 
  filter(abs(cor) < 0.5)
library(ggstatsplot)
exprSet=as.data.frame(t(tpm_tcga))
names(exprSet)[21921] = "rtCDI"
ggscatterstats(data = exprSet,
               y = PDCD1,
               x = rtCDI,
               #type = "pearson",
               centrality.para = "mean",
               margins = "both",
               xfill = "#88c4e8",
               yfill = "#c74546",
               marginal.type = "histogram",
               title = "Relationship between PDCD1 and rtCDI")

## risk ICB                ##################################
rm(list=ls())
options(stringsAsFactors = F)

gene_set_group <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/gene_set_group.txt", header=FALSE)
table(gene_set_group$V1)

geneset_gen = unique(gene_set_group[gene_set_group$V1=="ICG","V2"])
desired_order <- c("HLA-A","HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DPA1","HLA-DQB1","HLA-DRA","MICA","MICB","ICAM1","ITGB2","CD28","BTN3A1","BTN3A2","PDCD1LG2","SLAMF7","VTCN1","CD28","ICOSLG","CCL5","CX3CL1","CXCL10","CXCL9","IL12A","IL1B","TGFB1","TNF","TNFSF4","TNFSF9","VEGFA","VEGFB","ENTPD1","GZMA","HMGB1","PRF1","CD27","ADORA2A","CD40","PDCD1","TNFRSF14","TNFRSF4","EDNRB","LAG3","TLR4","ICOS","TNFRSF9","TNFRSF18")
# Stimulating = c("ICOS","CD134","CD27","CD137","CD40","CD278","TNFRSF9","CD134","TNFSF4","CD154","TNFSF5","DNAM1","2B4","DC-SIGN","DR3")
# Inhibitory = c("CTLA4","CD274","TIGIT","LMTK3","PD-1","LAG3","TIGIT","BTLA","CD279","CD152","CD223","TIM3","AA2R","CEACAM1","CD200R","CD28","GITR","LIGHT")
#table(geneset_gen %in% Inhibitory)
# load("G:/fangshe/cell_death_qj3/cell_death_qj3/phen_sva_lgg.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
# exp=as.data.frame(t(exp[rownames(exp) %in% Stimulating,]))
exp =tpm_sva[,substr(colnames(tpm_sva),1,4) %in% "TCGA"]
# exp=as.data.frame(t(exp[rownames(exp) %in% Inhibitory,]))
exp=as.data.frame(t(exp))
# 使用select()函数重新排列列的顺序
exp <- exp %>% select(desired_order)
load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
exp$Risk_Score=tcga_risk_exp[rownames(exp),"Risk_Score"]
exp$Risk_Level=tcga_risk_exp[rownames(exp),"Risk_Level"]
exp$Risk_Level=factor(exp$Risk_Level,levels = c("H","L"),labels = c("High","Low"))
exp=na.omit(exp)
exp1=reshape2::melt(exp,id.vars = c("Risk_Score","Risk_Level"))
colnames(exp1)
gen_sort=do.call(rbind,lapply(unique(exp1$variable), function(x){data.frame(tcga_type=x,gene=median(exp1[exp1$variable==x,"value"]))}))

gen_sort=na.omit(gen_sort)
#gen_sort=gen_sort[order(gen_sort$gene),]
exp1$Gene=factor(exp1$variable,levels = gen_sort$tcga_type)

names(exp1)[2] <- "rtCDI"
# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("High" = "#c74546", "Low" = "#88c4e8")
pdf("risk_ICB_boxplot.pdf",10,5)
# 设置 "High" 和 "Low" 的颜色
ggplot(exp1, aes(x = Gene, y = value, fill = rtCDI)) +
  geom_boxplot() +
  xlab(label = "") +
  stat_compare_means(aes(group = rtCDI), label = "p.signif", na.rm = TRUE, hide.ns = TRUE) +
  ylab("Gene expression") +
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(text=element_text(size=17))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

pdf("risk_ICBStimulating_boxplot.pdf",10,6)
# 设置 "High" 和 "Low" 的颜色
ggplot(exp1, aes(x = Gene, y = value, fill = rtCDI)) +
  geom_boxplot() +
  xlab(label = "") +
  stat_compare_means(aes(group = rtCDI), label = "p.signif", na.rm = TRUE, hide.ns = TRUE) +
  ylab("Gene expression") +
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
   theme(text=element_text(size=17))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

pdf("risk_ICBInhibitory_boxplot.pdf",10,6)
# 设置 "High" 和 "Low" 的颜色
ggplot(exp1, aes(x = Gene, y = value, fill = rtCDI)) +
  geom_boxplot() +
  xlab(label = "") +
  stat_compare_means(aes(group = rtCDI), label = "p.signif", na.rm = TRUE, hide.ns = TRUE) +
  ylab("Gene expression") +
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  theme(text=element_text(size=17))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

## risk IPS                ################################
rm(list=ls())
options(stringsAsFactors = F)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_ips_score.RData")
rownames(ips_score)=ips_score$sample
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_tcga_risk_exp.RData")

load("G:/fangshe/cell_death_qj2/death_pathway3/phen_lgg.RData")
ips_score$Risk_Score=tcga_risk_exp[ips_score$sample,"Risk_Score"]
ips_score$Risk_Level=tcga_risk_exp[ips_score$sample,"Risk_Level"]
ips_score$Risk_Level=factor(ips_score$Risk_Level,levels = c("H","L"),labels = c("High","Low"))
##您提供的这行代码正在为 ips_score 数据帧中名为 "Risk_Level "的列分配因子级别。它将 "Risk_Level "列中的值 "H "和 "L "转换为因子，并标注相应的 "高 "和 "低"。这样做通常是为了确保在创建图表或进行统计分析时能正确解释因子水平。
names(ips_score)[9] <- "rtCDI"
ips_score=ips_score[!is.na(ips_score$rtCDI),]

# 在这里，你可以替换 "blue"、"green" 和 "red" 为你喜欢的颜色
fill_colors <- c("High" = "#c74546", "Low" = "#88c4e8")
pdf("lgg_risk_ips.pdf",4,4)
# 设置 "High" 和 "Low" 的颜色
ggplot(ips_score, aes(rtCDI, IPS)) +
  geom_boxplot(aes(fill=rtCDI))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=rtCDI),hide.ns = T,label.y = 11)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("IPS score")
dev.off()

## risk TIDE               ################################
rm(list=ls())
options(stringsAsFactors = F)

#load("lgg_sva_risk_exp.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_tcga_risk_exp.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/phen_lgg.RData")

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_tide.RData")
lgg_tide =lgg_tide[substr(lgg_tide$Patient,1,4) %in% "TCGA",]
rownames(lgg_tide)=lgg_tide$Patient
table(rownames(lgg_tide) %in% substr(rownames(phen),1,15))
table(rownames(lgg_tide) %in% substr(rownames(tcga_risk_exp),1,15))

lgg_tide$Risk_Score=tcga_risk_exp[lgg_tide$Patient,"Risk_Score"]
lgg_tide$Risk_Level=tcga_risk_exp[lgg_tide$Patient,"Risk_Level"]
lgg_tide$Risk_Level=factor(lgg_tide$Risk_Level,levels = c("H","L"),labels = c("High","Low"))
lgg_tide=lgg_tide[!is.na(lgg_tide$Risk_Level),]
names(lgg_tide)[17] = "rtCDI"

fill_colors <- c("High" = "#c74546", "Low" = "#88c4e8")
pdf("lgg_risk_tide.pdf",4,4)
ggplot(lgg_tide, aes(rtCDI, TIDE)) +
  geom_boxplot(aes(fill=rtCDI))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=rtCDI),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("TIDE score")

ggplot(data=lgg_tide, aes(y=Risk_Score, x=TIDE))+
  geom_point(color="black")+ theme(legend.position = "none")+ stat_smooth(method="glm",se=T)+
  stat_cor(method = "pearson",aes(color = "#CA0020"),label.sep = "\n", p.digits = 2)
dev.off()

#载入所需R包：
library(stringr)
library(ggplot2)
library(patchwork)
#小提琴图展示结果：
#1.TIDE小提琴图：
my_comparisons <- list(c("High", "Low"))
pdf("TIDE_Dysfunction_Exclusion_MSI.pdf",4,4)
  #添加比较分组
p1 <- ggviolin(lgg_tide, x = 'rtCDI', y = 'TIDE', fill = 'rtCDI',
               palette = c("#c74546","#88c4e8"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p1
#method可选：t.test、wilcox.test、anova、kruskal.test
#2.Dysfunction小提琴图：
p2 <- ggviolin(lgg_tide, x = 'rtCDI', y = 'Dysfunction', fill = 'rtCDI',
               palette = c("#c74546","#88c4e8"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p2

#Exclusion小提琴图：
p3 <- ggviolin(lgg_tide, x = 'rtCDI', y = 'Exclusion', fill = 'rtCDI',
               palette = c("#c74546","#88c4e8"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p3

#MSI小提琴图：
colnames(lgg_tide)[6] <- c('MSI') #简化一下列名
p4 <- ggviolin(lgg_tide, x = 'rtCDI', y = 'MSI', fill = 'rtCDI',
               palette = c("#c74546","#88c4e8"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p4
dev.off()
## risk IC50               ################################
rm(list=ls())
options(stringsAsFactors = F)
##您提供的代码似乎是分析基因表达（或其他测量值）与癌症患者药物敏感性预测之间相关性的管道的一部分。它还能生成相关结果的图形表示。下面是代码的详细说明：
load("TCGA_oncoPredict_gdsc2.RData")##这个是啥数据
load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
dim(DrugPredictions_GDSC2)
DrugPredictions_GDSC2[1:4,1:4]
##qj
tcga_risk_exp$V1=rownames(tcga_risk_exp)
dat=merge(tcga_risk_exp[,c(11:13)],DrugPredictions_GDSC2,by="V1")
rownames(dat) = dat$V1
dat$V1=NULL
##您可以根据共同标识符合并 risk_exp 和 DrugPredictions_GDSC2 数据框中的相关列。

batch_cor <- function(gene){
  y <- as.numeric(dat[,gene])
  do.call(rbind,future_lapply(colnames(dat)[3:200], function(x){
    dd  <- cor.test(as.numeric(dat[,x]),y,type="spearman")
    data.frame(all_gene=x,cor=dd$estimate,p_value=dd$p.value )
  }))
}
library(future.apply)
#plan(multiprocess)
all_cor=batch_cor("Risk_Score")
all_cor$all_gene=as.character(all_cor$all_gene)
all_cor$fdr=p.adjust(p = all_cor$p_value,method = "fdr")
save(all_cor,file="risk_ic50_cor.RData")
##您定义了一个函数 batch_cor()，用于计算特定基因的表达与所有药物的药物敏感性预测之间的斯皮尔曼相关性和 p 值。您使用 future_lapply() 并行计算多个进程。

load("risk_ic50_cor.RData")
table(all_cor$fdr<0.01 & abs(all_cor$cor)>0.2)
table(all_cor$fdr<0.05)
all_cor=all_cor[all_cor$fdr<0.01 & abs(all_cor$cor)>0.2,]
colnames(all_cor)
all_cor=all_cor[order(all_cor$cor),]


all_cor <- subset(all_cor, row.names(all_cor) != "cor129")
all_cor$all_gene=str_split(all_cor$all_gene,"_",simplify = T)[,1]
all_cor$all_gene=factor(all_cor$all_gene,levels = all_cor$all_gene)

#修改列明
##加载之前保存的 risk_ic50_cor.RData 数据后，通过计算符合特定显著性和相关系数标准的基因来分析相关性。

pdf("lgg_sva_risk_ic50.pdf",12,6)
ggplot(data = all_cor, mapping = aes(x = all_gene, y = cor , fill = -log10(fdr))) +
  geom_bar(stat = 'identity',width = 0.7,position = position_dodge(0.9))+ labs(title = "",x="",y="Rs of drug sensitivity and rtCDI") +
  theme(axis.text.x = element_text(angle = 90, size=14,hjust = 1))+
   
  scale_fill_gradient2( low = "#88c4e8",mid = "#88c4e8", high = "#c74546",midpoint = 10)

dev.off()

##您可以使用 ggplot2 软件包创建 PDF 图，以直观显示相关性。x 轴代表基因，y 轴代表相关系数，颜色代表 FDR 的负对数。颜色梯度表示相关性的显著性。

load("risk_ic50_cor.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/drug_pathway_gdsc2.RData")
drug_pathway$all_gene=str_split(drug_pathway$all_gene,"_",simplify = T)[,1]
all_cor$all_gene=str_split(all_cor$all_gene,"_",simplify = T)[,1]

all_cor$all_gene=as.character(all_cor$all_gene)

all_cor$all_gene=gsub("[.]|-","_",all_cor$all_gene)
table(drug_pathway$all_gene %in% all_cor$all_gene)

dat=merge(all_cor,drug_pathway[,c(7,4)],by="all_gene",all=T)
table(dat$fdr<0.05 & abs(dat$cor) > 0.2)
dat=dat[dat$fdr<0.05 & abs(dat$cor) > 0.2,]
dat = na.omit(dat)
pdf("lgg_sva_risk_ic50_pathway.pdf",12,6)
ggplot(dat,aes(x=all_gene,y=pathway_name))+
  geom_point(aes(size=-log10(fdr),color=cor))+
  scale_color_gradient(high="#c74546",low="#88c4e8")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 13),axis.text.y = element_text(size = 13))+
  xlab("")+ylab("")
dev.off()

####
load("TCGA_oncoPredict_gdsc2.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/risk_ic50_cor.RData")
res = DrugPredictions_GDSC2
res$sample_group <- tcga_risk_exp$Risk_Level[match(res$V1, rownames(tcga_risk_exp))]
colnames(res)=str_split(colnames(res),"_",simplify = T)[,1]
res = na.omit(res)
# a = res[,c("BI-2536","Saracatinib","EX-527
# ","VX-680","Tivozanib","HG-5-113-01
# ","Midostaurin","GW-2580","XAV939","Paclitaxel","Salubrinal","AKT inhibitor VIII","T0901317")]
# b = res[,c("Linsitinib","LFM-A13","STF-62247","BMS-754807","Bryostatin-1","KIN001-266","BMS-509744","TAK-715","GSK1904529A","Cetuximab","MP470","TL-1.85","QS11")]
c("Linsitinib","BMS-754807","BI-2536","XAV939","V1","Vorinostat","Dasatinib") %in% colnames(res)

c = res[,c("Linsitinib","Entospletinib","Rapamycin","ERK","V1","Vorinostat","Dasatinib")]
rownames(c) = c$V1
c$V1 = NULL
# c$sample_group <- tcga_risk_exp$Risk_Level[match(rownames(c), rownames(tcga_risk_exp))]
c = na.omit(c)

# 定义样本组的颜色映射
sample_group_colors <- c("L" = "#88c4e8", "H" = "#c74546")
sample_group = tcga_risk_exp$Risk_Level

pdf("lgg_HL_ic50.pdf",10,5)
c %>%
  bind_cols(sample_group = sample_group) %>%
  pivot_longer(1:6, names_to = "drugs", values_to = "ic50") %>%
  ggplot(., aes(sample_group, ic50, fill = sample_group)) +
  geom_boxplot() +
  # 使用手动指定的颜色映射
  scale_fill_manual(values = sample_group_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  facet_wrap(vars(drugs), scales = "free_y", nrow = 2) +
  stat_compare_means()
dev.off()


## risk EMT                #################################
rm(list=ls())

load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
EMT_gene <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/EMT_gene.txt", header=FALSE)

library(limma)
tpm_log=scale(tpm_sva)
boxplot(tpm_log[,1:30])
tpm_log = tpm_log[,substr(colnames(tpm_log),1,4) %in% "TCGA"]
gen_e=unique(EMT_gene[EMT_gene$V2=="E",1])
gen_e=gen_e[gen_e %in% rownames(tpm_log)]
gen_m=unique(EMT_gene[EMT_gene$V2=="M",1])
gen_m=gen_m[gen_m %in% rownames(tpm_log)]
dat=as.data.frame(t(tpm_log[rownames(tpm_log) %in% EMT_gene$V1,]))
dat$emt_score=rowSums(dat[,gen_m])/length(gen_m)-rowSums(dat[,gen_e])/length(gen_e)

load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_tcga_risk_exp.RData")
dat$Risk_Score=tcga_risk_exp[rownames(dat),"Risk_Score"]
dat$Risk_Level=tcga_risk_exp[rownames(dat),"Risk_Level"]
dat$Risk_Level=factor(dat$Risk_Level,levels = c("H","L"),labels = c("High","Low"))
dat=na.omit(dat)
save(dat,file="cluster_risk_emt.RData")

names(dat)[78] = "rtCDI"
fill_colors <- c("High" = "#c74546", "Low" = "#88c4e8")
pdf("emt_risk.pdf",4,4)
ggplot(data=dat, aes(y=emt_score, x=Risk_Score))+
  geom_point()+ theme(legend.position = "none")+  
  stat_smooth(method="glm",se=T)+ylab("EMT score")+xlab("Risk_Score")+
  stat_cor(method = "pearson",aes(color="#CA0020"),label.sep = "\n", p.digits = 2)+ 
  theme_bw()+theme(panel.spacing.y=unit(0, "lines"),legend.position = "none")

ggplot(dat, aes(rtCDI, emt_score)) +
  geom_boxplot(aes(fill=rtCDI))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=rtCDI),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("EMT score")
dev.off()

## risk Estimate###############
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)

load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
colnames(tpm_sva) = gsub("-",".",colnames(tpm_sva))
tpm_TCGA = tpm_sva[,colnames(tpm_sva) %in% rownames(tcga_risk_exp)]
dim(tpm_TCGA)#[1] 21920   279
Grouplist = ifelse(as.numeric(substr(colnames(tpm_TCGA),14,15)) < 10,'tumor','normal')
#normal  0 tumor  534
estimate <- function(dat,pro){
  input.f <-  paste0(pro,'_estimate_input.txt')
  output.f <- paste0(pro,'_estimate_gene.gct')
  output.ds <- paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  ###转换表达矩阵格式
  filterCommonGenes(input.f=input.f,output.f=output.f ,id="GeneSymbol")
  ###计算Score
  estimateScore(input.ds = output.f,output.ds=output.ds,platform="illumina") ## 注意platform
  scores <- read.table(output.ds,skip = 2,header = T)
  rownames(scores) <- scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  return(scores)
}


TPMscores = estimate(tpm_TCGA,pro = "LGG_TPM")
head(TPMscores)
#[1] "Merged dataset includes 9911 genes (501 mismatched)."
#[1] "1 gene set: StromalSignature  overlap= 136"
#[1] "2 gene set: ImmuneSignature  overlap= 140"
TPMscores = as.data.frame(TPMscores)
save(TPMscores,file = "Estimate_TCGA_scores.Rdata")

load("G:/fangshe/cell_death_qj3/cell_death_qj3/Estimate_TCGA_scores.Rdata")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
colnames(TPMscores)
rownames(tcga_risk_exp) = gsub("-",".",rownames(tcga_risk_exp))
TPMscores$group = tcga_risk_exp[rownames(TPMscores) %in% rownames(tcga_risk_exp), "Risk_Level"]
TPMscores$sample <- rownames(TPMscores)
TPMscores =  melt(TPMscores)

## Using group, sample as id variables

colnames(TPMscores)=c("group","sample","status","score")  #设置行名
head(TPMscores)

##      group          sample       status      score
## 1    Tumor TCGA.CV.6943.01 StromalScore   906.2923
## 2    Tumor TCGA.CV.6959.01 StromalScore  -352.6656
## 3 Nontumor TCGA.CV.7438.11 StromalScore -1183.4705
## 4 Nontumor TCGA.CV.7242.11 StromalScore -1067.1461
## 5    Tumor TCGA.CV.7432.01 StromalScore -1234.5253
## 6 Nontumor TCGA.CV.6939.11 StromalScore   424.8381

# 5.3 计算误差线
ESTI_Data_summary <- summarySE(TPMscores, measurevar="score", groupvars=c("group","status"))
head(ESTI_Data_summary)

##      group        status  N       score        sd        se       ci
## 1 Nontumor  StromalScore 43  -918.91119  837.0492 127.64881 257.6057
## 2 Nontumor   ImmuneScore 43  -210.02121  504.0363  76.86481 155.1195
## 3 Nontumor ESTIMATEScore 43 -1128.93240 1138.1268 173.56271 350.2637
## 4    Tumor  StromalScore 43  -517.35577  659.5859 100.58590 202.9906
## 5    Tumor   ImmuneScore 43   -67.37634  638.1790  97.32139 196.4025
## 6    Tumor ESTIMATEScore 43  -584.73211 1138.8398 173.67145 350.4832
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}

ESTI_split_violin <- ggplot(TPMscores,aes(x= status,y= score,fill= group))+
  geom_split_violin(trim= F,color="white",scale = "area") + #绘制分半的小提琴图
  geom_point(data = ESTI_Data_summary,aes(x= status, y= score),pch=19,
             position=position_dodge(0.4),size= 1)+ #绘制均值为点图
  geom_errorbar(data = ESTI_Data_summary,aes(ymin = score-ci, ymax= score+ci), 
                width= 0.05, 
                position= position_dodge(0.4), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c("#c74546", "#88c4e8"))+ 
  labs(y=(""),x=NULL,title = "") + 
  theme_bw()+ mytheme +
  scale_x_discrete(labels=c("Stromal","Immune","ESTIMATE")) +
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox",
                     label.y = max(TPMscores$score),
                     hide.ns = T)
ESTI_split_violin

ggsave(ESTI_split_violin,filename = "./Output/ESTIMATE_plot.pdf", height = 10,width = 10,units = "cm")




## risk mRNAsi             ##################################
rm(list=ls())
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_gtex_all_tumor_mrnasi.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_tcga_risk_exp.RData")
mrnasi$Risk_Score=tcga_risk_exp[rownames(mrnasi),"Risk_Score"]
mrnasi$Risk_Level=tcga_risk_exp[rownames(mrnasi),"Risk_Level"]
mrnasi$Risk_Level=factor(mrnasi$Risk_Level,levels = c("H","L"),labels = c("High","Low"))
mrnasi=mrnasi[!is.na(mrnasi$Risk_Score),]
save(mrnasi,file="cluster_risk_mrnasi.RData")

mrnasi=mrnasi[,c(4,8,12:13)]
colnames(mrnasi)[1]="mRNAsi"
colnames(mrnasi)[2]="mDNAsi"
colnames(mrnasi)[4]="rtCDI"
fill_colors <- c("High" = "#c74546", "Low" = "#88c4e8")

pdf("mRNAsi_risk.pdf",4,4)
ggplot(data=mrnasi, aes(y=mRNAsi, x=Risk_Score))+
  geom_point()+ theme(legend.position = "none")+  
  stat_smooth(method="glm",se=T)+ylab("mRNAsi score")+xlab("Risk_Score")+
  stat_cor(method = "pearson",aes(color="#CA0020"),label.sep = "\n", p.digits = 2)+ 
  theme_bw()+theme(panel.spacing.y=unit(0, "lines"),legend.position = "none")

ggplot(mrnasi, aes(rtCDI, mRNAsi)) +
  geom_boxplot(aes(fill=rtCDI))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=rtCDI),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("mRNAsi score")
dev.off()

pdf("mDNAsi_risk.pdf",4,4)
ggplot(data=mrnasi, aes(y=mDNAsi, x=Risk_Score))+
  geom_point()+ theme(legend.position = "none")+  
  stat_smooth(method="glm",se=T)+ylab("mDNAsi score")+xlab("Risk_Score")+
  stat_cor(method = "pearson",aes(color="#CA0020"),label.sep = "\n", p.digits = 2)+ 
  theme_bw()+theme(panel.spacing.y=unit(0, "lines"),legend.position = "none")
ggplot(mrnasi, aes(rtCDI, mDNAsi)) +
  geom_boxplot(aes(fill=rtCDI))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=rtCDI),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("mDNAsi score")
dev.off()


## stratum                 ###########################################

rm(list=ls())
load("phen_lgg.RData")
phen$Age=ifelse(phen$Age>=60,">=60","<60")
colnames(phen)
phen_sva=unique(phen[,c(8,10,3,2,6,7)])

library(reshape2)
colnames(phen)
phen$link=1
risk_exp1=melt(phen,id.vars = "link")
variable <- summary(risk_exp1$variable)
risk_exp1$flow <- rep(1:variable[1], length(variable))

mycol <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462','#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8',
           '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5',
           '#6181BD', '#F34800', '#64A10E', '#FF00FF', '#c7475b', '#049a0b', '#BEAED4',
           '#FDC086', '#FFFF99', '#386CB0', '#F0027F', '#4253ff', '#ff4308', '#D8D155','#64495D', '#7CC767')

library(ggalluvial)
pdf("stratum_plot.pdf",16,8)
ggplot(risk_exp1, aes(x = variable, y = 1,stratum = value, alluvium = flow, fill = value)) +
  geom_stratum() + geom_flow(aes.flow = 'forward') + scale_fill_manual(values = mycol) + 
  geom_text(stat = 'stratum', infer.label = TRUE, size = 2.5) +
  labs(x = '', y = '') + theme(legend.position = 'none', panel.background = element_blank(),line = element_blank(), axis.text.y = element_blank())
dev.off()

## risk score pheno        #######################################
rm(list=ls())

load("phen_lgg.RData")
phen$Age=ifelse(phen$Age>=60,">=60","<60")
phen$ldh1_mutation=ifelse(phen$ldh1_mutation %in% "0","non_Mut",ifelse(phen$ldh1_mutation %in% "1","Mut",NA))


pdf("lgg_sva_riskscore_pheno.pdf",4,4)
fill_colors <- c("FEMALE" = "#c74546", "MALE" = "#88c4e8")
ggplot(phen[!is.na(phen$Gender),], aes(Gender, Risk_Score))+
  geom_boxplot(aes(fill=Gender))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=Gender),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("Risk_Score")

fill_colors <- c(">=60" = "#c74546", "<60" = "#88c4e8")
ggplot(phen[!is.na(phen$Age),], aes(Age, Risk_Score))+
  geom_boxplot(aes(fill=Age))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=Age),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("Risk_Score")

fill_colors <- c("G2" = "#c74546", "G3" = "#88c4e8")
ggplot(phen[!is.na(phen$Grade) & phen$Grade != "",], aes(Grade, Risk_Score)) +
  geom_boxplot(aes(fill=Grade))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=Grade),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("Risk_Score")

fill_colors <- c("non_Mut" = "#c74546", "Mut" = "#88c4e8")
ggplot(phen[!is.na(phen$ldh1_mutation),], aes(ldh1_mutation, Risk_Score)) +
  geom_boxplot(aes(fill=ldh1_mutation))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
  stat_compare_means(aes(group=ldh1_mutation),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = fill_colors) +  # 手动设置填充颜色
  labs(title = "")+xlab(label = "")+ylab("Risk_Score")

dev.off()


## cibersotrX反卷积可视化#####################
CIBERSORTx_no_Results <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/CIBERSORTx_no_Results.txt", row.names=1)
CIBERSORTx_yes_Results <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/CIBERSORTx_yes_Results.txt", row.names=1)
CGGA_matrix_no <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/CGGA_matrix_no.txt", row.names=1)
CGGA_matrix_yes <- read.delim("G:/fangshe/cell_death_qj2/death_pathway3/CGGA_matrix_yes.txt", row.names=1)

res = rbind(CIBERSORTx_no_Results,CIBERSORTx_yes_Results)
dt = res[,-(11:13)]

#宽数据转换为长数据(ggplot2绘图格式)：
dt <- dt %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type, value = value, -sample)
head(dt)

#将value根据样本转换为百分比形式(新增一列)：
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)

#指定绘图顺序：
dtt$cell_type <- factor(dtt$cell_type,levels = unique(dt$cell_type))

#自定义主题：
mytheme <- theme(axis.title = element_text(size = 12),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.text = element_text(size = 10),
                 legend.position = "bottom")

#配色挑选
d_palettes <- palettes_d_names #查看离散型色板(paletteer包)
col <- paletteer_d("khroma::smoothrainbow",n=28)

#col <- paletteer_d("RColorBrewer::Paired")

p <- ggplot(dtt,
            aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = col) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p

#重新指定箱线图排序（按相对丰度中位数从大到小）：
dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)

dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))


#重新绘图(代码相同)：
p1 <- ggplot(dtt,
             aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = col) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p1

#新增分组列：

dtt$group =ifelse(dtt$sample %in% rownames(CIBERSORTx_yes_Results),"radiation",ifelse(dtt$sample %in% rownames(CIBERSORTx_no_Results),"NO_radiation",NA))

dtt = na.omit(dtt)

#分组箱线图展示：
p2 <- ggplot(dtt,
             aes(x = cell_type,y = proportion,fill = group)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = c("#4979b6","#d9352a")) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme + theme(axis.text.x = element_text(angle = 45))
p2

# Assuming you have a data frame called dtt with columns cell_type, sample, and proportion


#使用t test或wilcox test进行两两比较(T检验为例)：
t <- t.test(proportion ~ group, data = dtt)
tj <- adjust_pvalue(t, method = 'fdr') #p值矫正；
tj

#根据p.adj添加显著性标记符号；
tj <- add_significance(tj, 'p.value')
tj
#在图表中添加 p 值或者显著性标记；
lab <- add_xy_position(tj, x = 'cell_type', dodge = 0.65)
pdf("LGGGBM_ssgsea_lq.pdf",13,9)

###免疫细胞相关性分析
#计算相关性系数：
res = rbind(CIBERSORTx_no_Results,CIBERSORTx_yes_Results)
res = res[,-(11:14)]
rescor <- cor(res)
View(rescor)
library(corrplot)
corrplot(rescor,
         method = "square",
         order = "hclust",
         tl.cex = 0.6,
         tl.col = "black")
#添加显著性标记:
resorp <- cor.mtest(rescor, conf.level = .95) #使用cor.mtest做显著性检验;

#提取p值矩阵；
p.mat <- resorp$p
View(p.mat)

#相关性热图中展示显著性标记：
corrplot(rescor,
         method = "color",
         order = "hclust",
         tl.cex = 0.6,
         tl.col = "black",
         p.mat = resorp$p, sig.level = c(.001, .01, .05),outline="white",
         insig = "label_sig",pch.cex = 0.6, pch.col = "white")

##兴趣基因与各类免疫细胞的相关性分析
tpm = cbind(CGGA_matrix_no,CGGA_matrix_yes)
tpm1=log(tpm+1)
#"SFRP2"    "MELK"     "ATP6V1G2"
#构建目的基因集：
genes <- c("SFRP2","MELK","ATP6V1G2")
#从原表达矩阵中取子集：
exp_genes <- tpm1[genes,]

#将免疫细胞丰度矩阵和目的基因表达量矩阵合并：
rb <- cbind(res,t(exp_genes))
rownames(rb)

#重新计算相关性系数（步骤同上）：
rbcor <- cor(rb)

#计算显著性差异（步骤同上）：
rbcorp <- cor.mtest(rbcor, conf.level = .95)

#提取p值矩阵（步骤同上）：
p.mat2 <- rbcorp$p
res = t(res)
#切割相关性矩阵：
split <- rbcor[1:nrow(res), #行取免疫细胞所在的行
               (ncol(rb)-length(genes)+1):ncol(rb)] #列取目的基因所在的行
View(split) #行名免疫细胞，列名目的基因；

#切割p值矩阵:
splitp <- p.mat2[1:nrow(res), #行取免疫细胞所在的行
                 (ncol(rb)-length(genes)+1):ncol(rb)] #列取目的基因所在的行
View(splitp)

# Check for missing values in splitp
missing_values <- is.na(splitp)
# Create a vector of marks with a default value of "NA"
marks <- rep("NA", length(splitp))
# Set marks based on non-missing values
marks[!missing_values & splitp < 0.05] <- "***"
marks[!missing_values & splitp < 0.1] <- "**"
marks[!missing_values & splitp < 0.5] <- "*"

# Convert the marks vector into a matrix
mark <- matrix(marks, nrow = nrow(splitp))


# Create the mark matrix with default marks for missing values
mark <- matrix(case_when(
  splitp < 0.05 ~ "***",
  splitp < 0.1 ~ "**",
  splitp < 0.5 ~ "*",
  TRUE ~ "NA"),  # Default mark for missing values
  nrow = nrow(splitp))

#免疫细胞-基因相关性热图：
col2<-colorRampPalette(c("blue","white", "red"))(95)

library(pheatmap)
pheatmap(t(split),
         display_numbers = t(mark),
         number_color = "black",
         fontsize_number = 13,
         color = col2,
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_col = 8,
         fontsize_row = 8,
         angle_col = 90)

## 火山图     #########################################################

load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_rr_rs_de.RData")
res$logFC=-res$logFC
p_v=0.05   
cor_v=0.5 
de_gene=rownames(res[res$adj.P.Val< 0.05 & abs(res$logFC)>0.5,])

#画火山图
df = res
df$group = 'ns' #添加一组

df$group[which((res$adj.P.Val<0.05)&(res$logFC > 0.5))]='Up'#上调基因
df$group[which((res$adj.P.Val<0.05)&(res$logFC< -0.5))]='Down'

write.table(df,"tcga_lgg_rr_rs_de.csv",sep=",")
#下调基因
table(df$group)
df$padj = -log10(res$adj.P.Val)#将P转化为-Logp

# 将需要标记的基因放置在label列(logFC >= 5)
library(ggrepel)
df$label=""

gene <- c(rownames(multi_cox_coefficients))

df$label[match(gene,rownames(df))]<- gene #将up,down基因名添加在Label这一列
table(df$label)
p <- ggplot(
  # 数据、映射、颜色
  df, aes(x = logFC, y = -log10(adj.P.Val), colour = group)) +
  geom_point(alpha = 0.9, size = 1) +
  scale_color_manual(values = c("#4d97cd", "#bdc3d2", "#c74546")) +
  # 辅助线
  geom_vline(xintercept = c(-0.5, 0.5), lty = 5, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.01), lty = 5, col = "black", lwd = 0.8) +
  # 添加基因标签
  geom_text_repel(aes(label = df$label),size=4,direction="both",segment.alpha=0.6,max.overlaps =200,nudge_x = 0.2,nudge_y=0.2) +
  # 坐标轴
  labs(x = "log2(fold change)",
       y = "-log10 (adj.P.Val)") +
  theme_bw() +
  theme(text=element_text(size=20))+
  # 图例
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank())
p



#####veen图#####
rm(list=ls())
load("G:/fangshe/cell_death_qj2/death_pathway3/gene_list.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/tcga_lgg_rr_rs_de.RData")
#相关R包安装与载入：
library(VennDiagram)
library(eulerr)
install.packages("UpSetR")
library(UpSetR)
library(RColorBrewer)

p_v=0.05   
km_p=0.01   
cor_v=0.5 

All_diffSig = res[res$adj.P.Val<p_v & abs(res$logFC)>cor_v,]
all_diffsig = rownames(All_diffSig)

#指定统计的分组列，并设置作图颜色、字体样式等
venn_list <- list(DEGs= rownames(All_diffSig),PCD = gene_list$Gene)

pdf("lgg_veen.pdf")
venn.plot = venn.diagram(venn_list,
                         ext.text = F,
                         resolution = 20, 
                         imagetype = "pdf", 
                         alpha=c(0.7,0.7),
                         fill=c("#c74546","#88c4e8"),
                         cat.fontface = 2,
                         fontfamily = "serif",
                         cex = 2,
                         cat.cex = 2, # 字体大小
                         main="",
                         main.cex = 2, 
                         main.fontface = 2, 
                         main.fontfamily = 2,
                         filename = NULL)

# 显示 Venn 图
grid.draw(venn.plot)

grid.draw(venn.plot)
dev.off() # 关闭当前的绘图设备
grid.newpage() # 清除图形

file.remove(list.files(pattern = "*log"))
dev.off()
#取交集
inter <- get.venn.partitions(venn_list)
interset_mRNA <- as.data.frame(inter$..values..[1])
All_diffSig_veen =  All_diffSig[interset_mRNA$X1,]##提取交集第四列第一行
save(inter,file = "venninter.Rdata")
save(interset_mRNA,file = "interset_gene.Rdata")
write.csv(interset_mRNA,file = "36rtPCDrelated_genes.csv")

#####热图####
rm(list=ls())
library(pheatmap)
library(gplots)
library(colorRamps)
library(RColorBrewer)
library(ComplexHeatmap)


load("phen_lgg.RData")
rownames(phen) <- substr(phen$Sample, 1, 15)
phen$Age=ifelse(phen$Age>=60,">=60",ifelse(phen$Age<60,"<60",NA))
phen$ldh1_mutation=ifelse(phen$ldh1_mutation %in% "0","non_Mut",ifelse(phen$ldh1_mutation %in% "1","Mut",NA))
colnames(phen)[14] = 'rtCDI'


load("G:/fangshe/cell_death_qj2/death_pathway3/tpm_sva_lgg_tpm.RData")
table(colnames(tpm_sva) %in% rownames(phen))

load("multi_cox_coefficients.RData")
multi_cox_gene=rownames(multi_cox_coefficients)
exp=tpm_sva[multi_cox_gene,colnames(tpm_sva) %in% rownames(phen)]

phen=phen[colnames(exp),c(8,14,2,3,5)]
phen=phen[order(phen$rtCDI),]
colnames(phen)

ann_colors = list(
  rtCDI=c( H ="#c74546", L="#88c4e8"),
  dataset=c(CGGA325="#8DD3C7", CGGA693="#FB8072",TCGA="#BEBADA"),
  Gender=c(FEMALE="#FB9A99",MALE="#A65628"),
  Age=c(`>=60`="#4DAF4A",`<60`="#984EA3"),
  ldh1_mutation=c(non_Mut="#0571B0",Mut="#CA0020"),
  Grade=c(G2="#FFFFB3", G3="#ffad60")
)
pdf("model_heatmap.pdf",7,5)
ComplexHeatmap::pheatmap(
  as.matrix(exp[,rownames(phen)]),
                 scale = "row",
                 cluster_cols = T,
                 show_rownames = T,
                 show_colnames = F,
                 annotation_colors = ann_colors,
                 col = colorRampPalette(c("navy","white","firebrick3"))(100),
  row_labels = rownames(exp),
  annotation_col = phen,
  column_split=phen$rtCDI,
  column_title = NULL,row_title = NULL)

dev.off()


##感兴趣基因与TMB相关性计算及可视化############
load("D:/R/maftools/paad_tmb.Rdata")

load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
geneName <- c("Your Gene")  #select your interesting gene or genes
geneExp <- as.data.frame(t(paad.tpm[geneName,]))
#geneExp<- log2(geneExp) #if need log2()
dim(geneExp)
paad_tmb$Tumor_Sample_Barcode
comSample <- intersect(paad_tmb$Tumor_Sample_Barcode,rownames(geneExp))
length(comSample)

geneExp <- geneExp %>%
  rownames_to_column(var = "Tumor_Sample_Barcode") %>%
  inner_join(paad_tmb,by="Tumor_Sample_Barcode")
## or___##
df <- inner_join(geneExp,paad_tmb,by="Tumor_Sample_Barcode")
df <- merge(geneExp,paad_tmb,by="Tumor_Sample_Barcode")
##
cor(geneExp$Your Gene,geneExp$TMB)
cor.test(geneExp$Your Gene,geneExp$TMB)
str(geneExp)
geneExp <- as.data.frame(geneExp)




##多种免疫细胞浸润分析########

 load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_timer_score.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_xcell.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_MCPcounterScore.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_quantiseq_score.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_cib_raw_result.RData")
load("G:/fangshe/cell_death_qj2/death_pathway3/lgg_sva_epic_score.RData")
load("G:/fangshe/cell_death_qj3/cell_death_qj3/lgg_tcga_risk_exp.RData")
#加载需要的R包
library(pheatmap)
library(RColorBrewer)
#library(immunedeconv)
library(tidyverse)
library(gplots)
#添加免疫浸润算法名称
xcell_result =xcell_result[,substr(colnames(xcell_result),1,4) %in% "TCGA"]
xcell_result$cell_type = rownames(xcell_result)
xcell_result$cell_type<-paste0(xcell_result$cell_type,"_XCELL")
#添加新的列名
xcell_result$Immunemethod<-"XCELL"

epic_score=as.data.frame(t(epic_score$cellFractions))
epic_score =epic_score[,substr(colnames(epic_score),1,4) %in% "TCGA"]
#添加免疫浸润算法名称
epic_score$cell_type = rownames(epic_score)
epic_score$cell_type<-paste0(epic_score$cell_type,"_EPIC")
#添加新的列名
epic_score$Immunemethod<-"EPIC"

MCPcounterScore=as.data.frame(MCPcounterScore)
MCPcounterScore =MCPcounterScore[,substr(colnames(MCPcounterScore),1,4) %in% "TCGA"]
#添加免疫浸润算法名称
MCPcounterScore$cell_type = rownames(MCPcounterScore)
MCPcounterScore$cell_type<-paste0(MCPcounterScore$cell_type,"_MCPcounter")
#添加新的列名
MCPcounterScore$Immunemethod<-"MCPcounter"


cib_result =cib_result[,substr(colnames(cib_result),1,4) %in% "TCGA"]
#添加免疫浸润算法名称
cib_result$cell_type = rownames(cib_result)
cib_result$cell_type<-paste0(cib_result$cell_type,"_CIBERSORT")
#添加新的列名
cib_result$Immunemethod<-"CIBERSORT"


timer_score =timer_score[,substr(colnames(timer_score),1,4) %in% "TCGA"]
#添加免疫浸润算法名称
timer_score$cell_type = rownames(timer_score)
timer_score$cell_type<-paste0(timer_score$cell_type,"_TIMER")
#添加新的列名
timer_score$Immunemethod<-"TIMER"



a = quantiseq_score$cell_type
quantiseq_score =quantiseq_score[,substr(colnames(quantiseq_score),1,4) %in% "TCGA"]
quantiseq_score$cell_type = a
#添加免疫浸润算法名称
quantiseq_score$cell_type<-paste0(quantiseq_score$cell_type,"_QUANTISEQ")
#添加新的列名
quantiseq_score$Immunemethod<-"QUANTISEQ"

#合并六种算法免疫细胞矩阵
Immune<-rbind(xcell_result,epic_score,MCPcounterScore,cib_result,quantiseq_score,timer_score)

#提取免疫细胞矩阵
final<-select(Immune,-Immunemethod)
#把列名改成行名
rownames(final) <- NULL
final<-column_to_rownames(final,var="cell_type")
#匹配高低风险组样本信息
Immunesample<-match(rownames(tcga_risk_exp),colnames(final))
final<-final[,Immunesample]
# 定义颜色
#3.绘制免疫细胞热图
methods.col <- brewer.pal(n = length(unique(Immune$Immunemethod)),name = "Paired")
# 创建注释
# 列注释，位于热图顶端
annCol <- data.frame(RiskType = tcga_risk_exp$Risk_Level,
                     # 以上是risk score和risk type两种注释，可以按照这样的格式继续添加更多种类的注释信息，记得在下面的annColors里设置颜色
                     row.names = rownames(tcga_risk_exp),
                     stringsAsFactors = F)

annCol$RiskType=factor(annCol$RiskType,levels = c("H","L"),labels = c("H","L"))
# 行注释，位于热图左侧
annRow<-data.frame(Methods=factor(Immune$Immunemethod,levels=unique(Immune$Immunemethod)),
                   row.names = Immune$cell_type,
                   stringsAsFactors = F)

# 为各注释信息设置颜色
annColors <- list(Methods = c("XCELL" = methods.col[1],
                              "EPIC" = methods.col[2],
                              "MCPcounter" = methods.col[3],
                              "QUANTISEQ" =methods.col[4],
                              "TIMER"=methods.col[5],
                              "CIBERSORT"=methods.col[6]),
                   
                  # 下面是列注释的颜色，可依此设置更多注释的颜色
                  "RiskType" = c("H" = "#c74546","L" = "#88c4e8"))

# 数据标准化
indata <- final[,colSums(final) > 0] # 确保没有富集全为0的细胞
#确定免疫细胞含量范围
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
plotdata <- standarize.fun(indata,halfwidth = 2)
order <- order(annCol$RiskType)
sorted_plotdata <- plotdata[, order]

# 创建一个空的数据框来存储差异分析结果
differential_results <- data.frame(CellType = rownames(sorted_plotdata), p_value = numeric(nrow(sorted_plotdata)))

# 循环逐行计算t检验并存储结果
for (i in 1:nrow(sorted_plotdata)) {
  group1 <- sorted_plotdata[i, annCol$RiskType == "H"]
  group2 <- sorted_plotdata[i, annCol$RiskType == "L"]
  t_test_result <- wilcox.test(group1, group2)
  
  # 存储p值和细胞类型名称
  differential_results[i, "p_value"] <- t_test_result$p.value
}

# 根据显著性水平（例如，0.05）筛选具有显著性差异的细胞类型
significantly_different_celltypes <- differential_results[differential_results$p_value < 0.05,]

# pheatmap绘图
pdf("immune_heatmap_by_pheatmap.pdf", 10,20)

pheatmap::pheatmap(mat = as.matrix(sorted_plotdata), 
                   border_color = NA,
                   color = bluered(64),
                   cluster_rows = F,
                   cluster_cols = F,
                   show_rownames = T,
                   show_colnames = F,
                   annotation_col = annCol,
                   annotation_row = annRow,
                   annotation_colors = annColors,
                   gaps_row = cumsum(table(annRow$Methods)),
                   gaps_col = table(annCol$RiskType),
                   cellwidth = 0.8,
                   cellheight = 10)

dev.off()



