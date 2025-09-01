setwd("C:/Users/css/Desktop/标志物建模分析/FWBC")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(limma)
library(imputeLCMD)
library(pheatmap)
library(ggrepel)
library(survival)
library(survminer)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(DT)
library(tidyr)
library(tibble)
library(stringr)
library(limma)
library(GSVA)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggvenn)

volcano.plot=function(data.diff,fdr=0.01,fc=1.5,topn="TP53",base_family='sans',title =title){
  data.diff.s=merge(data.frame(Var1=c("down-regulated","unchanged","up-regulated"),color=c("steelblue","grey","red")),
                    as.data.frame(table(data.diff$Regulation)),by="Var1",all.x=T) %>% mutate(label=paste(Var1,Freq,sep=": "))
  data.diff$label=replace(data.diff$`Gene name`,!(data.diff$`Gene name` %in% topn),values = NA)
  p2=ggplot(data.diff,aes(x=log2FC,y=-log10(adj.P.Val),color=Regulation))+
    geom_point(size=1)+theme_classic(base_size = 14,base_family=base_family)+
    scale_color_manual(values = setNames(data.diff.s$color,data.diff.s$Var1),labels=data.diff.s$label)+
    geom_hline(yintercept = -log10(fdr),linetype=2)+geom_vline(xintercept = c(-log2(fc),log2(fc)),linetype=2)+
    geom_text_repel(aes(label=label),max.overlaps = 50) + labs(title = title)
  return(p2)
}
volcano.plot1=function(data.diff,pval=0.01,fc=1.5,topn="TP53",base_family='sans',title =title){
  data.diff.s=merge(data.frame(Var1=c("down-regulated","unchanged","up-regulated"),color=c("steelblue","grey","red")),
                    as.data.frame(table(data.diff$Regulation)),by="Var1",all.x=T) %>% mutate(label=paste(Var1,Freq,sep=": "))
  data.diff$label=replace(data.diff$`Gene name`,!(data.diff$`Gene name` %in% topn),values = NA)
  p2=ggplot(data.diff,aes(x=log2FC,y=-log10(P.Value),color=Regulation))+
    geom_point(size=1)+theme_classic(base_size = 14,base_family=base_family)+
    scale_color_manual(values = setNames(data.diff.s$color,data.diff.s$Var1),labels=data.diff.s$label)+
    geom_hline(yintercept = -log10(pval),linetype=2)+geom_vline(xintercept = c(-log2(fc),log2(fc)),linetype=2)+
    geom_text_repel(aes(label=label),max.overlaps = 50) + labs(title = title)
  return(p2)
}
volcano.plot.t=function(data.diff,fdr=0.01,fc=1.5,topn="TP53",base_family='sans', title = title){
  data.diff.s=merge(data.frame(Var1=c("down-regulated","unchanged","up-regulated"),color=c("steelblue","grey","red")),
                    as.data.frame(table(data.diff$Regulation)),by="Var1",all.x=T) %>% mutate(label=paste(Var1,Freq,sep=": "))
  data.diff$label=replace(data.diff$`Gene name`,!(data.diff$`Gene name` %in% topn),values = NA)
  p2=ggplot(data.diff,aes(x=log2FC,y=-log10(adj.P.Val),color=Regulation))+
    geom_point(size=1)+theme_classic(base_size = 14,base_family=base_family)+
    scale_color_manual(values = setNames(data.diff.s$color,data.diff.s$Var1),labels=data.diff.s$label)+
    geom_hline(yintercept = -log10(fdr),linetype=2)+geom_vline(xintercept = c(-log2(fc),log2(fc)),linetype=2)+
    geom_point(data = data.diff %>% filter(`Gene name` %in% topn),
               aes(log2FC,-log10(adj.P.Val)),color = "gold",size = 2) + 
    geom_text_repel(aes(label=label),max.overlaps = 10, color = "black") + labs(title = title)
  return(p2)
}
###------------------- TCGA & CPTAC data exploration ------------------###
#### TCGA #### 
projects <- getGDCprojects()$project_id
results <- GDCquery(project = "TCGA-BRCA",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
# GDCdownload(results,method = "api")
# expdat <- GDCprepare(query = results)
# save(expdat,file = "TCGA/TCGA-BRCA.Rdata")
load("TCGA/TCGA-BRCA.Rdata")
rowdata <- rowData(expdat)
metadata <- colData(expdat)
mrna <- expdat[rowdata$gene_type == "protein_coding"]
count <- assay(mrna,"unstranded")
gene_name <- rowData(mrna)$gene_name
# count.matrix <- merge(as.data.frame(rowdata[,c("gene_id","gene_name")]),as.data.frame(count),by.x ="gene_id",by.y="row.names") %>% select(-gene_id) %>% column_to_rownames(var="gene_name")
anno <- data.frame(ID= metadata$barcode, tissue_type = metadata$tissue_type,sample_type = metadata$sample_type) %>% filter(sample_type != "Metastatic") %>% 
  mutate(sample_type = case_when(sample_type == "Primary Tumor"~"Tumor", sample_type=="Solid Tissue Normal" ~ "Adjacent_Normal"))
anno$sample_type = factor(anno$sample_type, levels = c("Tumor","Adjacent_Normal")) 
count = count[-which(rowSums(count)<50),] # 过滤低表达基因
# QC
tcga.data.long = merge(anno,melt(count), by.x="ID", by.y = "Var2") %>% distinct()
ggplot(distinct(tcga.data.long),aes(x=sample_type,y=log2(value+1), fill=sample_type))+ geom_violin(position = position_dodge(width = 1.5)) +geom_boxplot(position = position_dodge(width=0.5)) + theme_bw()

# DEseq2差异分析
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count[,anno$ID]), colData = anno, design=~ sample_type)
res <- DESeq(dds)
# saveRDS(res,file = "TCGA/TCGA_BRCA_DESeq2.res.RData")
res = read_rds("TCGA/TCGA_BRCA_DESeq2.res.RData")
tcga.diff <- results(res,contrast =c("sample_type","Tumor","Adjacent_Normal"))
tcga.diff <- as.data.frame(tcga.diff)
tcga.diff <- merge(as.data.frame(rowdata[,c("gene_id","gene_name")]),tcga.diff,by.x ="gene_id",by.y="row.names") %>% mutate(Regulation =case_when(padj <0.05 & log2FoldChange > log2(1.5)~ "up-regulated", padj < 0.05 & log2FoldChange < -log2(1.5) ~ "down-regulated", TRUE ~ "unchanged"))
table(tcga.diff$Regulation)
colnames(tcga.diff)<- c("gene_id","Gene name", "baseMean","log2FC", "lfcSE","stat","pvalue","adj.P.Val","Regulation")
volcano.plot(tcga.diff,fdr = 0.05,fc=1.5, topn=NA, title="TCGA BRCA Tumor-vs-Normal_ajdacent")

#### CPTAC #### 
pro.data = fread("CPTAC/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv")
pro.data = pro.data %>% select(1,grep("Unshared Log Ratio",colnames(pro.data))) %>% filter(!Gene %in% c("Median","Mean","StdDev"))

cptac.samples = read_excel("CPTAC/S039_Breast_Cancer_Prospective_Collection_Specimens_r1.xlsx") %>% na.omit(Label)
colnames(cptac.samples) = c("Group","Order", "Label","Participants ID","TMT experiments","TMT channels")
s=intersect(cptac.samples$Label, sub(" Unshared Log Ratio","",colnames(pro.data)))
pro.data = pro.data %>% select(1, paste0(s," Unshared Log Ratio"))
colnames(pro.data)[2:ncol(pro.data)] = s
cptac.samples=cptac.samples %>% filter(Label %in% s)
# QC
data.long=melt(pro.data,"Gene",na.rm = T) %>% mutate(variable=sub(variable,pattern=" Unshared Log Ratio",replacement="")) %>% as.data.frame()  %>% merge(cptac.samples[,c("Label","Group")],by.x="variable",by.y="Label")
ggplot(data.long,aes(Group,value))+geom_boxplot()+ylim(-5,5)

# Limma 
MA=subset(pro.data, rowSums(is.na(pro.data))<nrow(cptac.samples)/2,select = c("Gene",as.character(cptac.samples$Label))) %>% column_to_rownames(var = "Gene")
cptac.samples$Group=factor(cptac.samples$Group, levels=c("Tumor","Adjacent_Normal"))
design <- model.matrix(~0+cptac.samples$Group)
compare=c("Tumor","Adjacent_Normal")
colnames(design) <- compare
fit <- lmFit(MA, design)
contrast.matrix <- makeContrasts(contrasts=paste(compare[1],compare[2],sep="-"), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
cptac.diff=topTable(fit2,n=Inf,adjust.method = "BH") %>% dplyr::rename(log2FC="logFC") %>% rownames_to_column(var ="Gene name") %>% 
  mutate(Regulation=case_when(log2FC > log2(1.5) & adj.P.Val < 0.05 ~ "up-regulated", log2FC < -log2(1.5) & adj.P.Val < 0.05 ~ "down-regulated", TRUE ~"unchanged"))
cptac.diff$log2FC=round(cptac.diff$log2FC,3)
cptac.diff$P.Value=signif(cptac.diff$P.Value,3)
cptac.diff$adj.P.Val=signif(cptac.diff$adj.P.Val,3)

volcano.plot(cptac.diff,fdr = 0.05,fc=1.5, topn=NA, title="CPTAC BRCA Tumor-vs-Normal_ajdacent")
list1 = list(TCGA.UP =tcga.diff$`Gene name`[tcga.diff$Regulation=="up-regulated"], 
             CPTAC.UP =cptac.diff$`Gene name`[cptac.diff$Regulation=="up-regulated"])
ggvenn(list1,fill_color = brewer.pal(3,"Set1"))

intersect1 = intersect(tcga.diff$`Gene name`[tcga.diff$Regulation=="up-regulated"], cptac.diff$`Gene name`[cptac.diff$Regulation=="up-regulated"])

# t = intersect(intersect_genes, data.diff$Gene[data.diff$Regulation!="unchanged"])
# mat = pro.data %>% filter(Gene %in% c("CXCL10","VCAN")) %>% column_to_rownames(var = "Gene")
# a = cptac.samples %>% select(Label, Group) %>% arrange(Group) %>% column_to_rownames(var = "Label")
# pheatmap(mat[,rownames(a)], show_colnames = F, cluster_cols = F, annotation_col = a)
# 
# merge(melt(pro.data, "Gene"),cptac.samples[,c("Label","Group")], by.x= "variable", by.y ="Label") %>% filter(Gene %in% c("CXCL10","VCAN")) %>% ggplot(aes(x = Group, y= value, fill=Group))+geom_boxplot() + facet_grid(~Gene) + theme_bw()
# merge(melt(mat %>% rownames_to_column(var="Gene")), a, by.x = "variable",by.y = "row.names") %>% ggplot(aes(x = Group, y = value)) + geom_boxplot() + facet_wrap(~Gene)
# 
# volcano.plot.t(data.diff, fdr=0.05, fc = 2,topn = intersect_genes,title = "Tumor-vs-Normal_adjacent")
# 
# 
# t = data.diff %>% filter(`Gene name` %in% intersect_genes) 
# length(8:ncol(t))-1

#### 导入血浆蛋白数据 #### 
data = fread("FWMS20240721_withoutQC/1-MS_identified_summary/Table 1a. MS identified protein detail.xls")
sample=fread("samples_withoutQCs.txt") %>% mutate(Group = factor(Group, levels = c("Con","BT","MT")))
data.sub=subset(data,select=c("Protein accession","Gene name","Protein description",sample$ID))

DEP.MT.f = fread("FWMS20240721_withoutQC/2-Differentially_expressed_protein/Table 2. DEPs list of MT-vs-Con.xls")
DEP.BT.f = fread("FWMS20240721_withoutQC/2-Differentially_expressed_protein/Table 2. DEPs list of BT-vs-Con.xls")
DEP.MB.f = fread("FWMS20240721_withoutQC/2-Differentially_expressed_protein/Table 2. DEPs list of MT-vs-BT.xls") %>% mutate(Regulation=case_when(log2FC>0 &P.Value < 0.05~ "up-regulated", log2FC<0&P.Value<0.05~"down-regulated", TRUE~"unchanged"))

volcano.plot1(DEP.MT.f, pval=0.05, fc = 1.3,topn = NA,title = "MT-vs-Con")
volcano.plot1(DEP.BT.f, pval=0.05, fc = 1.3,topn = NA,title = "BT-vs-Con")
volcano.plot1(DEP.MB.f, pval=0.05, fc = 1.1,title = "MT-vs-BT")

ggarrange(
  volcano.plot1(DEP.BT.f, pval=0.05, fc = 1.3,topn = NA,title = ""),
  volcano.plot1(DEP.MT.f, pval=0.05, fc = 1.3,topn = NA,title = ""),
  volcano.plot1(DEP.MB.f, pval=0.05, fc = 1.1,title = ""),
  nrow = 3,labels = c("BT vs Con","MT vs Con","MB vs BT"),
  label.x = 0.2
)

#---------------------- GSEA ----------------------#
kegg.data=fread("C:/Users/css/Desktop/Rstudio/Rshiny/DIA/DB/Homo_sapiens_9606_proteome_sp_20240620/KEGG_annotation.tsv")
kegg.gtm=lapply(split(kegg.data,kegg.data$Terms),function(x) x$`Protein accession`)

enrich.test=function(x){
  a=x[1]
  b=x[2]-a
  c=x[3]-a
  d=x[4]-a-b-c
  p.value=fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")$p.value
}
enrich.demo=function(go.data,diff.list){
  
  go.data.diff=merge(go.data,subset(diff.list,select=c("Protein accession","Gene name")),by="Protein accession")
  temp.s=go.data.diff %>% group_by(Terms) %>% mutate(`Gene names`=paste(`Gene name`,collapse = ";"),`Proteins`=paste(`Protein accession`,collapse = ';')) %>%
    ungroup() %>% distinct(Terms,`Gene names`,Proteins)
  
  mapping=data.frame(table(go.data.diff$Terms))
  background=data.frame(table(go.data$Terms))
  go.data.state=merge(mapping,background,by="Var1")
  colnames(go.data.state)=c("Terms","mapping","backgroupd")
  go.data.state$`All mapping`=length(unique(go.data.diff$`Protein accession`))
  go.data.state$`All background`=length(unique(go.data$`Protein accession`))
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$backgroupd/go.data.state$`All mapping`*go.data.state$`All background`,2)
  go.data.state=subset(go.data.state,mapping>1&`Fold enrich`>1)
  go.data.state$p.value=signif(apply(go.data.state[2:5],1,enrich.test),3)
  go.data.state$fdr=signif(p.adjust(go.data.state$p.value),3)
  go.data.state=subset(go.data.state,select=c("Terms","mapping","Fold enrich","p.value","fdr"))
  go.data.state=merge(go.data.state,temp.s,by="Terms")
  return(go.data.state)
}

enrich.out=rbindlist(
  lapply(list(BT=DEP.BT.f %>% filter(Regulation!='unchanged'),
              MT=DEP.MT.f %>% filter(Regulation!='unchanged'),
              MB=DEP.MB.f %>% filter(Regulation!='unchanged')
  ),function(x){
    rbindlist(lapply(split(x,x$Regulation),function(y){
      enrich.demo(kegg.data %>% filter(`Protein accession` %in% data.sub$`Protein accession`),y)
    }),idcol="Regulation")
  }),idcol="Type"
) %>% mutate(Type=factor(Type,levels=c("BT","MT","MB")))

ggplot(enrich.out %>% filter(p.value<0.01,stringr::str_length(Terms)<40),aes(Type,Terms))+geom_point(aes(size=mapping,color=p.value))+facet_grid(.~Regulation,scales = 'free')+
  scale_colour_material("red",reverse = T)+theme_bw(base_size = 14)

#---------------------- Intersect  ----------------------#
list2 = list(
             BT_UP=DEP.BT.f$`Protein accession`[DEP.BT.f$Regulation=='up-regulated'],
             MT_UP=DEP.MT.f$`Protein accession`[DEP.MT.f$Regulation=='up-regulated'],
             MB_UP=DEP.MB.f$`Protein accession`[DEP.MB.f$Regulation=='up-regulated'])

ggvenn(list2, fill_color = brewer.pal(3, "Set1"))

intersect2=intersect(
  intersect(DEP.BT.f$`Protein accession`[DEP.BT.f$Regulation=='up-regulated'],DEP.MT.f$`Protein accession`[DEP.MT.f$Regulation=='up-regulated']),
  DEP.MB.f$`Protein accession`[DEP.MB.f$Regulation=='up-regulated']
)
intersect_genes = data$`Gene name`[data$`Protein accession` %in% intersect2]

list3 = list(TCGA_CPTAC.intersect = intersect1, BT_MT_MB.intersect = intersect_genes)
ggvenn(list3, fill_color = brewer.pal(3, "Set2"))
intersect3 = intersect(intersect_genes, intersect1)
sample$Group = factor(sample$Group, levels = c("Con","BT","MT"))
ggarrange(
  ggplot(merge(sample,melt(data,id.vars=1:8) %>% filter(`Gene name` %in% intersect3),by.x="ID",by.y="variable"))+
    geom_tile(aes(ID,`Gene name`,fill=value))+scale_fill_gradientn(colors = c("steelblue","white","red"),limits=c(-4,4),oob=scales::squish)+
    facet_grid(.~Group,space = 'free',scales = 'free')+theme_test()+theme(axis.text.x=element_blank()),
  ggplot(merge(sample,melt(data,id.vars=1:8) %>% filter(`Gene name` %in% intersect3),by.x="ID",by.y="variable"))+
    geom_boxplot(aes(Group,value,fill=Group))+facet_wrap(.~`Gene name`,scales = 'free')+theme_test(),
  nrow = 1
)

# MT.up=DEP.MT.f$`Gene name`[DEP.MT.f$P.Value<0.05 & DEP.MT.f$log2FC>0]
# BT.up=DEP.BT.f$`Gene name`[DEP.BT.f$P.Value<0.05 & DEP.BT.f$log2FC>0]
# MB.up=DEP.MB.f$`Gene name`[DEP.MB.f$P.Value<0.05 & DEP.MB.f$log2FC>0]

## 生存分析 ## 
data.rna=fread("HiSeqV2")
data.rna.os=fread("data_bcr_clinical_data_patient_brca.txt") %>%
  mutate(`Overall Survival Status`=case_when(`Overall Survival Status`=='0:LIVING'~0,`Overall Survival Status`=='1:DECEASED'~1)) %>%
  filter(!is.na(`Overall Survival Status`)) %>% select(`Overall Survival Status`,`Overall Survival (Months)`,`Patient Identifier`) %>%
  mutate(`Overall Survival (Months)`=as.numeric(`Overall Survival (Months)`))
data.rna.long=melt(subset(data.rna,sample %in% list),id.vars = "sample") %>% filter(grepl(variable,pattern="-01$")) %>% mutate(`Patient Identifier`=sub(variable,pattern="-01$",replacement="")) %>%
  merge(data.rna.os,by="Patient Identifier")
data.rna.long.cox=rbindlist(lapply(split(data.rna.long, data.rna.long$sample),function(t){
  as.data.frame(summary(coxph(Surv(time = `Overall Survival (Months)`,event = `Overall Survival Status`)~value,data = t))$coef)
}),idcol = "Gene name")
ggplot(merge(sample,melt(data,id.vars=1:8) %>% filter(`Gene name` %in% list2),by.x="ID",by.y="variable"))+
  geom_tile(aes(ID,`Gene name`,fill=value))+scale_fill_gradientn(colors = c("steelblue","white","red"),limits=c(-4,4),oob=scales::squish)+
  facet_grid(.~Group,space = 'free',scales = 'free')+theme_test()+theme(axis.text.x=element_blank())
ggplot(data.rna.long.cox %>% filter(`Pr(>|z|)`<0.05))+geom_point(aes(coef,`Gene name`))
ggsurvplot(surv_fit(Surv(time = `Overall Survival (Months)`,event = `Overall Survival Status`)~CXCL10,data=data.rna.long %>% filter(sample=='CXCL10') %>% mutate(CXCL10=case_when(value>median(value)~'High',TRUE~'Low'))), pval=T, pval.method = T)
ggsurvplot(surv_fit(Surv(time = `Overall Survival (Months)`,event = `Overall Survival Status`)~VCAN,data=data.rna.long %>% filter(sample=='VCAN') %>% mutate(VCAN=case_when(value>median(value)~'High',TRUE~'Low'))), pval=T, pval.method = T)
 

## 临床信息 
meta = read_excel("乳腺癌蛋白组样本明细-云启系统-20240813.xlsx", sheet = 1)
colnames(meta) = c("Name","ID","Group","Number","Sample Type", "Info")

clinical = read_excel("乳腺癌项目临床信息20240226.xlsx", sheet=1) %>% filter(`姓名` %in% meta$Name)
health = read_excel("乳腺癌项目收样信息-健康人群20240809update(1).xlsx")

## 筛选到的标志物做机器学习 --------------------------------------------
markers = fread("FWMS20240721_withoutQC/1-MS_identified_summary/Table 1a. MS identified protein detail.xls") %>% filter(`Gene name` %in% intersect3) %>% 
  column_to_rownames(var = "Gene name") %>% select(sample$Sample) %>%  t() %>% as.matrix()
markers[is.na(markers)] = 0
markers = merge(sample[,c("Sample","Group")], as.data.frame(markers) %>% rownames_to_column(var="Sample"),by = "Sample") %>% column_to_rownames(var = "Sample")

write.table(markers, file = "Machine Learning/Markers.txt",sep = "\t", row.names = F)

markers.info = data.sub %>% filter(`Gene name` %in% intersect3) %>% select(`Protein accession`, `Gene name`,`Protein description`)
bt = DEP.BT.f %>% filter(`Gene name` %in% intersect3) %>% select(`Protein accession`, log2FC, P.Value) 
colnames(bt)[2:3]= paste0("BT-vs-Con.", colnames(bt)[2:3])
mt = DEP.MT.f %>% filter(`Gene name` %in% intersect3) %>% select(`Protein accession`, log2FC, P.Value) 
colnames(mt)[2:3]= paste0("MT-vs-Con.", colnames(mt)[2:3])
mb = DEP.MB.f %>% filter(`Gene name` %in% intersect3) %>% select(`Protein accession`, log2FC, P.Value) 
colnames(mb)[2:3]= paste0("MB-vs-Con.", colnames(mb)[2:3])
markers.info = merge(merge(merge(markers.info, bt, by="Protein accession"),mt,by="Protein accession"),mb,by="Protein accession")
write.table(markers.info, file ="markers.xls", sep ="\t",row.names = F)

# test1 = pro.data %>% filter(Gene %in% list2) %>% column_to_rownames(var="Gene") %>% t() %>% as.matrix()
# test1 = merge(cptac.samples[,c("Label","Group")], as.data.frame(test1) %>% rownames_to_column(var="Label"), by="Label") %>% column_to_rownames(var ="Label") 
# write.table(test1, file = "Machine Learning/test_data_from_CPTAC.txt", sep="\t", row.names = F)

count.tmp = as.data.frame(count.matrix) %>% filter(gene_name %in% n) 
t = count.tmp[,-1] %>% column_to_rownames(var = "gene_name") 
a = anno %>% column_to_rownames(var = "ID") %>% select(sample_type) %>% arrange(sample_type)
pheatmap(log2(t[,rownames(a)]+1),scale = "row",show_colnames = F, annotation_col = a, cluster_col=F)

merge(a, melt(t %>% rownames_to_column(var = "Gene name")), by.x ="row.names", by.y = "variable") %>% ggplot(aes(x = `Gene name`, y = log2(value+1), fill = sample_type)) + geom_boxplot()


## 11.26 再看一下早晚期有没有差异，然后看看模型鉴别 良性 vs 早期，健康+良性 vs 早期的性能 -------------------
library(ggpubr)
library(ggplot2)

clinical_data = read_excel("萧一乳腺癌项目临床信息.xlsx")
meta = read_excel("乳腺癌蛋白组样本明细-云启系统-20240813.xlsx", sheet = 1)
colnames(meta) = c("Name","ID","Group","Number","Sample Type", "Info")

meta1 = merge(meta, subset(clinical_data, select=c("姓名","TNM分期")), by.x = "Name", by.y ="姓名", all.x = T) %>% 
  mutate(Tumor_stage = case_when(`TNM分期`=="/"~"Benign", `TNM分期`%in% c("I","II")~ "Early",`TNM分期`%in% c("III","IV")~"Late", TRUE ~"HC"))
pro.data = fread("FWMS20240721_withoutQC/1-MS_identified_summary/Table 1a. MS identified protein detail.xls")
samples = fread("samples_withoutQCs.txt") %>% merge(subset(meta1, select=c("ID","Tumor_stage")), by = "ID") %>% dplyr::select(-ID)

pro.exp = pro.data %>% column_to_rownames(var = "Protein accession") %>% dplyr::select(samples$Sample) %>% na.omit()
annot = samples %>%  column_to_rownames(var="Sample")
annot_color = list(Group=pal_npg()(length(unique(annot$Group))), Tumor_stage=pal_d3("category20c")(length(unique(annot$Tumor_stage)))) 
names(annot_color$Tumor_stage)=c("Benign","Early","HC","Late")
names(annot_color$Group)=unique(annot$Group)
p1 = pheatmap(cor(pro.exp,use="pairwise.complete.obs",method = 'spearman'),annotation_col = annot,border_color = 'white',
              show_rownames =F,annotation_colors = annot_color,silent = T,show_colnames = F)
data.pca=prcomp(t(na.omit(pro.exp)), scale=T, center = T)

p2 = ggplot(merge(samples,as.data.frame(data.pca$x) %>% rownames_to_column(var="Sample"),by="Sample"),aes(PC1,PC2))+
  geom_point(aes(color=Tumor_stage),size=3)+scale_color_d3("category20c")+theme_classic() +geom_text_repel(aes(label = Sample))
pro.qc.plot = ggarrange(p1$gtable,p2,nrow = 1,labels = "AUTO")
pro.qc.plot
ggsave(filename = "Late-vs-Early/PCC and PCA plot.pdf", plot =pro.qc.plot, width = 12, height =6)
ggsave(filename = "Late-vs-Early/PCC and PCA plot.png", plot =pro.qc.plot, width = 12, height =6, bg="white")

limma.filter=function(data.quant,s,adj="BH"){
  data.diff=as.data.frame(subset(data.quant,select = c("Protein accession","Gene name","Protein description",s$ID)))
  MA=subset(data.diff,rowSums(is.na(data.diff))<nrow(s)/2,select=as.character(s$ID))
  design <- model.matrix(~ 0+s$Type)
  compare=levels(s$Type)
  colnames(design) <- compare
  fit <- lmFit(MA, design)
  contrast.matrix <- makeContrasts(contrasts=paste(compare[1],compare[2],sep="-"), levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  out=topTable(fit2,n=Inf,adjust.method = adj) %>% dplyr::rename(log2FC="logFC")
  data.diff=merge(data.diff,subset(out,select=c("log2FC","P.Value","adj.P.Val")),by="row.names")
  data.diff$log2FC=round(data.diff$log2FC,3)
  data.diff$P.Value=signif(data.diff$P.Value,3)
  data.diff$adj.P.Val=signif(data.diff$adj.P.Val,3)
  data.diff=data.diff[-1]
  return(list(data=data.table(data.diff),sample=s))
}

## 晚期vs早期
s = samples %>% filter(Tumor_stage %in% c("Early","Late")) %>% select(-Group) %>%  mutate(Type=factor(Tumor_stage, levels =c("Late","Early"))) %>% select(-Tumor_stage)
colnames(s) = c("ID","Type")
s = s[order(Type),]
diff.res = limma.filter(pro.data, s, "BH") 
diff.data = diff.res$data %>% mutate(Regulation = case_when(log2FC>log2(1.3)&P.Value<0.05~ "Up", log2FC< -log2(1.3)&P.Value<0.05~"Down",TRUE~"Unchanged")) 
table(diff.data$Regulation)
markers = fread("Machine Learning/Markers.txt")
intersect = intersect(diff.data$`Gene name`[diff.data$Regulation!="Unchanged"], colnames(markers))
write.table(diff.data, file = "Late-vs-Early/Differentially expressed proteins in Late-vs-Early.xls", sep = "\t", row.names = F)
p = diff.data %>% mutate(label = case_when(`Gene name` %in% intersect ~ `Gene name`, TRUE ~"")) %>% ggplot(aes(x=log2FC, y= -log10(P.Value), color = Regulation)) + geom_point() + 
  scale_color_manual(values = c("Up"="#E64B35FF", "Down"= "#4DBBD5FF", "Unchanged"="grey")) + theme_classic() + 
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)), linetype=2) +
  geom_point(data=diff.data %>% filter(`Gene name` %in% intersect),aes(log2FC, -log10(P.Value)), color = "gold", size = 3.5) +
  geom_text(aes(label = label), color="black", hjust = 1, vjust =1 , size = 3.5)
ggsave(filename = "Late-vs-Early/Volcano plot of Late-vs-Early.pdf", width = 6, height = 6, plot = p)
ggsave(filename = "Late-vs-Early/Volcano plot of Late-vs-Early.png", width = 6, height = 6, plot = p)

## 良性vs早期
s = samples %>% filter(Tumor_stage %in% c("Benign","Early")) %>% select(-Group) %>%  mutate(Type=factor(Tumor_stage, levels =c("Benign","Early"))) %>% select(-Tumor_stage)
colnames(s) = c("ID","Type")
s = s[order(Type),]
diff.res = limma.filter(pro.data, s, "BH") 
diff.data = diff.res$data %>% mutate(Regulation = case_when(log2FC>log2(1.3)&P.Value<0.05~ "Up", log2FC< -log2(1.3)&P.Value<0.05~"Down",TRUE~"Unchanged")) 
diff.data %>% ggplot(aes(x=log2FC, y= -log10(P.Value), color = Regulation)) + geom_point() + scale_color_manual(values = c("Up"="#E64B35FF", "Down"= "#4DBBD5FF", "Unchanged"="grey")) + theme_classic()
table(diff.data$Regulation)
intersect= intersect(diff.data$`Gene name`[diff.data$Regulation!="Unchanged"], colnames(markers))

write.table(diff.data, file = "Late-vs-Early/Differentially expressed proteins in Benign-vs-Early.xls", sep = "\t", row.names = F)
p = diff.data %>% mutate(label = case_when(`Gene name` %in% intersect ~ `Gene name`, TRUE ~"")) %>% ggplot(aes(x=log2FC, y= -log10(P.Value), color = Regulation)) + geom_point() + 
  scale_color_manual(values = c("Up"="#E64B35FF", "Down"= "#4DBBD5FF", "Unchanged"="grey")) + theme_classic() + 
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)), linetype=2) +
  geom_point(data=diff.data %>% filter(`Gene name` %in% intersect),aes(log2FC, -log10(P.Value)), color = "gold", size = 3.5) +
  geom_text(aes(label = label), color="black", hjust = 1, vjust =1 , size = 3.5)
ggsave(filename = "Late-vs-Early/Volcano plot of Benign-vs-Early.pdf", width = 6, height = 6, plot = p)
ggsave(filename = "Late-vs-Early/Volcano plot of Benign-vs-Early.png", width = 6, height = 6, plot = p)

## 健康+良性vs早期
s = samples %>% filter(Tumor_stage %in% c("HC","Benign","Early")) %>% select(-Group) %>% mutate(Tumor_stage=case_when(Tumor_stage %in% c("HC","Benign")~"BH",Tumor_stage=="Early"~"Early")) %>%  
  mutate(Type=factor(Tumor_stage, levels =c("BH","Early"))) %>% select(-Tumor_stage)
colnames(s) = c("ID","Type")
s = s[order(Type),]
diff.res = limma.filter(pro.data, s, "BH") 
diff.data = diff.res$data %>% mutate(Regulation = case_when(log2FC>log2(1.3)&P.Value<0.05~ "Up", log2FC< -log2(1.3)&P.Value<0.05~"Down",TRUE~"Unchanged")) 
diff.data %>% ggplot(aes(x=log2FC, y= -log10(P.Value), color = Regulation)) + geom_point() + scale_color_manual(values = c("Up"="#E64B35FF", "Down"= "#4DBBD5FF", "Unchanged"="grey")) + theme_classic()
table(diff.data$Regulation)
intersect=intersect(diff.data$`Gene name`[diff.data$Regulation!="Unchanged"], colnames(markers))
write.table(diff.data, file = "Late-vs-Early/Differentially expressed proteins in BH-vs-Early.xls", sep = "\t", row.names = F)
p = diff.data %>% mutate(label = case_when(`Gene name` %in% intersect ~ `Gene name`, TRUE ~"")) %>% ggplot(aes(x=log2FC, y= -log10(P.Value), color = Regulation)) + geom_point() + 
  scale_color_manual(values = c("Up"="#E64B35FF", "Down"= "#4DBBD5FF", "Unchanged"="grey")) + theme_classic() + 
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)), linetype=2) +
  geom_point(data=diff.data %>% filter(`Gene name` %in% intersect),aes(log2FC, -log10(P.Value)), color = "gold", size = 3.5) +
  geom_text(aes(label = label), color="black", hjust = 1, vjust =1 , size = 3.5)
ggsave(filename = "Late-vs-Early/Volcano plot of BH-vs-Early.pdf", width = 6, height = 6, plot = p)
ggsave(filename = "Late-vs-Early/Volcano plot of BH-vs-Early.png", width = 6, height = 6, plot = p)

## 根据早晚期做机器学习的建模
selected_markers = colnames(markers)[-1]
markers = fread("FWMS20240721_withoutQC/1-MS_identified_summary/Table 1a. MS identified protein detail.xls") %>% filter(`Gene name` %in% selected_markers) %>% 
  column_to_rownames(var = "Gene name") %>% select(samples$Sample) %>%  t() %>% as.matrix()
markers[is.na(markers)] = 0
markers = merge(samples[,c("Sample","Tumor_stage")], as.data.frame(markers) %>% rownames_to_column(var="Sample"),by = "Sample") %>% column_to_rownames(var = "Sample")
colnames(markers)[1] = "Group"
write.table(markers, file = "Machine Learning/Markers2.txt",sep = "\t", row.names = F)

### 6.24 10个marker做富集分析，以及热图展示 -----------------------------------
library(ggh4x)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(readxl)

dir.create("enrich_res")
markers = fread("markers.xls")
entrezid = bitr(markers$`Gene name`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

kegg.res = enrichKEGG(entrezid$ENTREZID,organism = "hsa",pvalueCutoff = 1)
kegg.res = setReadable(kegg.res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg.res@result$GeneRatio = apply(kegg.res@result, 1,function(x){
  round(eval(parse(text=x[["GeneRatio"]])),2)
})
write.table(kegg.res@result,file = paste0("enrich_res/KEGG enrich results.xls"), sep="\t", row.names = F)

go.res = enrichGO(entrezid$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 1, readable = T)
go.res@result$GeneRatio = apply(go.res@result, 1,function(x){
  round(eval(parse(text=x[["GeneRatio"]])),2)
})
write.table(go.res@result, file = paste0("enrich_res/GO BP enrich results.xls"), sep = "\t", row.names = F)

rea.res = enrichPathway(gene = entrezid$ENTREZID, pvalueCutoff = 1, readable=TRUE)
rea.res@result$GeneRatio = apply(rea.res@result, 1,function(x){
  round(eval(parse(text=x[["GeneRatio"]])),2)
})
write.table(rea.res@result,file = paste0("enrich_res/Reactome enrich results.xls"), sep ="\t", row.names = F)

go.ids = c("GO:0001649","GO:0045669","GO:0033627","GO:0060562","GO:0001503","GO:0045667","GO:0045785","GO:0009612","GO:0006909","GO:0042692")
kegg.ids = c("hsa04514","hsa04623","hsa04512","hsa04061","hsa04670","hsa04668","hsa04650","hsa04390","hsa04141","hsa04310")
rea.ids = c("R-HSA-3772470","R-HSA-1474244","R-HSA-449147","R-HSA-216083","R-HSA-6785807","R-HSA-198933","R-HSA-202733","R-HSA-201681","R-HSA-195721","R-HSA-6798695")

go = go.res@result %>% filter(ID %in% go.ids) %>% mutate(Category = "Biological process")
kegg = kegg.res@result %>% filter(ID %in% kegg.ids) %>% dplyr::select(-c(category, subcategory)) %>% mutate(Category = "KEGG")
rea = rea.res@result %>% filter(ID %in% rea.ids) %>% mutate(Category = "Reactome")
cols = enrich.res$Color
pathways=c("osteoblast differentiation","epithelial tube morphogenesis","ossification","Cell adhesion molecules","Cytosolic DNA-sensing pathway",
           "ECM-receptor interaction","Hippo signaling pathway","Wnt signaling pathway","Extracellular matrix organization","Signaling by Interleukins")
enrich.res = rbind(go, rbind(kegg,rea)) %>% arrange(Category,GeneRatio) %>% mutate(Color = case_when(Description %in% pathways~"#E41A1C" ,TRUE~"#377EB8")) 
enrich.res$Description=factor(enrich.res$Description, levels = enrich.res$Description)
# enrich.res$Color = factor(enrich.res$Color,levels = unique(enrich.res$Color))
strip =strip_themed(background_x = elem_list_rect(fill = brewer.pal(3,"Set2")))
n.color=enrich.res$Color
names(n.color)=enrich.res$Description
p1 = ggplot(enrich.res,aes(x =GeneRatio, y =Description, color = pvalue, size = Count)) + geom_point() +
  facet_wrap2(~Category, nrow = 3,scales = "free", strip = strip) + 
  scale_color_material(palette = "red",reverse = T) +
  theme_minimal(base_size = 12) + ylab("Pathways") + 
  scale_size(range = c(3,8)) + 
  scale_y_discrete(
    labels = function(x) str_wrap(x, width = 45)  # 调整 width 控制换行位置
  ) + 
  theme(axis.text.y=element_text(color = enrich.res$Color))

enrich = enrich.res %>% filter(Description %in% pathways) %>% separate_rows(geneID, sep = "\\/") %>% dplyr::select(Category,Description,geneID) 
colors = enrich.res %>% filter(Description %in% pathways) %>% 
  mutate(Color = case_when(Category=="Biological process"~"#66C2A5", Category=="KEGG"~"#FC8D62",TRUE~"#8DA0CB"))
sample=fread("samples_withoutQCs.txt") %>% mutate(Group = case_when(Group=="Con"~"HC",Group=="BT"~"BBD",Group=="MT"~"BC")) %>% mutate(Group = factor(Group, levels = c("HC","BBD","BC")))
pro.data = fread("FWMS20240721_withoutQC/1-MS_identified_summary/Table 1a. MS identified protein detail.xls") %>% filter(`Gene name` %in% markers$`Gene name`) %>% 
  column_to_rownames(var = "Gene name") %>% dplyr::select(sample$Sample)

data.long= fread("FWMS20240721_withoutQC/1-MS_identified_summary/Table 1a. MS identified protein detail.xls") %>% filter(`Gene name` %in% markers$`Gene name`) %>% 
  dplyr::select(`Gene name`, sample$Sample) %>% melt(measure.vars = sample$Sample, variable.name = "Sample") %>% merge(sample, by="Sample") %>% 
  group_by(`Gene name`) %>% mutate(value=value-mean(value,na.rm=T)) 
data.long.new = merge(data.long, enrich,by.x = "Gene name", by.y = "geneID") 
strip1 = strip_themed(background_x = elem_list_rect(fill = pal_npg()(3)),
                      background_y = elem_list_rect(fill = colors$Color))
p2 = ggplot(data.long.new,aes(x = Sample,y = `Gene name`,fill = value))+
  geom_raster()+
  facet_grid2(Description~Group,scales = "free",strip = strip1,space = 'free')+
  scale_fill_gradientn(colors = c(pal_material("blue",reverse = T)(10),'white',
                                  pal_material("red")(10)),limits=c(-2,2),oob=scales::squish,na.value = "white")+
  theme_bw(base_size = 12)+
  theme(panel.border = element_blank(),
        strip.text.y.right =element_text(angle = 0),
        axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(size = 9,angle = 75,hjust = 1,vjust = 1),
        legend.position = "right",
        axis.ticks.y = element_blank())
p = ggarrange(p1,p2, labels = c("A","B"))
ggsave(filename = "enrich_res/Enrich plot and Heatmap.pdf", width = 16,height =9, plot =p)






















