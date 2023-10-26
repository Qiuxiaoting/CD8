#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
devtools::install_github("YuLab-SMU/DOSE")
devtools::install_github("YuLab-SMU/HDO.db")   
devtools::install_github('YuLab-SMU/clusterProfiler')

#???ð?
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="merge.txt"         #?????????ļ?
riskFile="risk.TCGA.txt"     #?????ļ?
gmtFile="c2.cp.kegg.Hs.symbols.gmt"     #???????ļ?
setwd("/Users/apple/Desktop/brca t cell/25.GSEA")     #???ù???Ŀ¼

#??ȡ?ļ?,?????????ļ?????????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#分组信息
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2", "1", group)
#data=data[,group==0]
#data=t(data)
#rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
#data=t(avereps(data))

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(risk), colnames(data))
risk=risk[sameSample,]
data=data[,sameSample]

#?ߵͷ??ձȽϣ??õ?logFC
dataL=data[,row.names(risk[risk[,"Risk"]=="low",])]
dataH=data[,row.names(risk[risk[,"Risk"]=="high",])]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)

#??ȡ???????ļ?
gmt=read.gmt(gmtFile)

#????GSEA????????
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1, minGSSize=15, maxGSSize=500)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]      #?Ը????Ľ??????й???
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)
	
#?????߷????鸻????ͼ??
termNum=5     #????չʾͨ·????Ŀ??չʾǰ5??????????????ͨ·
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]      #??ȡչʾͨ·??????
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in high risk group")
	pdf(file="GSEA.highrisk.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}

#?????ͷ??ո?????ͼ??
termNum=5     #????չʾͨ·????Ŀ??չʾǰ5??????????????ͨ·
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]      #??ȡչʾͨ·??????
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low risk group")
	pdf(file="GSEA.lowrisk.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

