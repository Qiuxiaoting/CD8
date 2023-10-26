#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#???ð?
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

riskFile="risk.TCGA.txt"     #?????ļ?
expFile="merge.txt"         #?????????ļ?
gene="CD274"                 #?????ı?׼????
setwd("/Users/apple/Desktop/brca t cell/28.riskGene")      #???ù???Ŀ¼

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
	
#ɾ????????Ʒ
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2", "1", group)
#data=data[,group==0]
	
#??ȡĿ??????????��
data=rbind(data, gene=data[gene,])
exp=t(data[c("gene",gene),])
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
	
#??ȡ?????????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
#?ϲ?????
sameexp=t(exp)
Sample=intersect(row.names(exp), row.names(risk))
exp=exp[sameSample,]
exp=log2(exp+1)
risk=risk[sameSample,]
data=cbind(as.data.frame(exp), as.data.frame(risk))
	
#???ñȽ???
data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
	
#????????ͼ
boxplot=ggboxplot(data, x="Risk", y="gene", color="Risk",
			      xlab="",
			      ylab=paste0(gene, " expression"),
			      legend.title="",
			      palette = c("blue", "red"),
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)	
#????????ͼ
pdf(file="boxplot.pdf", width=4.5, height=4.25)
print(boxplot)
dev.off()
	
	
#?????Է???
x=as.numeric(data[,"riskScore"])
y=as.numeric(data[,gene])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		xlab("Risk score") + ylab(paste0(gene, " expression"))+
		geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
#??????????ͼ??
pdf(file="cor.pdf", width=4.5, height=4.25)
print(p2)
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

