#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)              #???ð?
expFile="GEO.share.txt"     #?????????ļ?
cliFile="time.xlsx"          #?????????ļ?
setwd("D:\\Nature\\CD8\\CD8\\18.geoMergeTime")     #???ù???Ŀ¼

#??ȡ?????ļ????????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)# �����ݽ���ת�ã���ԭ�����б���У�ԭ�����б���С�

#??ȡ?????????ļ?
#cli=readxl:("cliFile", header=T, sep="\t", check.names=F, row.names=1)
cli=readxl::read_xlsx("time.xlsx",col_names = TRUE,)
cli=as.matrix(cli)
rownames(cli)=cli[,1]
cli<-gsub("-",".",cli)# ��ʱ�������еĶ̺����滻Ϊ�㡣
#???ݺϲ?
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)

#?????ϲ????Ľ???
out=cbind(id=row.names(out),out)
write.table(out,file="GEO.expTime.txt",sep="\t",row.names=F,quote=F)


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?�???ʦ????: seqbio@foxmail.com
######?�???ʦ΢??: eduBio

