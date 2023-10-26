#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)              #???冒?
expFile="GEO.share.txt"     #?????????募?
cliFile="time.xlsx"          #?????????募?
setwd("D:\\Nature\\CD8\\CD8\\18.geoMergeTime")     #???霉???目录

#??取?????募????????????募?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)# 将数据进行转置，将原来的行变成列，原来的列变成行。

#??取?????????募?
#cli=readxl:("cliFile", header=T, sep="\t", check.names=F, row.names=1)
cli=readxl::read_xlsx("time.xlsx",col_names = TRUE,)
cli=as.matrix(cli)
rownames(cli)=cli[,1]
cli<-gsub("-",".",cli)# 将时间数据中的短横线替换为点。
#???莺喜?
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)

#?????喜????慕???
out=cbind(id=row.names(out),out)
write.table(out,file="GEO.expTime.txt",sep="\t",row.names=F,quote=F)


######??????学??: https://www.biowolf.cn/
######?纬?链??1: https://shop119322454.taobao.com
######?纬?链??2: https://ke.biowolf.cn
######?纬?链??3: https://ke.biowolf.cn/mobile
######?饪???师????: seqbio@foxmail.com
######?饪???师微??: eduBio

