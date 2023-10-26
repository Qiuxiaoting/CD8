#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)               #???冒?
expFile="TCGA.share.txt"     #?????????募?
cliFile="time.txt"           #?????????募?
setwd("D:\\Nature\\CD8\\CD8\\17.mergeTime")      #???霉???目录

# 从 expFile 中读取基因表达数据并进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F) # 从文件中读取数据，包括表头，使用制表符分隔。
rt=as.matrix(rt)
rownames(rt)=rt[,1]# 设置矩阵的行名为第一列的数据。
exp=rt[,2:ncol(rt)]# 提取基因表达数据，去除第一列的标识。
dimnames=list(rownames(exp),colnames(exp))# 设置数据的行名和列名。
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)# 将数据转换为数值矩阵
data=avereps(data) # 对数据进行平均处理。

#鍘婚櫎姝ｅ父鏍锋湰 閮芥槸tnbc鏍锋湰 璇曡瘯涓嶅幓闄?
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2","1",group)
#data=data[,group==0]
#colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
#data=t(data)
#data=avereps(data) 

#??取?????????募?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

# 找到共有的样本数据并整合时间数据和基因表达数据。
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli, data)# 将时间数据和基因表达数据按列合并。

# 添加样本标识列，并将整合后的数据写入文件 "TCGA.expTime.txt"。
out=cbind(id=row.names(out), out)# 添加样本标识列。
write.table(out, file="TCGA.expTime.txt", sep="\t", row.names=F, quote=F)


######??????学??: https://www.biowolf.cn/
######?纬?链??1: https://shop119322454.taobao.com
######?纬?链??2: https://ke.biowolf.cn
######?纬?链??3: https://ke.biowolf.cn/mobile
######?饪???师????: seqbio@foxmail.com
######?饪???师微??: eduBio

