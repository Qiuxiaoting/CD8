#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#???冒?
library(limma)
library(sva)
# 定义TCGA和GEO基因表达数据文件的文件名
tcgaExpFile="merge.txt"       #TCGA?????????募?
geoExpFile="selectedtnbc_expression.txt"     #GEO?????????募?
geneFile="interGene.txt"       #?????????斜??募?
setwd("D:\\Nature\\CD8\\CD8\\16.intersect")     #???霉???目录

# 读取和处理TCGA基因表达数据
rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1] 杩欎釜鍛戒护浼氬皢琛屽悕娑堝け锛屽彲鑳芥槸璇诲彇鏂囦欢闂锛屼繚璇佽鍚嶆槸鍩哄洜鍚嶅氨琛?
exp=rt[,2:ncol(rt)]#[, 2:ncol(rt)] 表示选择 rt 数据框中的所有行（空的第一个逗号）和从第2列到最后一列的所有列
dimnames=list(rownames(exp),colnames(exp))#rownames(exp) 返回 exp 矩阵的行名，即行索引的标签; colnames(exp) 返回 exp 矩阵的列名，即列索引的标签
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)#avereps(tcga) 函数将检测并去除矩阵 tcga 中的重复行，并计算重复行的平均值。这样，每个唯一的行都只会出现一次，而重复的行会被合并为一个平均值
tcga=log2(tcga+1)#log2(tcga+1) 是对 tcga 矩阵的每个元素执行对数变换的操作。log2 表示以2为底的对数变换，这是常用的对数变换之一。tcga+1 首先在每个元素上加1，以避免出现0值，因为对数不能应用于0或负数。

# 读取和处理GEO基因表达数据
rt=read.table(geoExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]

exp=rt[,3:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)

# 对GEO数据进行处理，包括log2转换和规范化
qx=as.numeric(quantile(geo, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))#c(0, 0.25, 0.5, 0.75, 0.99, 1.0)：表示要计算的分位数  v
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )#qx[5]>100 检查第99百分位数（qx[5]）是否大于100，qx[6]-qx[1]>50 检查第99百分位数和第0百分位数之间的差是否大于50，qx[2]>0 检查第一四分位数（qx[2]）是否大于0
if(LogC){
    geo[geo<0]=0# 将 geo 数据矩阵中小于0的值设为0，即对数据进行截尾，以确保所有数据都不小于0
    geo=log2(geo+1)}
geo=normalizeBetweenArrays(geo)

# 取TCGA和GEO数据共同的基因，进行进一步处理
sameGene=intersect(row.names(tcga),row.names(geo))
tcgaOut=tcga[sameGene,]
geoOut=geo[sameGene,]

# 合并TCGA和GEO数据
all=cbind(tcgaOut,geoOut)
batchType=c(rep(1,ncol(tcgaOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType, par.prior=TRUE)
tcgaOut=outTab[,colnames(tcgaOut)]
#tcgaOut[tcgaOut<0]=0
geoOut=outTab[,colnames(geoOut)]
#geoOut[geoOut<0]=0

# 保存处理后的数据
tcgaTab=rbind(ID=colnames(tcgaOut), tcgaOut)#rbind 函数用于将两个数据框或矩阵按行合并，即将它们上下堆叠在一起。ID=colnames(tcgaOut) 创建了一个行名为 "ID" 的向量，该向量的值是 tcgaOut 数据矩阵的列名（即样本的标识）。
write.table(tcgaTab, file="TCGA.normalize.txt", sep="\t", quote=F, col.names=F)
geoTab=rbind(ID=colnames(geoOut), geoOut)
write.table(geoTab,file="GEO.normalize.txt",sep="\t",quote=F,col.names=F)

# 读取感兴趣基因列表并取交集
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(tcgaOut))#intersect() 函数用于找到两个向量（或集合）之间的交集，即它会返回两者都包含的元素。
tcgaShareExp=tcgaOut[sameGene,]
geoShareExp=geoOut[sameGene,]

# 保存共同基因的处理后数据
tcgaShareExp=rbind(ID=colnames(tcgaShareExp),tcgaShareExp)
write.table(tcgaShareExp,file="TCGA.share.txt",sep="\t",quote=F,col.names=F)
geoShareExp=rbind(ID=colnames(geoShareExp),geoShareExp)
write.table(geoShareExp,file="GEO.share.txt",sep="\t",quote=F,col.names=F)


######??????学??: https://www.biowolf.cn/
######?纬?链??1: https://shop119322454.taobao.com
######?纬?链??2: https://ke.biowolf.cn
######?纬?链??3: https://ke.biowolf.cn/mobile
######?饪???师????: seqbio@foxmail.com
######?饪???师微??: eduBio

