#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#???ð?
library(limma)
library(sva)
# ����TCGA��GEO������������ļ����ļ���
tcgaExpFile="merge.txt"       #TCGA?????????ļ?
geoExpFile="selectedtnbc_expression.txt"     #GEO?????????ļ?
geneFile="interGene.txt"       #?????????б??ļ?
setwd("D:\\Nature\\CD8\\CD8\\16.intersect")     #???ù???Ŀ¼

# ��ȡ�ʹ���TCGA�����������
rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1] 这个命令会将行名消失，可能是读取文件问题，保证行名是基因名就�?
exp=rt[,2:ncol(rt)]#[, 2:ncol(rt)] ��ʾѡ�� rt ���ݿ��е������У��յĵ�һ�����ţ��ʹӵ�2�е����һ�е�������
dimnames=list(rownames(exp),colnames(exp))#rownames(exp) ���� exp ��������������������ı�ǩ; colnames(exp) ���� exp ��������������������ı�ǩ
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)#avereps(tcga) ��������Ⲣȥ������ tcga �е��ظ��У��������ظ��е�ƽ��ֵ��������ÿ��Ψһ���ж�ֻ�����һ�Σ����ظ����лᱻ�ϲ�Ϊһ��ƽ��ֵ
tcga=log2(tcga+1)#log2(tcga+1) �Ƕ� tcga �����ÿ��Ԫ��ִ�ж����任�Ĳ�����log2 ��ʾ��2Ϊ�׵Ķ����任�����ǳ��õĶ����任֮һ��tcga+1 ������ÿ��Ԫ���ϼ�1���Ա������0ֵ����Ϊ��������Ӧ����0������

# ��ȡ�ʹ���GEO�����������
rt=read.table(geoExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]

exp=rt[,3:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)

# ��GEO���ݽ��д���������log2ת���͹淶��
qx=as.numeric(quantile(geo, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))#c(0, 0.25, 0.5, 0.75, 0.99, 1.0)����ʾҪ����ķ�λ��  v
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )#qx[5]>100 ����99�ٷ�λ����qx[5]���Ƿ����100��qx[6]-qx[1]>50 ����99�ٷ�λ���͵�0�ٷ�λ��֮��Ĳ��Ƿ����50��qx[2]>0 ����һ�ķ�λ����qx[2]���Ƿ����0
if(LogC){
    geo[geo<0]=0# �� geo ���ݾ�����С��0��ֵ��Ϊ0���������ݽ��н�β����ȷ���������ݶ���С��0
    geo=log2(geo+1)}
geo=normalizeBetweenArrays(geo)

# ȡTCGA��GEO���ݹ�ͬ�Ļ��򣬽��н�һ������
sameGene=intersect(row.names(tcga),row.names(geo))
tcgaOut=tcga[sameGene,]
geoOut=geo[sameGene,]

# �ϲ�TCGA��GEO����
all=cbind(tcgaOut,geoOut)
batchType=c(rep(1,ncol(tcgaOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType, par.prior=TRUE)
tcgaOut=outTab[,colnames(tcgaOut)]
#tcgaOut[tcgaOut<0]=0
geoOut=outTab[,colnames(geoOut)]
#geoOut[geoOut<0]=0

# ���洦���������
tcgaTab=rbind(ID=colnames(tcgaOut), tcgaOut)#rbind �������ڽ��������ݿ������кϲ��������������¶ѵ���һ��ID=colnames(tcgaOut) ������һ������Ϊ "ID" ����������������ֵ�� tcgaOut ���ݾ�����������������ı�ʶ����
write.table(tcgaTab, file="TCGA.normalize.txt", sep="\t", quote=F, col.names=F)
geoTab=rbind(ID=colnames(geoOut), geoOut)
write.table(geoTab,file="GEO.normalize.txt",sep="\t",quote=F,col.names=F)

# ��ȡ����Ȥ�����б���ȡ����
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(tcgaOut))#intersect() ���������ҵ������������򼯺ϣ�֮��Ľ����������᷵�����߶�������Ԫ�ء�
tcgaShareExp=tcgaOut[sameGene,]
geoShareExp=geoOut[sameGene,]

# ���湲ͬ����Ĵ���������
tcgaShareExp=rbind(ID=colnames(tcgaShareExp),tcgaShareExp)
write.table(tcgaShareExp,file="TCGA.share.txt",sep="\t",quote=F,col.names=F)
geoShareExp=rbind(ID=colnames(geoShareExp),geoShareExp)
write.table(geoShareExp,file="GEO.share.txt",sep="\t",quote=F,col.names=F)


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?�???ʦ????: seqbio@foxmail.com
######?�???ʦ΢??: eduBio
