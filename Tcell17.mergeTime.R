#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)               #???ð?
expFile="TCGA.share.txt"     #?????????ļ?
cliFile="time.txt"           #?????????ļ?
setwd("D:\\Nature\\CD8\\CD8\\17.mergeTime")      #???ù???Ŀ¼

# �� expFile �ж�ȡ����������ݲ����д���
rt=read.table(expFile, header=T, sep="\t", check.names=F) # ���ļ��ж�ȡ���ݣ�������ͷ��ʹ���Ʊ����ָ���
rt=as.matrix(rt)
rownames(rt)=rt[,1]# ���þ��������Ϊ��һ�е����ݡ�
exp=rt[,2:ncol(rt)]# ��ȡ����������ݣ�ȥ����һ�еı�ʶ��
dimnames=list(rownames(exp),colnames(exp))# �������ݵ�������������
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)# ������ת��Ϊ��ֵ����
data=avereps(data) # �����ݽ���ƽ��������

#去除正常样本 都是tnbc样本 试试不去�?
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2","1",group)
#data=data[,group==0]
#colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
#data=t(data)
#data=avereps(data) 

#??ȡ?????????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

# �ҵ����е��������ݲ�����ʱ�����ݺͻ���������ݡ�
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli, data)# ��ʱ�����ݺͻ���������ݰ��кϲ���

# ����������ʶ�У��������Ϻ������д���ļ� "TCGA.expTime.txt"��
out=cbind(id=row.names(out), out)# ����������ʶ�С�
write.table(out, file="TCGA.expTime.txt", sep="\t", row.names=F, quote=F)


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?�???ʦ????: seqbio@foxmail.com
######?�???ʦ΢??: eduBio

