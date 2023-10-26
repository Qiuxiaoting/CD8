#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)      #引用包
setwd("C:\\Users\\lexb\\Desktop\\Tcell\\27.maftools")     #设置工作目录

#读取风险文件，得到瀑布图的注释信息
risk=read.table("risk.TCGA.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#读取基因突变的文件
geneNum=20     #设置展示基因的数目
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#定义注释的颜色
ann_colors=list()
col=c("#0088FF", "#FF5555")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

#绘制低风险组瀑布图
pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#绘制高风险组瀑布图
pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio

