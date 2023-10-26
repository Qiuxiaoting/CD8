#install.packages("VennDiagram")


library(VennDiagram)      #???冒?
setwd("D:\\Nature\\CD8\\CD8\\14.venn")     #???霉???目录
geneList=list()

#??取T细???????????慕????募?
rt=read.table("CD8T_diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,5])              #将数据框或矩阵中的第五列数据转换为一个向量，并将这个向量存储在名为 geneNames 的变量中。这个向量通常包含基因的名称或标识符，因此 geneNames 可能代表基因名称的向量。
geneNames=gsub("^ | $","",geneNames)     #将向量 geneNames 中的每个元素，如果其以空格开头或以空格结尾，则将这些空格删除，最终得到的向量中的元素将没有空格开头或结尾的部分。
uniqGene=unique(geneNames)               #找出 geneNames 中所有不重复的元素，将它们存储在 uniqGene 中，这样 uniqGene 中的每个元素都是独一无二的，不会出现重复
geneList[["CD8T_diff"]]=uniqGene         #??T细???牟??????????氲??????斜???

#??取???呋??????斜??募?
rt=read.table("gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #??取???呋?????????
geneNames=gsub("^ | $","",geneNames)     #去????????尾?目崭?
uniqGene=unique(geneNames)               #?????呋???取unique
geneList[["Immune_gene"]]=uniqGene       #?????????氐幕??????氲??????斜???

##使用给定的基因列表或集合创建一个维恩图，并指定图形的颜色、标签位置和其他参数。生成的维恩图可以用于可视化不同基因集合之间的重叠关系。
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)#通过执行 grid.draw(venn.plot)，你会在R图形设备上看到一个维恩图，这个图形表示不同集合之间的重叠关系
dev.off()

#geneList 包含三个列表（A、B、C），那么首先会计算列表A和列表B的交集，然后再计算这个交集与列表C的交集，最终得到包含了所有列表的交集。所以，interGenes 将包含 geneList 中所有列表的交集，这表示其中包含了在所有列表中都存在的基因或元素。这种操作通常用于生物学或数据分析中，以查找多个数据集之间的共同元素或交集。
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGene.txt", sep="\t", quote=F, col.names=F, row.names=F)


######??????学??: https://www.biowolf.cn/
######?纬?链??1: https://shop119322454.taobao.com
######?纬?链??2: https://ke.biowolf.cn
######?纬?链??3: https://ke.biowolf.cn/mobile
######?饪???师????: seqbio@foxmail.com
######?饪???师微??: seqBio

