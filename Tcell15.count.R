countFilter=20      #???????????еļ???ֵ???ڻ????? 20
inputFile="string_interactions_short.tsv"            #?????ļ?(??????ϵ?ļ?)
setwd("D:\\Nature\\CD8\\CD8\\15.PPI")     #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=FALSE, comment.char="")
tb=table(c(as.vector(rt[,1]), as.vector(rt[,2])))#tb ????һ??????��?еı????????е?һ?а?????һ?????ݵ?Ψһֵ???ڶ??а????ڶ??????ݵ?Ψһֵ??????ʾ???ǵĳ??ִ???
tb=sort(tb, decreasing=T)#??֮ǰ?????ı??? tb ???????????????????ս??򣨵ݼ?˳???????????б???????

#???????????Ļ???
outTab=as.data.frame(tb)#as.data.frame() ??R?????е?һ???????????ڽ????????????ͣ??????񡢾????ȣ?ת??Ϊ???ݿ?
outTab=outTab[outTab[,2]>=countFilter,]#[outTab[, 2] >= countFilter, ] ??һ??????ɸѡ??????????Ŀ????ѡ?????ڶ????е?ֵ???ڻ????? countFilter ???С?
colnames(outTab)=c("Gene","Count")#colnames(outTab) = c("Gene", "Count") ?ⲿ?ִ??뽫֮ǰ??ȡ???????滻Ϊ?µ?????????????һ?е????ƽ???????Ϊ "Gene"???ڶ??е????ƽ???????Ϊ "Count"
write.table(outTab,file="hubGene.txt",sep="\t",quote=F,row.names=F)
showNum=nrow(outTab)#nrow() ??R?????е?һ???????????ڻ?ȡ???ݿ????????????????????ݿ??еļ?¼??��??,ͨ??ִ?????д??룬???????? showNum ??��?еõ? outTab ???ݿ?????????
n=as.matrix(tb)[1:showNum,]#as.matrix(tb) ?????? tb ת??Ϊһ??????????????Ϊ?????;????ǲ?ͬ?????ݽṹ????ʱ??Ҫ??????ת??Ϊ??ͬ?????ݽṹ?Խ????ض??Ĳ???

#????һ??ˮƽ????ͼ??horizontal bar plot??
pdf(file="barplot.pdf", width=8, height=8)
par(mar=c(5,6,1,2),xpd=T)
bar=barplot(n,horiz=TRUE,col="skyblue",names=FALSE,xlim=c(0,ceiling(max(n)/5)*5),xlab="Number of adjacent nodes")
text(x=n*0.935,y=bar,n)     #??ͼ???м????ڽӽڵ???Ŀ
text(x=-0.2,y=bar,label=names(n),xpd=T,pos=2)     #??ͼ???м??ϻ?????????
dev.off()


######??????