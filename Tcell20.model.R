#install.packages("glmnet")
#install.packages("survival")


#???ð?
library("glmnet")
library("survival")

trainFile="TCGA.uniSigExp.txt"    #TCGA???ݿⵥ?????????????ı????ļ?
testFile="GEO.expTime.txt"        #GEO???????ݺ????????ݺϲ????ļ?
setwd("/Users/apple/Desktop/168.Tcell资料/20.model")     #???ù???Ŀ¼

#??ȡtrain???????ļ?
rt=read.table(trainFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime[rt$futime<=0]=0.003

#????lasso?ع?ģ??
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
#????lasso?ع???ͼ??
pdf("lasso.lambda.pdf")
plot(fit, xvar="lambda", label=TRUE)
dev.off()
#???ƽ?????֤??ͼ??
cvfit=cv.glmnet(x, y, family="cox", maxit=1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

#?ҵ???????֤??????С?ĵ㣬????????ģ?͹?ʽ
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

#????TCGA???ݿ??ķ????ļ?
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
Risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.TCGA.txt",sep="\t",quote=F,row.names=F)

#????GEO???ݿ??ķ????ļ?
rt=read.table(testFile, header=T, sep="\t", row.names=1)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO.txt",sep="\t",quote=F,row.names=F)


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

