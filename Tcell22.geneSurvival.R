#install.packages("survival")
#install.packages("survminer")


#???ð?
library(survival)
library(survminer)

inputFile="risk.TCGA.txt"      #?????ļ?
setwd("/Users/apple/Desktop/brca t cell/22.geneSur")     #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[,1:(ncol(rt)-2)]

#?Ի???????ѭ??
outTab=data.frame()
for(gene in colnames(rt)[3:ncol(rt)]){
	if(sd(rt[,gene])<0.001){next}
	data=rt[,c("futime", "fustat", gene)]
	colnames(data)=c("futime", "fustat", "gene")
	
	#??ȡ????cutoff
	res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
	res.cat=surv_categorize(res.cut)
	fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
	
	#?Ƚϸߵͱ???????????????
	diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
	pValue=1-pchisq(diff$chisq, df=1)
	outVector=cbind(gene, res.cut$cutpoint[1], pValue)
	outTab=rbind(outTab,outVector)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
		
	#????????????
	surPlot=ggsurvplot(fit,
					data=res.cat,
					conf.int=F,
					pval=pValue,
					pval.size=6,
					legend.title=gene,
					legend.labs=c("high","low"),
					xlab="Time(years)",
					break.time.by=2,
					palette=c("red", "blue"),
					risk.table=F,
					risk.table.title="",
					risk.table.height=.25)
	
	#????ͼ??
	pdf(file=paste0("Survival.",gene,".pdf"), width=4.5, height=4, onefile=FALSE)
	print(surPlot)
	dev.off()
}
#???????????ƺ?pֵ?ı????ļ?
write.table(outTab,file="survival.result.txt",sep="\t",row.names=F,quote=F)


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

