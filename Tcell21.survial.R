#install.packages("survival")
#install.packages("survminer")


#???ð?
library(survival)
library(survminer)
setwd("/Users/apple/Desktop/brca t cell/21.survival")     #???ù???Ŀ¼

#?????????????ĺ???
bioSurvival=function(inputFile=null, outFile=null){
	#??ȡ?????ļ?
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#?Ƚϸߵͷ??????????????죬?õ??????????Ե?pvalue
	diff=survdiff(Surv(futime, fustat) ~ Risk, data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
		
	#????????????
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           ylab="Overall survival",
		           break.time.by = 2,
		           palette=c("red", "blue"),
		           risk.table=TRUE,
		           risk.table.title="",
		           risk.table.height=.25)
	#????ͼ??
	pdf(file=outFile, width=6, height=5, onefile=FALSE)
	print(surPlot)
	dev.off()
}

#???ú???,????????????
bioSurvival(inputFile="risk.TCGA.txt", outFile="survival.TCGA.pdf")
bioSurvival(inputFile="risk.GEO.txt", outFile="survival.GEO.pdf")


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

