#install.packages('survival')
#install.packages("survminer")


#???ð?
library(survival)
library(survminer)#�ṩ���������������survival analysis�����������ߣ�survival curves�����ӻ��Ĺ��ߺͺ������������ͨ�������о��¼�������ʱ�䣬�综���������������¼����Լ�����Щ�¼���ص����ء�

coxPfilter=0.05                 # ����������ˮƽ
inputFile="TCGA.expTime.txt"     #?????ļ?
setwd("D:\\Nature\\CD8\\CD8\\19.uniCox")     #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365# �����ʱ��ת��Ϊ��

# ����һ���յ����ݿ����ڴ洢���
outTab=data.frame()
sigGenes=c("futime","fustat")# �����������survival analysis���У�fustat ͨ������������״̬���¼�״̬��������˵��fustat ��һ�������Ʊ�����ͨ��ʹ�����·�ʽ���룺1: ����"�¼�"��"ʧ��"��ͨ����ָ�о��й�ע�Ĳ����¼��������յ㣬�������������������ȡ����ʾ���ض�ʱ��㷢�����¼���0: ����"���¼�"��"���"�����ʾ���ض�ʱ���û�з����¼������߸�����Ȼ��
#����pֵ�Ľ��������ȷ����Щ�����ڴ���������п�������Ҫ��Ԥ�����ӡ�
for(i in colnames(rt[,3:ncol(rt)])){#ͨ��forѭ������һ�����ӵ�3�е����ݼ�rt�����һ�У�ncol(rt)�������ݼ�rt����������
	#ģ�͵�Ŀ����Ԥ������ʱ��
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)#Surv(futime, fustat) ~ rt[,i] �����ǻع�ģ�͵��趨��Surv(futime, fustat) ����������ʱ�䣨����ʱ�䣩���¼�״̬���ɹ���ʧ�ܣ��Ĺ�ϵ��rt[,i] �Ƕ������������ͱ����������������ݼ��е�һ���ض��С�Surv(futime, fustat) ~ rt[,i] �����ǻع�ģ�͵��趨��Surv(futime, fustat) ����������ʱ�䣨����ʱ�䣩���¼�״̬���ɹ���ʧ�ܣ��Ĺ�ϵ��rt[,i] �Ƕ������������ͱ����������������ݼ��е�һ���ض��С�
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	# ���pֵ�����С����ֵ�������ӵ������
	if(coxP<coxPfilter){#������������ coxP��Ҳ����Cox�ع�ģ���ж����������¼��������յ�pֵ���Ƿ�С��Ԥ��� coxPfilter ��ֵ��coxPfilter ͨ����һ��������ˮƽ������ 0.05������ȷ����������յĹ�ϵ�Ƿ�������
	    sigGenes=c(sigGenes,i)#���pֵС����ֵ������ǰ�ı���i���ӵ���ΪsigGenes���������б��У��Ա�����������ͳ�����������ġ�
		outTab=rbind(outTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			        )
	}
}

#?????????صĽ???
write.table(outTab,file="TCGA.uniCox.txt",sep="\t",row.names=F,quote=F)

#??ȡ???????????????ı???????
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="TCGA.uniSigExp.txt",sep="\t",row.names=F,quote=F)


# ����һ���������ڴ���ɭ��ͼ
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  # ��ȡCox�ع����ļ�
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")# "%.3f" �Ǹ�ʽ˵��������ʾҪ����ֵ��ʽ��Ϊ������λС���ĸ�������rt$"HR" �е���ֵ���������ʽת��Ϊ�ַ�����
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		

	height=nrow(rt)/12.5+5
	pdf(file=forestFile, width=7, height=height)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	# ����ɭ��ͼ����ಿ��
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		
	# ����ɭ��ͼ���Ҳಿ��
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))# ����x�᷶Χ
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])# ����Hazard ratio��ֵ������ɫ
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)# ���Ƶ�
	axis(1)# ����x��
	dev.off()
}

# ���ú�������ɭ��ͼ
bioForest(coxFile="TCGA.uniCox.txt", forestFile="forest.pdf", forestCol=c("red","green"))



