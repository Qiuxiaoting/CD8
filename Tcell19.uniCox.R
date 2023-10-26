#install.packages('survival')
#install.packages("survminer")


#???冒?
library(survival)
library(survminer)#提供了用于生存分析（survival analysis）和生存曲线（survival curves）可视化的工具和函数。生存分析通常用于研究事件发生的时间，如患病、死亡或其他事件，以及与这些事件相关的因素。

coxPfilter=0.05                 # 设置显著性水平
inputFile="TCGA.expTime.txt"     #?????募?
setwd("D:\\Nature\\CD8\\CD8\\19.uniCox")     #???霉???目录

#??取?????募?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365# 将存活时间转换为年

# 创建一个空的数据框用于存储结果
outTab=data.frame()
sigGenes=c("futime","fustat")# 在生存分析（survival analysis）中，fustat 通常代表着生存状态或事件状态。具体来说，fustat 是一个二进制变量，通常使用以下方式编码：1: 代表"事件"或"失败"，通常是指研究中关注的不良事件或生存终点，例如死亡、疾病复发等。这表示在特定时间点发生了事件。0: 代表"无事件"或"存活"，这表示在特定时间点没有发生事件，或者个体仍然存活。
#根据p值的结果，可以确定哪些变量在此生存分析中可能是重要的预测因子。
for(i in colnames(rt[,3:ncol(rt)])){#通过for循环，逐一遍历从第3列到数据集rt的最后一列，ncol(rt)代表数据集rt的总列数。
	#模型的目标是预测生存时间
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)#Surv(futime, fustat) ~ rt[,i] 部分是回归模型的设定。Surv(futime, fustat) 定义了生存时间（或存活时间）和事件状态（成功或失败）的关系。rt[,i] 是独立变量（解释变量），它代表数据集中的一个特定列。Surv(futime, fustat) ~ rt[,i] 部分是回归模型的设定。Surv(futime, fustat) 定义了生存时间（或存活时间）和事件状态（成功或失败）的关系。rt[,i] 是独立变量（解释变量），它代表数据集中的一个特定列。
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	# 检查p值，如果小于阈值，则添加到结果中
	if(coxP<coxPfilter){#这个条件语句检查 coxP，也就是Cox回归模型中独立变量与事件发生风险的p值，是否小于预设的 coxPfilter 阈值。coxPfilter 通常是一个显著性水平，例如 0.05，用于确定变量与风险的关系是否显著。
	    sigGenes=c(sigGenes,i)#如果p值小于阈值，将当前的变量i添加到名为sigGenes的向量或列表中，以标记这个变量在统计上是显著的。
		outTab=rbind(outTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			        )
	}
}

#?????????氐慕???
write.table(outTab,file="TCGA.uniCox.txt",sep="\t",row.names=F,quote=F)

#??取???????????????谋???????
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="TCGA.uniSigExp.txt",sep="\t",row.names=F,quote=F)


# 定义一个函数用于创建森林图
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  # 读取Cox回归结果文件
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")# "%.3f" 是格式说明符，表示要将数值格式化为带有三位小数的浮点数。rt$"HR" 中的数值将按这个格式转换为字符串。
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
		
	# 绘制森林图的左侧部分
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		
	# 绘制森林图的右侧部分
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))# 设置x轴范围
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])# 根据Hazard ratio的值设置颜色
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)# 绘制点
	axis(1)# 添加x轴
	dev.off()
}

# 调用函数生成森林图
bioForest(coxFile="TCGA.uniCox.txt", forestFile="forest.pdf", forestCol=c("red","green"))




