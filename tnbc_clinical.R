rm(list=ls())
options(stringsAsFactors = F)
setwd("D:\\Nature\\CD8\\CD8\\06.getClinical")
merge<-read.table("merge.txt", header=T, sep="\t")
xlsx_1<-read.xlsx("D:\\Nature\\CD8\\CD8\\06.getClinical\\clinical.xlsx", sheet = "clinical")
file.exists("clinical.xlsx")
format_from_ext("clinical.xlsx")
clinical<-read.table("clinical.txt",header = T,sep = "\t")
rt<-t(merge)
a<-rownames(rt)
a<-gsub("\\.", "-", a, ignore.case = FALSE, perl = FALSE,
        fixed = FALSE, useBytes = FALSE)
a=as.data.frame(a)
rownames(clinical)<-clinical[,1]
match<-clinical[rownames(clinical)%in%a]
match<-subset(clinical,rownames(clinical)%in%a)
write.table(match, file="tnbc_clinical.txt", sep="\t", row.names=FALSE)

match<-merge.data.frame(b,clinical,by=intersect(rownames(b),rownames(clinical)))
