#install.packages("VennDiagram")


library(VennDiagram)      #???ð?
setwd("D:\\Nature\\CD8\\CD8\\14.venn")     #???ù???Ŀ¼
geneList=list()

#??ȡTϸ???????????Ľ????ļ?
rt=read.table("CD8T_diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,5])              #�����ݿ������еĵ���������ת��Ϊһ��������������������洢����Ϊ geneNames �ı����С��������ͨ��������������ƻ��ʶ������� geneNames ���ܴ����������Ƶ�������
geneNames=gsub("^ | $","",geneNames)     #������ geneNames �е�ÿ��Ԫ�أ�������Կո�ͷ���Կո��β������Щ�ո�ɾ�������յõ��������е�Ԫ�ؽ�û�пո�ͷ���β�Ĳ��֡�
uniqGene=unique(geneNames)               #�ҳ� geneNames �����в��ظ���Ԫ�أ������Ǵ洢�� uniqGene �У����� uniqGene �е�ÿ��Ԫ�ض��Ƕ�һ�޶��ģ���������ظ�
geneList[["CD8T_diff"]]=uniqGene         #??Tϸ???Ĳ??????????�??????б???

#??ȡ???߻??????б??ļ?
rt=read.table("gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #??ȡ???߻?????????
geneNames=gsub("^ | $","",geneNames)     #ȥ????????β?Ŀո?
uniqGene=unique(geneNames)               #?????߻???ȡunique
geneList[["Immune_gene"]]=uniqGene       #?????????صĻ??????�??????б???

##ʹ�ø����Ļ����б��򼯺ϴ���һ��ά��ͼ����ָ��ͼ�ε���ɫ����ǩλ�ú��������������ɵ�ά��ͼ�������ڿ��ӻ���ͬ���򼯺�֮����ص���ϵ��
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)#ͨ��ִ�� grid.draw(venn.plot)�������Rͼ���豸�Ͽ���һ��ά��ͼ�����ͼ�α�ʾ��ͬ����֮����ص���ϵ
dev.off()

#geneList ���������б���A��B��C������ô���Ȼ�����б�A���б�B�Ľ�����Ȼ���ټ�������������б�C�Ľ��������յõ������������б��Ľ��������ԣ�interGenes ������ geneList �������б��Ľ��������ʾ���а������������б��ж����ڵĻ����Ԫ�ء����ֲ���ͨ����������ѧ�����ݷ����У��Բ��Ҷ�����ݼ�֮��Ĺ�ͬԪ�ػ򽻼���
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGene.txt", sep="\t", quote=F, col.names=F, row.names=F)


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?�???ʦ????: seqbio@foxmail.com
######?�???ʦ΢??: seqBio

