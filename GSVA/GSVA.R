

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")



library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
library(devtools)
library(IRanges)

expFile="normalize.txt"              
clusterFile="cluster.txt"           
gmtFile="c2.cp.kegg.symbols.gmt"     
setwd("D:\\bioinfor\\sepsis\\sepsis\\18.GSVA")     


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)


group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]


geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())


ssgseaScore=gsva(data, geneSets, method='gsva')

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=ssgseaScore[,nameC1,drop=F]
dataC2=ssgseaScore[,nameC2,drop=F]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
data=cbind(dataC1, dataC2)
Type=c(rep("C1",conNum), rep("C2",treatNum))


outTab=data.frame()
for(i in row.names(data)){
	test=t.test(data[i,] ~ Type)
	pvalue=test$p.value
	t=test$statistic
	if(pvalue<0.05){
		Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
		outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
	}
}


termNum=10     
outTab=outTab[order(outTab$t),]
outTab=outTab[c(1:termNum,(nrow(outTab)-termNum):nrow(outTab)),]
pdf(file="barplot_KEGG.pdf", width=14, height=14)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
		palette=c("#BA6DEC", "#6BD9F5"), sort.val = "asc", sort.by.groups = T,
		rotate=TRUE, legend="right", title="",
		xlab="Term", ylab="t value of GSVA score, C2 vs C1",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()






