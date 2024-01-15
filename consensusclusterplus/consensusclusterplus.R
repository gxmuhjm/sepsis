

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)      
expFile="diffGeneExp.txt"          
workDir="D:\\bioinfor\\sepsis\\sepsis\\14.cluster_K=9\\GSE66099"     
setwd(workDir)      


data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)


group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="Treat"]


maxK=9     
results=ConsensusClusterPlus(data,
              maxK=maxK,  
              reps=50,   
              pItem=0.8,  
              pFeature=1,  
              title=workDir,  
              clusterAlg="km",   
              distance="euclidean",  
              seed=123456,  
              plot="png")


calcICL(results, title="consensusScore", plot="png")


clusterNum=2       
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
cluster$Cluster=paste0("C", cluster$Cluster)
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)







