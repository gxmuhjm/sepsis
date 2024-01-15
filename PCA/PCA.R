

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")



library(limma)
library(ggplot2)

clusterFile="cluster.txt"      
setwd("D:\\bioinfor\\sepsis\\sepsis\\16.PCA")     


rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])

data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)


bioCol=c("#77BFA7","#F88D67","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]



veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), Cluster=g))
}


pdf(file="PCA.pdf", width=6.5, height=5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
	scale_colour_manual(name="Cluster", values =crgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()




