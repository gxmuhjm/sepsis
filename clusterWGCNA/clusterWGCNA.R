

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("GO.db", "preprocessCore", "impute","limma"))

#install.packages(c("gplots", "matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
#install.packages("WGCNA")



library(limma)
library(gplots)
library(WGCNA)

expFile="normalize.txt"       
clusterFile="cluster.txt"     
setwd("D:\\bioinfor\\sepsis\\sepsis\\20.clusterWGCNA\\GSE66099")     


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
selectGenes=names(tail(sort(apply(data,1,sd)), n=round(nrow(data)*0.25)))
data=data[selectGenes,]


cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=data[,nameC1,drop=F]
dataC2=data[,nameC2,drop=F]
conCount=ncol(dataC1)
treatCount=ncol(dataC2)
data=cbind(dataC1, dataC2)
datExpr0=t(data)


gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "cluster01.sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

abline(h = 20000, col = "red")
dev.off()


clust=cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples=(clust==1)
datExpr0=datExpr0[keepSamples, ]



traitData=data.frame(C1=c(rep(1,conCount),rep(0,treatCount)),
                     C2=c(rep(0,conCount),rep(1,treatCount)))
row.names(traitData)=colnames(data)
fpkmSamples=rownames(datExpr0)
traitSamples=rownames(traitData)
sameSample=intersect(fpkmSamples,traitSamples)
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]



sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="cluster02.sample_heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()


enableWGCNAThreads() 
powers = c(1:20)       
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="cluster03.scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
ͼ
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
ͼ
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


sft 
softPower =sft$powerEstimate    
adjacency = adjacency(datExpr0, power = softPower)
softPower


TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="cluster04.gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



minModuleSize = 100      
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="cluster05.Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()



MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="cluster06.Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()


moduleColors=dynamicColors
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
select = sample(nGenes, size=1000)
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method="average")
selectColors = moduleColors[select]

plotDiss=selectTOM^softPower
diag(plotDiss)=NA
myheatcol = colorpanel(250, "red", "orange", "lemonchiffon")    
pdf(file="cluster07.TOMplot.pdf", width=7, height=7)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)
dev.off()



moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="cluster08.Module_trait.pdf", width=6.5, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3.5, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


trait="C2"
traitColumn=match(trait,traitNames)  
for (module in modNames){
	column = match(module, modNames)
	moduleGenes = moduleColors==module
	if (nrow(geneModuleMembership[moduleGenes,]) > 1){
	    outPdf=paste("cluster09.", trait, "_", module,".pdf",sep="")
	    pdf(file=outPdf,width=7,height=7)
	    par(mfrow = c(1,1))
	    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
	                         xlab = paste("Module Membership in", module, "module"),
	                         ylab = paste("Gene significance for ",trait),
	                         main = paste("Module membership vs. gene significance\n"),
	                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
	    abline(v=0.8,h=0.5,col="red")
	    dev.off()
	}
}



probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "cluster.GS_MM.xls",sep="\t",row.names=F)



for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("cluster.module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}


geneSigFilter=0.2        
moduleSigFilter=0.6      
cli="GS.C2"
for(mol in unique(geneInfo$moduleColor)){
	geneInfoMol=geneInfo[geneInfo$moduleColor==mol,]
	mmi=paste0("MM", mol)
	geneInfoMol2=geneInfoMol[((abs(geneInfoMol[,mmi])>moduleSigFilter) & (abs(geneInfoMol[,cli])>geneSigFilter)),]
	write.table(geneInfoMol2[,1], file =paste0("cluster.hubGenes_",mmi,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}




