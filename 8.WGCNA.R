setwd("Work path") 

library(limma)
library(WGCNA)

data=read.csv("mRNA_TPM_new", header=T, row.names = 1)
# Convert to matrix.
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

# Modify sample names, keeping only the first 12 characters.
colnames(data)=substr(colnames(data),1,12)
colnames(data) = gsub('[.]', '-', colnames(data))
# Remove lowly expressed genes.
data=data[rowMeans(data)>0.5,]

# Read in clinical data, indicating whether the sample has cancer.
cli=read.csv("state.csv", header = TRUE, row.names = 1)

# Common samples.
sameSample=intersect(row.names(cli), colnames(data))
# Extract the common samples separately.
cli = cli[sameSample,]
data = data[,sameSample]
dim(data)

# Select the top 50% of genes with the highest variability for WGCNA.
selectGenes=names(tail(sort(apply(data,1,sd)), n=round(nrow(data)*0.5)))
data=data[selectGenes,]
dim(data)
# Transpose.
datExpr0=t(data)

# Check for missing values and outliers.
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster all samples and observe if there are any outliers or anomalies.
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot 1.
pdf(file = "1_sample_cluster.1.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

cutHeight = 100000
# Plot 2.
pdf(file = "1_sample_cluster.2.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cutHeight, col = "red")
dev.off()

# Cut.
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]
dim(data)
dim(datExpr0)

# Construct an automated network and detect modules.
# Scatter plot of power values, select the soft threshold.
enableWGCNAThreads()   
powers = c(1:20)       
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="2_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))

cex1 = 0.8
# Fit the index and scatter plot of power values, with the scale-free topology fitting index.
# Generally, select the first value that exceeds 0.9.
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=cex1,col="red")

# Scatter plot of average connectivity versus power values.
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Manually set the soft threshold.
softPower = 4
adjacency = adjacency(datExpr0, power = softPower)
softPower
# If there is no suitable soft threshold, you can choose based on the following conditions:
if (is.na(power)){ 
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18), 
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16), 
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14), 
                               ifelse(type == "unsigned", 6, 12))        
                 ) 
  ) 
}

# TOM matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Gene clustering.
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="3_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Dynamic tree cutting module identification.
minModuleSize=100 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="4_Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Identify clusters of similar modules.
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="5_Clustering_module.pdf",width=7,height=7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.3  # Cut height.
abline(h=MEDissThres, col = "red")
dev.off()

# Merge similar modules.
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="6_merged_dynamic.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Read in sample status information.
cli2=read.csv("state.csv", header=T, check.names=F, row.names=1)

sameSample2=intersect(row.names(cli), rownames(MEs))
# Extract the common samples separately.
MEs=MEs[sameSample2,]
datTraits=cli[sameSample2,]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# Correlation analysis.
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="7_Module_trait.pdf", width=6.5, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Module to which the gene belongs.
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "module_all.txt",sep="\t",row.names=F,quote=F)

# Output the genes for each module.
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

# Run the following code to visualize GS and MM.
module = "blue" # Module replacement operation.
Selectedclinical = "Normal"
Selectedclinical2 = "Tumor"

# Take MEturquoise and Grade as an example.
Selectedclinical = as.data.frame(datTraits[,Selectedclinical]);
names(Selectedclinical) = "Selectedclinical";
modNames = substring(names(MEs), 3)
datExpr1 = datExpr0[rownames(MEs),]
geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr1, Selectedclinical, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Selectedclinical), sep="");
names(GSPvalue) = paste("p.GS.", names(Selectedclinical), sep="");

# Plot the graph.
column = match(module, modNames)
moduleGenes = moduleColors==module
outPdf=paste(Selectedclinical2,"_", module, ".pdf",sep="")
pdf(file=outPdf,width=7,height=7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ", Selectedclinical2,sep=""),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# Output the hub genes of the module.
datMM=cbind(geneModuleMembership[,paste("MM", module, sep="")], geneTraitSignificance)
colnames(datMM)[1] = paste("MM", module, sep="")
# Importance and correlation.
geneSigFilter=0.1
moduleSigFilter=0.7
# Filter.
datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]
datMM=datMM[abs(datMM[,1])>moduleSigFilter,]
# Export.
write.table(row.names(datMM), file =paste0("hubGenes",module,".txt"),sep="\t",row.names=F,col.names=F,quote=F)

