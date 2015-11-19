#-*-r-*-
library(WGCNA);
options(stringsAsFactors = FALSE);
ccin <- read.table("ccin.exp",header=T,sep="\t");
datExpr = as.data.frame(t(ccin[,-c(1:3)]));
names(datExpr) = ccin$IU_GENE;
rownames(datExpr) = names(ccin)[-c(1:3)];

#net = blockwiseModules(datExpr, power = 6,
#  minModuleSize = 30,reassignThreshold = 0,
#  mergeCutHeight = 0.25,numericLabels = TRUE,
#  pamRespectsDendro = FALSE,saveTOMs = TRUE,
#  saveTOMFileBase = "Ccin_TOM",verbose = 3)

# Choose a set of soft-thresholding powers
#powers = c(c(1:10),seq(from = 12, to=20, by=2))
# Call the network topology analysis function
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
#pdf("softPower.pdf",width=9,height=5);

#par(mfrow = c(1,2));
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     xlab="Soft Threshold (power)",
#     ylab="Scale Free Topology Model Fit,signed R^2",
#     type="n",main = paste("Scale independence"));
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
## this line corresponds to using an R^2 cut-off of h
#abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
#     xlab="Soft Threshold (power)",ylab="Mean Connectivity",
#     type="n",main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 6

ADJ1=abs(cor(datExpr,use="p"))^softPower 
#k=softConnectivity(datE=datExpr,power=softPower) 
k=as.vector(apply(ADJ1,2,sum, na.rm=T)) 

dissADJ=1-ADJ1 
dissTOM=TOMdist(ADJ1) 

hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))

moduleColors=colorDynamicHybridTOM
# Select module probes
modules;

modules = c("turquoise","blue")
inModule = is.finite(match(moduleColors, modules)); 

probes = names(datExpr)
modProbes = probes[inModule];
modGenes = stress$Annotation[match(modProbes,stress$AffyID)]

# Select the corresponding Topological Overlap 

modTOM = ADJ1[inModule, inModule]; 
dimnames(modTOM) = list(modProbes, modProbes) 

cyt = exportNetworkToCytoscape(modTOM, 
edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), 
nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""), 
 weighted = TRUE, 
 threshold = 0.10, 
 nodeNames = modProbes, 
 altNodeNames = modGenes,
 nodeAttr = moduleColors[inModule] );
