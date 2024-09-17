##Reading data
library(WGCNA)
PaData = read.csv("Pa-gene-exp-TPM-20au-21sp.csv", header = TRUE, row.names = 1)
dim(PaData)
names(PaData)
PaData0=as.data.frame(t(log2(PaData+1)))
gsg = goodSamplesGenes(PaData0, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(PaData0), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 200, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)
keepSamples = (clust==1)
PaData = PaData0[keepSamples, ]
nGenes = ncol(PaData)
nSamples = nrow(PaData)
traitData = read.csv("Pa-DHW-Depth-symbD-Shannon-XS.csv")
dim(traitData)
names(traitData)
allTraits=traitData
dim(allTraits)
names(allTraits)
PaSample=rownames(PaData)
traitRows=match(PaSample, allTraits$Sample)
datTraits=allTraits[traitRows, -1]
rownames(datTraits)=allTraits[traitRows, 1]
collectGarbage()
sampleTree2 = hclust(dist(PaData), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
tiff(filename = "Fig. S4A.tiff", width = 15, height = 6, units = "in", bg = "white", res=600)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
while (!is.null(dev.list()))  dev.off()

##network construction
enableWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(PaData, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
net = blockwiseModules(PaData, power = 14, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "Pa-XS-20au-21sp-TOM", verbose = 3)
table(net$colors)
#     0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20 
# 11532  3902  1497  1003   625   598   287   277   261   239   231   189   183   133   131   130    70    52    43    40    34 
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

## association with traits
nGenes = ncol(PaData)
nSamples = nrow(PaData)
MEs0 = moduleEigengenes(PaData, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
tiff(filename = "Module-trait-relationships.tiff", width = 9, height = 9, units = "in", bg = "white", res=600)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(25), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
while (!is.null(dev.list()))  dev.off()

DHW = as.data.frame(datTraits$DHW)
names(DHW) = "DHW"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(PaData, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(PaData, DHW, use = "p")) 
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(DHW), sep="")
names(GSPvalue) = paste("p.GS.", names(DHW), sep="")

module = "green"
column = match(module, modNames)
moduleGenes = moduleColors==module
tiff(filename = "Fig. S4B.tiff", width = 6, height = 6, units = "in", bg = "white", res=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DHW",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
while (!is.null(dev.list()))  dev.off()

module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module
tiff(filename = "Fig. S4C.tiff", width = 6, height = 6, units = "in", bg = "white", res=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DHW",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
while (!is.null(dev.list()))  dev.off()

module = "tan"
column = match(module, modNames)
moduleGenes = moduleColors==module
tiff(filename = "Fig. S4D.tiff", width = 6, height = 6, units = "in", bg = "white", res=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DHW",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
while (!is.null(dev.list()))  dev.off()

module = "midnightblue"
column = match(module, modNames)
moduleGenes = moduleColors==module
tiff(filename = "Fig. S4E.tiff", width = 6, height = 6, units = "in", bg = "white", res=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DHW",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
while (!is.null(dev.list()))  dev.off()

module = "turquoise"
column = match(module, modNames)
moduleGenes = moduleColors==module
tiff(filename = "Fig. S4F.tiff", width = 6, height = 6, units = "in", bg = "white", res=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DHW",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
while (!is.null(dev.list()))  dev.off()

module = "lightcyan"
column = match(module, modNames)
moduleGenes = moduleColors==module
tiff(filename = "Fig. 4D.tiff", width = 6, height = 6, units = "in", bg = "white", res=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DHW",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "cyan")
while (!is.null(dev.list()))  dev.off()


## GO overrepresentation plot
library(stringr)
library(ggsci)
library(ggplot2)
GO_info = read.csv("Model-positive-GO-overrepresentation.csv", header = TRUE)
Pisitive <- ggplot(GO_info, aes(x=GeneRatio, y=Desc, color=Model)) + geom_point(aes(size=log10(Count)))+scale_y_discrete(labels=function(y) str_wrap(y, width=50), limits=c("sensory perception", "sensory perception of chemical stimulus", "cognition", "regulation of vasoconstriction", "response to acid", "regulation of tube size", "response to hyperoxia", "response to stimulus", "regulation of blood vessel size", "lipoxygenase pathway", "regulation of cellular metabolic process", "regulation of nitrogen compound metabolic process", "regulation of gene expression", "regulation of cellular component organization", "regulation of primary metabolic process", "regulation of signaling pathway", "mRNA processing", "regulation of nucleobase, nucleoside, nucleotide and nucleic acid metabolic process", "regulation of macromolecule biosynthetic process", "regulation of cellular process", "peptide catabolic process"))+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 4B.pdf", Pisitive, width = 6, height = 6, units = "in")
GO_info = read.csv("Model-negative-GO-overrepresentation.csv", header = TRUE)
Negative <- ggplot(GO_info, aes(x=GeneRatio, y=Desc, color=Model)) + geom_point(aes(size=log10(Count)))+scale_y_discrete(labels=function(y) str_wrap(y, width=50), limits=c("translation", "peptide biosynthetic process", "cotranslational protein targeting to membrane", "SRP-dependent cotranslational protein targeting to membrane", "amide biosynthetic process", "protein targeting to ER", "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay", "peptide metabolic process", "protein targeting to membrane", "translational initiation", "DNA replication", "DNA metabolic process", "small molecule metabolic process", "DNA-dependent DNA replication", "cellular nitrogen compound metabolic process", "heterocycle metabolic process", "DNA repair", "cellular biosynthetic process", "cellular aromatic compound metabolic process", "biosynthetic process"))+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 4C.pdf", Negative, width = 6, height = 6, units = "in")

