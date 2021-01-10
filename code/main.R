#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(WGCNA)
getwd()
joint_data = read.csv("joint_normalized.csv", row.names = 1)

datExpr = t(joint_data)

#=====================================================================================
#
#  Code chunk 2 : choose the power
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
cex1 = 0.9

pdf("soft_threshold.pdf", width=9, height=6)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft threshold (power)",ylab="Scale-free topology model fit, signed R²",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="blue")
dev.off() 


#=====================================================================================
#
#  Code chunk 3 : cal net
#                 in this part, we need to pay attention to parameters which are
#                 'power', 'minModuliSize', and 'mergeCutHeight'
#
#=====================================================================================
allowWGCNAThreads() #enableWGCNAThreads()

mergeCutHeight = 0.25
minModuleSize = 10
reassignThreshold = 0
power = 6
net <- blockwiseModules(datExpr, power = power,
                         TOMType = "unsigned", minModuleSize = minModuleSize,     #30,
                         reassignThreshold = reassignThreshold, mergeCutHeight =  mergeCutHeight,    # 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         saveTOMFileBase = "femaleMouseTOM")
netcolors = net$colors
matrixdata<- data.frame(cbind(t(datExpr), netcolors))


#=====================================================================================
#
#  Code chunk 4 : plot clustring
#
#=====================================================================================

mergedColors = labels2colors(net$colors)


pdf("dendrogram.pdf", width=9, height=6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = sprintf("Power = %d, minModuleSize = %d, mergeCutHeight = %.4f",
                                   power, minModuleSize, mergeCutHeight))
dev.off() 

geneCharVector <- matrix(0, nrow = 0, ncol = length(unique(netcolors))-1)
print("unique netcolors: ")
print(unique(netcolors))

eigengene <- matrix(0, nrow = length(unique(netcolors))-1, ncol = dim(t(datExpr))[2]) # Clusters * Samples
for (i in 1: (length(unique(netcolors))-1) ){
    geneID <- which(matrixdata$netcolors == i)
    # ===== Calculate Eigengene Start
    X <- t(datExpr)[geneID,]
    mu <- rowMeans(X)
    stddev <- rowSds(as.matrix(X), na.rm=TRUE) # standard deviation with 1/(n-1)
    #normalize X:
    XNorm <- sweep(X,1,mu)
    XNorm <- apply(XNorm, 2, function(x) x/stddev)
    SVD <- svd(XNorm, LINPACK = FALSE)
    eigengene[i,] <- t(SVD$v[,1])
    # ===== Calculate Eigengene Finished
    geneChar <- c(toString(i), labels2colors(i), rownames(matrixdata)[geneID])
    geneCharVector[i] <- list(geneChar)
}
colnames(eigengene) <- rownames(datExpr)

## Compute maximum length
max.length <- max(sapply(geneCharVector, length))
## Add NA values to list elements
geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
## Rbind

geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))
write.csv(geneCharVector2, 'coexpression_modules.csv')


for (i in 1:length(geneCharVector)){
  if ('mRNA*CLIC4|25932_calculated' %in% geneCharVector[i][[1]]){
    clusterid = i
  }
}
print(clusterid)
cluster = geneCharVector[clusterid][[1]]
cluster = cluster[3:length(cluster)]

targetcluster <- matrix(0, nrow = length(cluster), ncol = 2)
for (i in 1:length(cluster)){
  splitval = strsplit(cluster[i], '*', fixed = TRUE)[[1]]
  type = splitval[1]
  gene = strsplit(splitval[2], '|', fixed = TRUE)[[1]][1]
  targetcluster[i,1] = type
  targetcluster[i,2] = gene
}

targetcluster = data.frame(targetcluster)
colnames(targetcluster) = c('type', 'gene name')
write.csv(targetcluster, sprintf('cluster_%d.csv', clusterid))

