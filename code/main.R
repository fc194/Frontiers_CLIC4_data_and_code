library("RColorBrewer")
library("gplots")
library(ggplot2)
library(WGCNA)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rna<-read.table("73CN-AML-RNA-TCGA.csv",sep=',',header=T,row.name=1)
lbl = as.numeric(rna['CLIC4',] <= median(as.numeric(rna['CLIC4',])))
lbl[lbl==0] = 2
label = matrix(ncol=2,nrow=dim(rna)[2])
label[1:dim(rna)[2],1] = seq(1, dim(rna)[2])
label[1:dim(rna)[2],2] = lbl
label = data.frame(label)
write.table(label, "clic4_label.txt")


process_mRNA <- function(){
  #1 low; 2 high
  label<-read.table("clic4_label.txt",header=T)
  
  rna<-read.table("73CN-AML-RNA-TCGA.csv",sep=',',header=T,row.name=1)
  rna<-rna[-which(rowMeans(rna)<0.1),]
  rna<-log2(rna+1)
  
  low<-label[which(label[,2]==1),1]
  high<-label[which(label[,2]==2),1]
  
  rna_low<-rna[,low]
  rna_high<-rna[,high]
  
  rna_cmb<-cbind(rna_low,rna_high)
  rna_pval<-matrix(nrow=dim(rna_cmb)[1],ncol=2)
  for (i in 1:dim(rna_cmb)[1]){
    test_res<-t.test(rna_cmb[i,1:37],rna_cmb[i,38:73])
    log2fc<-log2(rowMeans(rna_cmb[i,38:73])/rowMeans(rna_cmb[i,1:37]))#high/low
    rna_pval[i,1]<-test_res$p.value
    rna_pval[i,2]<-log2fc
  }
  rownames(rna_pval)<-rownames(rna_cmb)
  colnames(rna_pval)<-c("pval","log2fc(high/low)")
  
  # padj<-p.adjust(rna_pval[,1],method="BY")
  # rna_pval_new<-cbind(rna_pval,padj)
  # colnames(rna_pval_new)<-c("pval","log2fc(high/low)","padj")
  
  sig<-rna_pval[which(rna_pval[,1]<0.01),]
  sig_label<-matrix(nrow=dim(rna_cmb)[1],ncol=1)
  sig_label[which(rna_pval[,1]>=0.01),1]<-"UnsigDEGs"
  sig_label[intersect (which(rna_pval[,1]<0.01),which(rna_pval[,2]>0)),1]<-"Up"
  sig_label[intersect (which(rna_pval[,1]<0.01),which(rna_pval[,2]<0)),1]<-"Down"
  rna_pval<-cbind(rna_pval,sig_label)
  colnames(rna_pval)<-c("pval","log2fc(high/low)","sig_label")
  sig_value<-rna_cmb[which(rna_pval[,1]<0.01),]
  
  write.table(rna_pval,"deg.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig,"sig_p_0.01.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(rna_cmb,"rna_cmb.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig_value,"sig_value_p_0.01.txt",quote=F,row.names=T,col.names=T,sep="\t")
  
  ###################KEGG
  a<-read.table("kegg_plot.txt",header=T,sep="\t")
  ggplot(data = a, aes(x = Term, y = -log10(a$PValue), fill = Label)) +  geom_bar(stat = "identity")+scale_fill_manual(values=c("#00A14B","#ED1C24")) +labs(x="",y="-log10(p-value)")+theme(axis.text.x=element_text(vjust = 1, hjust = 1,size=12),plot.margin=margin(l=30,r=5,t=5,b=5,unit="pt"))+ scale_x_discrete(limits=rev(a$Term))+coord_flip()
  
  ###################scatter
  a<-read.table("deg.txt",header=T,row.names=1,sep="\t")
  ggplot(a,aes(x=log2fc.high.low.,y=-log10(a$pval),colour=sig_label))+geom_point()+scale_color_manual(breaks=c("Up","Down","UnsigDEGs"),values=c("#00A14B","black","#ED1C24"))+labs(x="log2(Fold Change)",y="-log10(p-value)")+theme(panel.grid.major = element_blank())+theme_bw()+xlim(-4,4)
  
  ###################heatmap
  a<-read.table("sig_value_p_0.01.txt",header=T,row.names=1)
  g<-heatmap.2(as.matrix(a-rowMeans(a)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",margins=c(2,10),labCol=F,cexRow=0.1)
  
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
  
  ###################ANP32A
  rna_cmb<-read.table("rna_cmb.txt",row.names=1,header=T)
  clic4<-rna_cmb[which(rownames(rna_cmb) == "CLIC4|25932_calculated"),]
  plot(as.numeric(clic4[1,]),type="l")
}
process_mRNA()


process_miRNA <- function(){
  ### miRNA
  #1 low; 2 high
  label<-read.table("clic4_label.txt",header=T)
  
  rna<-read.table("73CN-AML-miRNA-TCGA.csv",sep=',',header=T,row.name=1)
  rna<-rna[-which(rowMeans(rna)<0.1),]
  rna<-log2(rna+1)
  
  low<-label[which(label[,2]==1),1]
  high<-label[which(label[,2]==2),1]
  rna_low<-rna[,low]
  rna_high<-rna[,high]
  
  rna_cmb<-cbind(rna_low,rna_high)
  
  a<-read.table("rna_cmb.txt",header=T,row.names=1)
  clic4<-a[which(rownames(a) == "CLIC4|25932_calculated"),]
  
  rna_pval<-matrix(nrow=dim(rna_cmb)[1],ncol=2)
  for (i in 1:dim(rna_cmb)[1]){
    test_res<-cor.test(as.numeric(clic4),as.numeric(rna_cmb[i,]))
    rna_pval[i,1]<-test_res$p.value
    rna_pval[i,2]<-test_res$estimate
  }
  rownames(rna_pval)<-rownames(rna_cmb)
  colnames(rna_pval)<-c("pval","corr1")
  
  sig<-rna_pval[which(rna_pval[,1]<0.05),]
  sig_value<-rna_cmb[which(rna_pval[,1]<0.05),]
  
  sig_label<-matrix(nrow=dim(rna_cmb)[1],ncol=1)
  sig_label[which(rna_pval[,1]>=0.05),1]<-"UnsigDEGs"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2]>0)),1]<-"Up"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2]<0)),1]<-"Down"
  rna_pval<-cbind(rna_pval,sig_label)
  colnames(rna_pval)<-c("pval","corr1","sig_label")
  
  write.table(rna_pval,"mirna_deg.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig,"mirna_sig_p_0.05.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(rna_cmb,"mirna_rna_cmb.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig_value,"mirna_sig_value.txt",quote=F,row.names=T,col.names=T,sep="\t")
  
  a<-read.table("mirna_deg.txt",header=T,row.names=1,sep="\t")
  ggplot(a,aes(x=pval,y=corr1,colour=sig_label))+geom_point()+scale_color_manual(breaks=c("Up","Down","UnsigDEGs"),values=c("#00A14B","black","#ED1C24"))+labs(x="Correlation test (P-value)",y="Correlation coeffients")+theme(panel.grid.major = element_blank())+theme_bw()
  
  a<-read.table("mirna_sig_value.txt",header=T,row.names=1)
  g<-heatmap.2(as.matrix(a-rowMeans(a)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",margins=c(2,10),labCol=F,cexRow=0.8)
  
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"mirna_heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"mirna_heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
}
process_miRNA()

process_mRNA_miRNA <- function(){
  rna_cmb<-read.table("rna_cmb.txt",header=T,row.names=1)
  mirna_cmb<-read.table("mirna_rna_cmb.txt",header=T,row.names=1)
  relations<-read.table("all_target_only_del_p.txt",header=F,sep="\t")
  
  corr1<-matrix(nrow=dim(relations)[1],ncol=2)
  for (i in 1:dim(relations)[1]){
    if(length(grep(pattern=paste("^",relations[i,1],"\\","|",sep=""),rownames(rna_cmb))) != 0){
      mrna<-rna_cmb[grep(pattern=paste("^",relations[i,1],"\\","|",sep=""),rownames(rna_cmb)),]
      if(dim(mrna)[1]>1){
        mrna <- colMeans(mrna)
      }
      mirna<-mirna_cmb[which(rownames(mirna_cmb) == relations[i,2]),]
      test_res<-cor.test(as.numeric(mrna),as.numeric(mirna))
      corr1[i,1]<-test_res$p.value
      corr1[i,2]<-test_res$estimate
    }
  }
  relations_new <- cbind(relations,corr1)
  colnames(relations_new)<-c("gene","mirna","p_val","corr1")
  
  sig1<-relations_new[which(relations_new[,3]<0.05),]
  sig2<-sig1[which(sig1[,4] > 0.4),]
  sig3<-sig1[which(sig1[,4] < -0.4),]
  sig<-rbind(sig2,sig3)
  sig_label<-matrix(nrow=dim(sig)[1],ncol=1)
  sig_label[which(sig[,4]>0),1]<-"Positive"
  sig_label[which(sig[,4]<0),1]<-"Negtive"
  sig_new<-cbind(sig,sig_label)
  colnames(sig_new)[5]<-c("sig_label")
  
  write.table(sig_new,"mirna_tar_sig_p_0.05_all_cor_0.4.txt",quote=F,row.names=F,col.names=T,sep="\t")
}
process_mRNA_miRNA()

process_meth <- function(){
  data<-read.csv("73CN-AML-cn_meth-TCGA.csv",row.names=1,header=T,sep=",")
  label<-read.table("clic4_label.txt",header=T)
  
  low<-label[which(label[,2]==1),1]
  high<-label[which(label[,2]==2),1]
  data.low<-data[,low]
  data.high<-data[,high]
  dnaM_cmb<-cbind(data.low,data.high)
  dnaM_cmb<-dnaM_cmb[-which(is.na(rowMeans(dnaM_cmb))),]
  
  rna_pval<-matrix(nrow=dim(dnaM_cmb)[1],ncol=2)
  for (i in 1:dim(dnaM_cmb)[1]){
    test_res<-t.test(dnaM_cmb[i,1:37],dnaM_cmb[i,38:73])
    log2fc<-log2(rowMeans(dnaM_cmb[i,38:73])/rowMeans(dnaM_cmb[i,1:37]))#high/low
    rna_pval[i,1]<-test_res$p.value
    rna_pval[i,2]<-log2fc
  }
  rownames(rna_pval)<-rownames(dnaM_cmb)
  colnames(rna_pval)<-c("pval","log2fc(high/low)")
  
  
  sig_label<-matrix(nrow=dim(dnaM_cmb)[1],ncol=1)
  sig_label[which(rna_pval[,1]>=0.05),1]<-"UnsigDEGs"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2] > 1)),1]<-"Up"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2] < -1)),1]<-"Down"
  sig_label[intersect(intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2] >= -1)),which(rna_pval[,2] <= 1)),1]<-"Marginal sigDEGs"
  
  rna_pval_new<-cbind(rna_pval,sig_label)
  colnames(rna_pval_new)<-c("pval","log2fc(high/low)","sig_label")
  
  sig<-rbind(rna_pval_new[which(rna_pval_new[,3] == "Up"),],rna_pval_new[which(rna_pval_new[,3] == "Down"),])
  sig_value<-rbind(dnaM_cmb[which(rna_pval_new[,3] == "Up"),],dnaM_cmb[which(rna_pval_new[,3] == "Down"),])
  
  write.table(rna_pval_new,"meth_deg.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig,"meth_sig_p_0.05_fc_1.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(dnaM_cmb,"meth_dnaM_cmb.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig_value,"meth_sig_value.txt",quote=F,row.names=T,col.names=T,sep="\t")
  
  
  ###################scatter
  a<-read.table("meth_deg.txt",header=T,row.names=1,sep="\t")
  ggplot(a,aes(x=log2fc.high.low.,y=-log10(a$pval),colour=sig_label))+geom_point()+scale_color_manual(breaks=c("Up","Down","Marginal sigDEGs","UnsigDEGs"),values=c("#00A14B","#2B3990","black","#ED1C24"))+labs(x="log2(Fold Change)",y="-log10(p-value)")+theme(panel.grid.major = element_blank())+theme_bw()+ggsave("meth_deg.png",width = 11, height = 8)
  
  ######################annotation
  #########clic4
  a<-read.table("meth_sig_p_0.05_fc_1.txt",header=T,sep="\t")
  anno<-read.csv("meth_sig_p_0.05_fc_1_anno.txt",header=T,sep=",")
  anno_select<-anno[,c(1,12,13,22,24,25,26,28,32)]
  merge_sig<-merge(a[,c(1,4)],anno_select,by.y="IlmnID",by.x="X",sort=F)
  merge_sig_up<-merge_sig[1:675,]
  merge_sig_down<-merge_sig[676:940,]
  table(merge_sig_up$Relation_to_UCSC_CpG_Island)
  table(merge_sig_down$Relation_to_UCSC_CpG_Island)
  
  island <-matrix(c(267, 173, 408, 92),nrow = 2)#267+408=675; 173+92=265
  nshelf <-matrix(c(16, 0, 659, 265),nrow = 2)
  nshore <-matrix(c(71, 42, 604, 223),nrow = 2)
  sshelf <-matrix(c(20, 1, 655, 264),nrow = 2)
  sshore <-matrix(c(79, 29, 596, 236),nrow = 2)
  other <-matrix(c(222, 20, 453, 245),nrow = 2)
  
  fisher.test(island)
  fisher.test(nshelf)
  fisher.test(nshore)
  fisher.test(sshelf)
  fisher.test(sshore)
  fisher.test(other)
  
  tss200 <-matrix(c(162, 89, 984, 449),nrow = 2)
  tss1500 <-matrix(c(172, 115, 974, 423),nrow = 2)
  body1 <-matrix(c(374, 138, 772, 400),nrow = 2)
  utr5 <-matrix(c(150, 107, 996, 431),nrow = 2)
  utr3 <-matrix(c(40, 2, 1106, 536),nrow = 2)
  exon <-matrix(c(78, 58, 1068, 480),nrow = 2)
  
  fisher.test(tss200)
  fisher.test(tss1500)
  fisher.test(body1)
  fisher.test(utr5)
  fisher.test(utr3)
  fisher.test(exon)
  
  write.table(merge_sig,"meth_sig_p_0.05_fc_1_anno_select.txt",quote=F,row.names=F,col.names=T,sep="\t")
  
  ###################heatmap
  a<-read.table("meth_sig_p_0.05_fc_1.txt",header=T,sep="\t")
  anno<-read.csv("meth_sig_p_0.05_fc_1_anno.txt",header=T,sep=",")
  anno_select<-anno[,c(1,12,13,22,24,25,26,28,32)]
  merge_sig<-merge(a[,c(1,4)],anno_select,by.y="IlmnID",by.x="X",sort=F)
  sig_value<-read.table("meth_sig_value.txt",header=T,row.names=1)
  
  merge_sig_island_value<-sig_value[which(merge_sig[,8] == "Island"),]
  merge_sig_promo_value<-sig_value[which(merge_sig[,10] == "Promoter_Associated"),]
  
  g<-heatmap.2(as.matrix(merge_sig_island_value-rowMeans(merge_sig_island_value)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",labCol=F,labRow=F)
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"meth_island_heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"meth_island_heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
  
  g<-heatmap.2(as.matrix(merge_sig_promo_value-rowMeans(merge_sig_promo_value)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",labCol=F,labRow=F)
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"meth_promo_heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"meth_promo_heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
  
}

process_meth()


#=====================================================================================
# WGCNA
#=====================================================================================

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

