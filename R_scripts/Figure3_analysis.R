#####################
# Figure 3 Analysis #
#####################
library(reshape2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(VennDiagram)
library(cowplot)
library(WGCNA)
library(flashClust)
source("/ifs/devel/michaelm/cgat/R/summarySE.R")
source("/ifs/devel/michaelm/cgat/R/clusterEigengenes.R")
source("crossCorrelation.R")

fo.gc <- read.table("/ifs/projects/proj036/jethro_bamfiles/deseq.dir/FovsGC_DESeq2_Jethro.tsv", h=T,
                    stringsAsFactors=F, sep="\t")
fo.gc$GeneID <- rownames(fo.gc)

cd40.intersect <- read.table("/ifs/projects/proj036/pipeline_nofilter_functional/FovsGC_compare.dir/CD40L-timepoint_intersection_genes.tsv",
                             h=T, stringsAsFactors=F, sep="\t")

lps.intersect <- read.table("/ifs/projects/proj036/pipeline_nofilter_functional/FovsGC_compare.dir/LPS-timepoint_intersection_genes.tsv",
                            h=T, stringsAsFactors=F, sep="\t")

igm.intersect <- read.table("/ifs/projects/proj036/pipeline_nofilter_functional/FovsGC_compare.dir/IgM-timepoint_intersection_genes.tsv",
                            h=T, stringsAsFactors=F, sep="\t")

# select the differentially expressed genes and lncRNAs with p <= 0.01
fo.gc.de.genes <- fo.gc$GeneID[fo.gc$padj <= 0.01]
fo.gc.de.genes <- fo.gc.de.genes[grepl(fo.gc.de.genes, pattern="ENSM")]

fo.gc.de.lncs <- fo.gc$GeneID[fo.gc$padj <= 0.01]
fo.gc.de.lncs <- fo.gc.de.lncs[grepl(fo.gc.de.lncs, pattern="LNC")]

length(fo.gc.de.genes)
length(fo.gc.de.lncs)

cd40.fo_gc.intersect.genes <- intersect(cd40.intersect$gene_id, fo.gc.de.genes)
igm.fo_gc.intersect.genes <- intersect(igm.intersect$gene_id, fo.gc.de.genes)
lps.fo_gc.intersect.genes <- intersect(lps.intersect$gene_id, fo.gc.de.genes)

# this get all of the genes differentially expressed over time and between
# GCs and follicular B cells.
fo.gc.intersect.all <- intersect(intersect(cd40.fo_gc.intersect.genes, lps.fo_gc.intersect.genes),
                                 igm.fo_gc.intersect.genes)

# get the TPM estimates of expression for all genes and lncRNAs to do the temporal
# correlation analysis
# the CD40L + IL-4 table contains the 0-time points
# these need to be inserted in to the beginning of the LPS and anti-IgM tables

cd40l.tpm <- read.table("/ifs/projects/proj036/pipeline_tpm_refcoding/tpm.dir/CD40L.tpm",
                        sep="\t", h=T, stringsAsFactors=F)

igm.tpm <- read.table("/ifs/projects/proj036/pipeline_tpm_refcoding/tpm.dir/IgM.tpm",
                        sep="\t", h=T, stringsAsFactors=F)
igm.tpm <- as.data.frame(append(igm.tpm, cd40l.tpm[, c(2:4)], after=1))

lps.tpm <- read.table("/ifs/projects/proj036/pipeline_tpm_refcoding/tpm.dir/LPS.tpm",
                        sep="\t", h=T, stringsAsFactors=F)
lps.tpm <- as.data.frame(append(lps.tpm, cd40l.tpm[, c(2:4)], after=1))

# We only want to cluster over the first 48 hours, so the first 21
# columns should contain 3 reps for each time point and 0, 1, 3, 6, 12, 24, 48
# aggregate over replicates (mean)

cd40l.fogc.genes <- cd40l.tpm[cd40l.tpm$Name %in% fo.gc.intersect.all, ]
rownames(cd40l.fogc.genes) <- cd40l.fogc.genes$Name
cd40l.fogc.genes <- cd40l.fogc.genes[, c(2:22)]
reps <- unique(unlist(lapply(strsplit(colnames(cd40l.fogc.genes),
                                      fixed=T, split="."), FUN=function(x) paste0(x[3]))))
times <-unique(unlist(lapply(strsplit(colnames(cd40l.fogc.genes),
                                      fixed=T, split="."), FUN=function(x) paste0(x[2]))))
cd40l.trans <- data.frame(t(cd40l.fogc.genes))
cd40l.trans$times <- times
cd40l.trans$reps <- reps
cd40l.melt <- melt(cd40l.trans, id.vars=c("times", "reps"))
cd40l.sum <- summarySE(cd40l.melt, measurevar="value",
                       groupvars=c("times", "variable"))
cd40l.agg <- data.frame(cd40l.sum$times, cd40l.sum$variable,
                        cd40l.sum$value)
colnames(cd40l.agg) <- c("times", "gene", "value")
cd40l.cast <- dcast(cd40l.agg, gene ~ times, value.var="value")
rownames(cd40l.cast) <- cd40l.cast$gene
cd40l.cast <- cd40l.cast[, -1]

# aggregate and munge IgM data
igm.fogc.genes <- igm.tpm[igm.tpm$Name %in% fo.gc.intersect.all, ]
rownames(igm.fogc.genes) <- igm.fogc.genes$Name
igm.fogc.genes <- igm.fogc.genes[, c(2:22)]

igm.trans <- data.frame(t(cd40l.fogc.genes))
igm.trans$times <- times
igm.trans$reps <- reps
igm.melt <- melt(igm.trans, id.vars=c("times", "reps"))
igm.sum <- summarySE(igm.melt, measurevar="value",
                       groupvars=c("times", "variable"))
igm.agg <- data.frame(igm.sum$times, igm.sum$variable,
                        igm.sum$value)
colnames(igm.agg) <- c("times", "gene", "value")
igm.cast <- dcast(igm.agg, gene ~ times, value.var="value")
rownames(igm.cast) <- igm.cast$gene
igm.cast <- igm.cast[, -1]

# aggregate and munge LPS data
lps.fogc.genes <- lps.tpm[lps.tpm$Name %in% fo.gc.intersect.all, ]
rownames(lps.fogc.genes) <- lps.fogc.genes$Name
lps.fogc.genes <- lps.fogc.genes[, c(2:22)]

lps.trans <- data.frame(t(cd40l.fogc.genes))
lps.trans$times <- times
lps.trans$reps <- reps
lps.melt <- melt(lps.trans, id.vars=c("times", "reps"))
lps.sum <- summarySE(lps.melt, measurevar="value",
                     groupvars=c("times", "variable"))
lps.agg <- data.frame(lps.sum$times, lps.sum$variable,
                      lps.sum$value)
colnames(lps.agg) <- c("times", "gene", "value")
lps.cast <- dcast(lps.agg, gene ~ times, value.var="value")
rownames(lps.cast) <- lps.cast$gene
lps.cast <- lps.cast[, -1]

lps.fogc.genes <- lps.tpm[lps.tpm$Name %in% fo.gc.intersect.all, ]
rownames(lps.fogc.genes) <- lps.fogc.genes$Name
lps.fogc.genes <- lps.fogc.genes[, c(2:22)]

hmcol <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))

###################################################
# Calculate CD40L + IL-4 cross-correlation matrix #
###################################################
cd40l.idx <- as.matrix(expand.grid(seq_len(nrow(cd40l.cast)),
                                   seq_len(nrow(cd40l.cast))))

cd40l.func <- function(x){
  x_val <- cd40l.cast[x[1], ]
  y_val <- cd40l.cast[x[2], ]
  crossCorrelate(x_val, y_val, lag=0)
}

cd40l.xcor <- apply(unlist(cd40l.idx), 1, FUN=cd40l.func)
cd40l.xcor.mat <- matrix(cd40l.xcor, nrow=nrow(cd40l.cast),
                         ncol=nrow(cd40l.cast))
cd40l.xcor.mat[is.na(cd40l.xcor.mat)] <- 0

png("manuscript_plots/CD40L-core_genes-cross_correlation-heatmap.png",
    height=12, width=12, units="in", res=300)
heatmap.2(cd40l.xcor.mat, trace="none", col=hmcol, labRow=F,
          labCol=F, density.info="none")
dev.off()

################################
# IgM cross-correlation matrix #
################################
igm.idx <- as.matrix(expand.grid(seq_len(nrow(igm.cast)),
                                 seq_len(nrow(igm.cast))))

igm.func <- function(x){
  x_val <- igm.cast[x[1], ]
  y_val <- igm.cast[x[2], ]
  crossCorrelate(x_val, y_val, lag=0)
}

igm.xcor <- apply(unlist(igm.idx), 1, FUN=igm.func)
igm.xcor.mat <- matrix(igm.xcor, nrow=nrow(igm.cast),
                         ncol=nrow(igm.cast))
igm.xcor.mat[is.na(igm.xcor.mat)] <- 0

png("manuscript_plots/IgM-core_genes-cross_correlation-heatmap.png",
    height=12, width=12, units="in", res=300)
heatmap.2(igm.xcor.mat, trace="none", col=hmcol, labRow=F,
          labCol=F, density.info="none")
dev.off()

################################
# LPS cross-correlation matrix #
################################
lps.idx <- as.matrix(expand.grid(seq_len(nrow(lps.cast)),
                                 seq_len(nrow(lps.cast))))

lps.func <- function(x){
  x_val <- lps.cast[x[1], ]
  y_val <- lps.cast[x[2], ]
  crossCorrelate(x_val, y_val, lag=0)
}

lps.xcor <- apply(unlist(lps.idx), 1, FUN=lps.func)
lps.xcor.mat <- matrix(lps.xcor, nrow=nrow(lps.cast),
                       ncol=nrow(lps.cast))
lps.xcor.mat[is.na(lps.xcor.mat)] <- 0

png("manuscript_plots/LPS-core_genes-cross_correlation-heatmap.png",
    height=12, width=12, units="in", res=300)
heatmap.2(lps.xcor.mat, trace="none", col=hmcol, labRow=F,
          labCol=F, density.info="none")
dev.off()

############################################
# Cluster the correlation matrices as      #
# 1 - corr using hierarchical clustering   #
# then cut the tree with dynamic tree      #
# cutting and assign co-expression modules #
############################################

cd40l.dist <- as.dist(abs(cd40l.xcor.mat - 1))
igm.dist <- as.dist(abs(igm.xcor.mat - 1))
lps.dist <- as.dist(abs(lps.xcor.mat - 1))

cd40l.clust <- flashClust(cd40l.dist, method="average")
igm.clust <- flashClust(igm.dist, method="average")
lps.clust <- flashClust(lps.dist, method="average")

cd40l.cut <- cutreeDynamic(dendro=cd40l.clust, method="tree",
                           minClusterSize=30, deepSplit=T)
cd40l.colors <- labels2colors(cd40l.cut)

igm.cut <- cutreeDynamic(dendro=igm.clust, method="tree",
                           minClusterSize=30, deepSplit=T)
igm.colors <- labels2colors(igm.cut)

lps.cut <- cutreeDynamic(dendro=lps.clust, method="tree",
                         minClusterSize=30, deepSplit=T)
lps.colors <- labels2colors(lps.cut)

png("manuscript_plots/CD40L-core_genes-TreeCut_clustering-dendrogram.png",
    height=9, width=9, units="in", res=300)
plotDendroAndColors(cd40l.clust, colors=cd40l.colors,
                    groupLabels="Dynamic tree cut",
                    dendroLabels=F, addGuide=T, guideHang=0.05,
                    hang=0.03, main="CD40L + IL-4\nCore Genes")
dev.off()

png("manuscript_plots/IgM-core_genes-TreeCut_clustering-dendrogram.png",
    height=9, width=9, units="in", res=300)
plotDendroAndColors(igm.clust, colors=igm.colors,
                    groupLabels="Dynamic tree cut",
                    dendroLabels=F, addGuide=T, guideHang=0.05,
                    hang=0.03, main="anti-IgM\nCore Genes")
dev.off()

png("manuscript_plots/LPS-core_genes-TreeCut_clustering-dendrogram.png",
    height=9, width=9, units="in", res=300)
plotDendroAndColors(lps.clust, colors=lps.colors,
                    groupLabels="Dynamic tree cut",
                    dendroLabels=F, addGuide=T, guideHang=0.05,
                    hang=0.03, main="LPS\nCore Genes")
dev.off()

########################################
# Compute eigenngenes for each cluster #
# CD40L + IL-4                         #
########################################
cd40l.clust.match <- data.frame(cbind(rownames(cd40l.cast), cd40l.colors))
colnames(cd40l.clust.match) <- c("genes", "cluster")

cd40l.eigenclust <- clusterPCA(cluster_frame=cd40l.clust.match,
                               expression_frame=cd40l.cast, n=as.numeric(times))

cd40l.eigenframe <- eigenExpress(cd40l.eigenclust, n=as.numeric(times))
cd40l.eigen.melt <- melt(cd40l.eigenframe, id.vars="cluster")
cd40l.eigen.melt$variable = as.numeric(as.character(cd40l.eigen.melt$variable))
cd40l.eigen.melt$value = as.numeric(cd40l.eigen.melt$value)

p_cd40l.eigen <- ggplot(cd40l.eigen.melt, aes(x=variable, y=value, colour=cluster)) +
  geom_line(size=2) + geom_point(size=2) + 
  theme_bw() + scale_x_continuous(limits=c(0, 50), breaks=cd40l.eigen.melt$variable) + 
  scale_y_continuous(limits=c(-20, 20)) +
  theme(axis.text=element_text(size=14, colour="black")) + 
  theme(axis.title=element_text(size=14, colour="black")) + 
  theme(legend.text=element_text(size=14, colour="black")) +
  theme(legend.title=element_text(size=14, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), legend.key=element_blank()) +
  labs(x="Time (hours)", y="Eigengene Expression") + 
  guides(colour=guide_legend(title="Coexpression\nCluster")) +
  scale_colour_manual(values=levels(cd40l.eigen.melt$cluster))

ggsave(p_cd40l.eigen, filename="manuscript_plots/CD40L-Core_clusters-eigengene_expression.png",
       height=8, width=12, dpi=300)

########################################
# Compute eigenngenes for each cluster #
# anti-IgM                             #
########################################
igm.clust.match <- data.frame(cbind(rownames(igm.cast), igm.colors))
colnames(igm.clust.match) <- c("genes", "cluster")

igm.eigenclust <- clusterPCA(cluster_frame=igm.clust.match,
                               expression_frame=igm.cast, n=as.numeric(times))

igm.eigenframe <- eigenExpress(igm.eigenclust, n=as.numeric(times))
igm.eigen.melt <- melt(igm.eigenframe, id.vars="cluster")
igm.eigen.melt$variable = as.numeric(as.character(igm.eigen.melt$variable))
igm.eigen.melt$value = as.numeric(igm.eigen.melt$value)

p_igm.eigen <- ggplot(igm.eigen.melt, aes(x=variable, y=value, colour=cluster)) +
  geom_line(size=2) + geom_point(size=2) + 
  theme_bw() + scale_x_continuous(limits=c(0, 50), breaks=igm.eigen.melt$variable) + 
  scale_y_continuous(limits=c(-20, 20)) +
  theme(axis.text=element_text(size=14, colour="black")) + 
  theme(axis.title=element_text(size=14, colour="black")) + 
  theme(legend.text=element_text(size=14, colour="black")) +
  theme(legend.title=element_text(size=14, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), legend.key=element_blank()) +
  labs(x="Time (hours)", y="Eigengene Expression") + 
  guides(colour=guide_legend(title="Coexpression\nCluster")) +
  scale_colour_manual(values=levels(igm.eigen.melt$cluster))

ggsave(p_igm.eigen, filename="manuscript_plots/IgM-Core_clusters-eigengene_expression.png",
       height=8, width=12, dpi=300)

########################################
# Compute eigenngenes for each cluster #
# LPS                                  #
########################################
lps.clust.match <- data.frame(cbind(rownames(lps.cast), lps.colors))
colnames(lps.clust.match) <- c("genes", "cluster")

lps.eigenclust <- clusterPCA(cluster_frame=lps.clust.match,
                             expression_frame=lps.cast, n=as.numeric(times))

lps.eigenframe <- eigenExpress(lps.eigenclust, n=as.numeric(times))
lps.eigen.melt <- melt(lps.eigenframe, id.vars="cluster")
lps.eigen.melt$variable = as.numeric(as.character(lps.eigen.melt$variable))
lps.eigen.melt$value = as.numeric(lps.eigen.melt$value)

p_lps.eigen <- ggplot(lps.eigen.melt, aes(x=variable, y=value, colour=cluster)) +
  geom_line(size=2) + geom_point(size=2) + 
  theme_bw() + scale_x_continuous(limits=c(0, 50), breaks=lps.eigen.melt$variable) + 
  scale_y_continuous(limits=c(-20, 20)) +
  theme(axis.text=element_text(size=14, colour="black")) + 
  theme(axis.title=element_text(size=14, colour="black")) + 
  theme(legend.text=element_text(size=14, colour="black")) +
  theme(legend.title=element_text(size=14, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), legend.key=element_blank()) +
  labs(x="Time (hours)", y="Eigengene Expression") + 
  guides(colour=guide_legend(title="Coexpression\nCluster")) +
  scale_colour_manual(values=levels(lps.eigen.melt$cluster))
ggsave(p_lps.eigen, filename="manuscript_plots/LPS-Core_clusters-eigengene_expression.png",
       height=8, width=12, dpi=300)

#############################################
# Cross-correlation of all lncRNAs with the #
# Core genes in each in vitro activation    #
# condition CD40L + IL-4                    #
#############################################

cd40l.lncs <- read.table("/ifs/projects/proj036/pipeline_sailfish/tpm.dir/CD40L.tpm",
                         sep="\t", h=T)
rownames(cd40l.lncs) <- cd40l.lncs$Name
cd40l.lncs <- cd40l.lncs[, c(2:22)]
reps <- unique(unlist(lapply(strsplit(colnames(cd40l.lncs),
                                      fixed=T, split="."), FUN=function(x) paste0(x[3]))))
times <-unique(unlist(lapply(strsplit(colnames(cd40l.lncs),
                                      fixed=T, split="."), FUN=function(x) paste0(x[2]))))
cd40l.trans.lncs <- data.frame(t(cd40l.lncs))
cd40l.trans.lncs$times <- times
cd40l.trans.lncs$reps <- reps
cd40l.melt.lncs <- melt(cd40l.trans.lncs, id.vars=c("times", "reps"))
cd40l.sum.lncs <- summarySE(cd40l.melt.lncs, measurevar="value",
                       groupvars=c("times", "variable"))
cd40l.agg.lncs <- data.frame(cd40l.sum.lncs$times, cd40l.sum.lncs$variable,
                        cd40l.sum.lncs$value)
colnames(cd40l.agg.lncs) <- c("times", "gene", "value")
cd40l.cast.lncs <- dcast(cd40l.agg.lncs, gene ~ times, value.var="value")
rownames(cd40l.cast.lncs) <- cd40l.cast.lncs$gene
cd40l.cast.lncs <- cd40l.cast.lncs[, -1]

# remove all zero-expression lncRNAs
cd40l.cast.lncs <- cd40l.cast.lncs[rowSums(cd40l.cast.lncs < 1) < 5,]

cd40l.idx.lncs <- as.data.frame(expand.grid(seq_len(nrow(cd40l.cast.lncs)),
                                            seq_len(nrow(cd40l.cast))))
                      
cd40l.func.lncs <- function(x){
  x_val <- cd40l.cast.lncs[x[1], ]
  y_val <- cd40l.cast[x[2], ]
  crossCorrelate(x_val, y_val, lag=0)
}

cd40l.lncs.xcor <- apply(cd40l.idx.lncs, 1, FUN=cd40l.func.lncs)
cd40l.lncs.xcor.mat <- matrix(cd40l.lncs.xcor, nrow=nrow(cd40l.cast.lncs),
                         ncol=nrow(cd40l.cast.lncs))
cd40l.lncs.xcor.mat[is.na(cd40l.lncs.xcor.mat)] <- 0

png("manuscript_plots/CD40L-lncRNA_vs_Coregenes-cross_correlation-heatmap.png",
    height=12, width=12, units="in", res=300)
heatmap.2(cd40l.lncs.xcor.mat, trace="none", col=hmcol, labRow=F,
          labCol=F, density.info="none")
dev.off()

#############################################
# Cross-correlation of all lncRNAs with the #
# Core genes in each in vitro activation    #
# condition anti-IgM                        #
#############################################

igm.lncs <- read.table("/ifs/projects/proj036/pipeline_sailfish/tpm.dir/IgM.tpm",
                         sep="\t", h=T)
igm.lncs <- as.data.frame(append(igm.lncs, cd40l.lncs[, c(1:3)], after=1))

rownames(igm.lncs) <- igm.lncs$Name
igm.lncs <- igm.lncs[, c(2:22)]
reps <- unique(unlist(lapply(strsplit(colnames(igm.lncs),
                                      fixed=T, split="."), FUN=function(x) paste0(x[3]))))
times <-unique(unlist(lapply(strsplit(colnames(igm.lncs),
                                      fixed=T, split="."), FUN=function(x) paste0(x[2]))))
igm.trans.lncs <- data.frame(t(igm.lncs))
igm.trans.lncs$times <- times
igm.trans.lncs$reps <- reps
igm.melt.lncs <- melt(igm.trans.lncs, id.vars=c("times", "reps"))
igm.sum.lncs <- summarySE(igm.melt.lncs, measurevar="value",
                            groupvars=c("times", "variable"))
igm.agg.lncs <- data.frame(igm.sum.lncs$times, igm.sum.lncs$variable,
                             igm.sum.lncs$value)
colnames(igm.agg.lncs) <- c("times", "gene", "value")
igm.cast.lncs <- dcast(igm.agg.lncs, gene ~ times, value.var="value")
rownames(igm.cast.lncs) <- igm.cast.lncs$gene
igm.cast.lncs <- igm.cast.lncs[, -1]

# remove all zero-expression lncRNAs
igm.cast.lncs <- igm.cast.lncs[rowSums(igm.cast.lncs < 1) < 5,]

igm.idx.lncs <- as.data.frame(expand.grid(seq_len(nrow(igm.cast.lncs)),
                                          seq_len(nrow(igm.cast))))

igm.func.lncs <- function(x){
  x_val <- igm.cast.lncs[x[1], ]
  y_val <- igm.cast[x[2], ]
  crossCorrelate(x_val, y_val, lag=0)
}

igm.lncs.xcor <- apply(igm.idx.lncs, 1, FUN=igm.func.lncs)
igm.lncs.xcor.mat <- matrix(igm.lncs.xcor, nrow=nrow(igm.cast.lncs),
                              ncol=nrow(igm.cast.lncs))
igm.lncs.xcor.mat[is.na(igm.lncs.xcor.mat)] <- 0

png("manuscript_plots/IgM-lncRNA_vs_Coregenes-cross_correlation-heatmap.png",
    height=12, width=12, units="in", res=300)
heatmap.2(igm.lncs.xcor.mat, trace="none", col=hmcol, labRow=F,
          labCol=F, density.info="none")
dev.off()

#############################################
# Cross-correlation of all lncRNAs with the #
# Core genes in each in vitro activation    #
# condition LPS                             #
#############################################

lps.lncs <- read.table("/ifs/projects/proj036/pipeline_sailfish/tpm.dir/LPS.tpm",
                       sep="\t", h=T)
lps.lncs <- as.data.frame(append(lps.lncs, cd40l.lncs[, c(1:3)], after=1))

rownames(lps.lncs) <- lps.lncs$Name
lps.lncs <- lps.lncs[, c(2:22)]
reps <- unique(unlist(lapply(strsplit(colnames(lps.lncs),
                                      fixed=T, split="."), FUN=function(x) paste0(x[3]))))
times <-unique(unlist(lapply(strsplit(colnames(lps.lncs),
                                      fixed=T, split="."), FUN=function(x) paste0(x[2]))))
lps.trans.lncs <- data.frame(t(lps.lncs))
lps.trans.lncs$times <- times
lps.trans.lncs$reps <- reps
lps.melt.lncs <- melt(lps.trans.lncs, id.vars=c("times", "reps"))
lps.sum.lncs <- summarySE(lps.melt.lncs, measurevar="value",
                          groupvars=c("times", "variable"))
lps.agg.lncs <- data.frame(lps.sum.lncs$times, lps.sum.lncs$variable,
                           lps.sum.lncs$value)
colnames(lps.agg.lncs) <- c("times", "gene", "value")
lps.cast.lncs <- dcast(lps.agg.lncs, gene ~ times, value.var="value")
rownames(lps.cast.lncs) <- lps.cast.lncs$gene
lps.cast.lncs <- lps.cast.lncs[, -1]

# remove all zero-expression lncRNAs
lps.cast.lncs <- lps.cast.lncs[rowSums(lps.cast.lncs < 1) < 5,]

lps.idx.lncs <- as.data.frame(expand.grid(seq_len(nrow(lps.cast.lncs)),
                                          seq_len(nrow(lps.cast))))

lps.func.lncs <- function(x){
  x_val <- lps.cast.lncs[x[1], ]
  y_val <- lps.cast[x[2], ]
  crossCorrelate(x_val, y_val, lag=0)
}

lps.lncs.xcor <- apply(lps.idx.lncs, 1, FUN=lps.func.lncs)
lps.lncs.xcor.mat <- matrix(lps.lncs.xcor, nrow=nrow(lps.cast.lncs),
                            ncol=nrow(lps.cast.lncs))
lps.lncs.xcor.mat[is.na(lps.lncs.xcor.mat)] <- 0

png("manuscript_plots/LPS-lncRNA_vs_Coregenes-cross_correlation-heatmap.png",
    height=12, width=12, units="in", res=300)
heatmap.2(lps.lncs.xcor.mat, trace="none", col=hmcol, labRow=F,
          labCol=F, density.info="none")
dev.off()
