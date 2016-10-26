################################################
# LNCGme01880 differential expression analysis #
################################################
library(DESeq2)
library(ggplot2)
library(gplots)
library(GGally)
library(RColorBrewer)
library(goseq)
library(cowplot)
library(org.Mm.eg.db)
library(pathview)
library(KEGG.db)
library(cowplot)
library(reactome.db)
library(biomaRt)
library(reshape2)

follic_counts <- read.table("/ifs/projects/proj036/pipeline_counts_LNCGme01880_KO/feature_counts.dir/Follicular-Ensembl72_all_plus_lncRNA-feature_counts.tsv.gz",
                            sep="\t", stringsAsFactors=F, h=T, row.names=1)
pre_counts <- read.table("/ifs/projects/proj036/pipeline_counts_LNCGme01880_KO/feature_counts.dir/preBcell-Ensembl72_all_plus_lncRNA-feature_counts.tsv.gz",
                         sep="\t", stringsAsFactors=F, h=T, row.names=1)

# need to set the condition factor levels to test WT as reference
follic.design <- data.frame("Replicate"=as.factor(rep(c(1, 2, 3), 2)),
                            "Condition"=c(rep("KO", 3), rep("WT", 3)),
                            "Cell"=rep("Follicular", 6))
rownames(follic.design) <- colnames(follic_counts)

pre.design <- data.frame("Replicate"=as.factor(rep(c(1, 2, 3),2)),
                         "Condition"=c(rep("KO", 3), rep("WT", 3)),
                         "Cell"=rep("pre_Bcell", 6))
rownames(pre.design) <- colnames(pre_counts)

follic_counts$GeneID <- rownames(follic_counts)
pre_counts$GeneID <- rownames(pre_counts)

merge_counts = merge(follic_counts, pre_counts, by="GeneID")
rownames(merge_counts) <- merge_counts$GeneID
merge_counts <- merge_counts[,c(2:13)]
head(merge_counts)

merge.design <- data.frame(do.call(rbind, list(follic.design, pre.design)))
rownames(merge.design) <- colnames(merge_counts)

# create DESeq data objects and do differential analysis
dds <- DESeqDataSetFromMatrix(countData = merge_counts,
                              colData = merge.design,
                              design= ~ Replicate + Cell + Condition)

# filter out non-expressed genes and lncRNAs
zeros = rowSums(counts(dds)) <= 5
dds <- dds[!zeros, ]

dds <- DESeq(dds)
res <- results(dds, alpha = 0.01, contrast=c("Condition", "KO", "WT"))
plotMA(res)

# how many significant hits at p<0.01?
sum(res$padj < 0.01, na.rm=TRUE)

# use log transformation with ridge value for visualisation
rlog.dds <- rlog(dds, blind=TRUE)

# pca and plot
rlog.pca <- plotPCA(rlog.dds, intgroup=c("Condition", "Cell", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(rlog.pca, "percentVar"))

all_pca <- ggplot(rlog.pca, aes(x=PC1, y=PC2, size=Cell, colour=Condition, shape=Replicate)) +
  geom_point() + 
  xlab(paste("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste("PC2: ", percentVar[2], "% variance")) + 
  theme_bw()

ggsave(all_pca, filename="/ifs/projects/proj036/R_sessions/LNCGme01880_KO-PCA.png",
       height=6, width=6)

# most of the variation is between cell types, need to do the analysis separately
follic_counts <- follic_counts[,c(1:6)]
dds.follic <- DESeqDataSetFromMatrix(countData = follic_counts,
                                     colData = follic.design,
                                     design = ~ Replicate + Condition)

zeros.follic = rowSums(counts(dds.follic)) <= 5
dds.follic <- dds.follic[!zeros.follic, ]

dds.follic <- DESeq(dds.follic)
res.follic <- results(dds.follic, alpha=0.01, contrast=c("Condition", "KO", "WT"))
plotMA(res.follic)
sum(res.follic$padj < 0.01, na.rm=TRUE)
res.follic[order(res.follic$padj, decreasing = F),]

# There are four differentially expressed genes in the follicular B cells,
# 3 of which are Heatshock proteins, the fourth is Tcirg1

# use log transformation with ridge value for visualisation
rlog.follic <- rlog(dds.follic, blind=TRUE)

# pca and plot
pca.follic <- plotPCA(rlog.follic, intgroup=c("Condition", "Replicate"), returnData=TRUE)
percentVar.follic <- round(100 * attr(pca.follic, "percentVar"))

pca_fo <- ggplot(pca.follic, aes(x=PC1, y=PC2, colour=Condition, size=Replicate)) +
  geom_point() + 
  xlab(paste("PC1: ", percentVar.follic[1], "% variance")) + 
  ylab(paste("PC2: ", percentVar.follic[2], "% variance")) + 
  theme_bw()

pca_fo
ggsave(pca_fo, filename="/ifs/projects/proj036/R_sessions/LNCGme01880-Fo_PCA.png", height=6,
       width=6)

# pre B cells
pre_counts <- pre_counts[, c(1:6)]
dds.pre <- DESeqDataSetFromMatrix(countData = pre_counts,
                                  colData = pre.design,
                                  design = ~ Replicate + Condition)

zeros.pre = rowSums(counts(dds.pre)) <= 5
dds.pre <- dds.pre[!zeros.pre, ]

dds.pre <- DESeq(dds.pre)
res.pre <- results(dds.pre, alpha=0.01, contrast=c("Condition", "KO", "WT"))

plotMA(res.pre)
sum(res.pre$padj < 0.05, na.rm=TRUE)
reordered <- data.frame(res.pre[order(res.pre$padj, decreasing = F),])
print(head(reordered[!is.na(reordered$padj<= 0.05),], n= 25))
# There are 68 differentially expressed genes in the pre-B cells
# with p<= 0.01, 179 with p<= 0.05
write.table(res.pre, file="/ifs/projects/proj036/R_sessions/preBcell-DE_results.tsv",
            sep="\t", quote=F)


# use log transformation with ridge value for visualisation
rlog.pre <- rlog(dds.pre, blind=TRUE)

pre.up <- res.pre[res.pre$log2FoldChange > 0, ]
pre.up$padj[is.na(pre.up$padj)] <- 1.0
pre.up <- pre.up[pre.up$padj <= 0.05, ]

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pre.up.cor <- cor(t(assay(rlog.pre)[rownames(pre.up),]))

png("/ifs/projects/proj036/R_sessions/LNCGme01880_KO-upreg-correlation.png",
    height=8, width=8, units="in", res=150)
par(cex.main=0.8)
heatmap.2(pre.up.cor, trace="none", col=hmcol, margins=c(10, 10),
          labRow=F, labCol=F,
          main="LNCGme01880 KO Up-regulated genes\nPearson Correlation")
dev.off()

pre.up$GeneID <- rownames(pre.up)
write.table(pre.up, file="/ifs/projects/proj036/R_sessions/preBcell-LNCGme01880_KO-upreg.tsv",
            quote=F, row.names=F, sep="\t")

pre.down <- res.pre[res.pre$log2FoldChange < 0, ]
pre.down$padj[is.na(pre.down$padj)] <- 1.0
pre.down <- pre.down[pre.down$padj <= 0.05, ]
pre.down.cor <- cor(t(assay(rlog.pre)[rownames(pre.down),]))

png("/ifs/projects/proj036/R_sessions/LNCGme01880_KO-downreg-correlation.png",
    height=8, width=8, units="in", res=150)
par(cex.main=0.8)
heatmap.2(pre.down.cor, trace="none", col=hmcol, margins=c(10, 10),
          labRow=F, labCol=F,
          main="LNCGme01880 KO Down-regulated genes\nPearson Correlation")
dev.off()

pre.down$GeneID <- rownames(pre.down)
write.table(pre.down, file="/ifs/projects/proj036/R_sessions/preBcell-LNCGme01880_KO-downreg.tsv",
            sep="\t", quote=F, row.names=F)

diff.genes <- rownames(res.pre[!is.na(res.pre$padj), ][res.pre[!is.na(res.pre$padj), ]$padj <= 0.05,])
rlog.df <- data.frame(assay(rlog.pre))
cor.df <- cor(t(rlog.df[diff.genes, ]))
cor.melt <- melt(cor.df)
colnames(cor.melt) <- c("lncRNA_id", "gene_id", "corr")
out.cor <- cor.melt[cor.melt$lncRNA_id == "ENSMUSG00000073155", ]
write.table(out.cor, file="/ifs/projects/proj036/pipeline_functional_LNCGme01880_KO/correlations.dir/LNCGme01880_pairs.tsv",
            sep="\t", quote=F, row.names=F)

# pca and plot
# pca.pre <- prcomp(assay(rlog.pre), center=T, scale=T)
# pcs.pre <- data.frame(pca.pre$rotation)
# pcs.pre$sample <- row.names(pcs.pre)
# pre.design$sample <- row.names(pre.design)
# pcs.pre <- merge(pcs.pre, pre.design, by="sample")

pca.pre <- plotPCA(rlog.pre, intgroup=c("Condition", "Replicate"), returnData=TRUE)
percentVar.pre <- round(100 * attr(pca.pre, "percentVar"))
# percentVar.pre <- round((pca.pre$sdev ** 2)/sum(pca.pre$sdev ** 2), digits = 4)

pre_pca <- ggplot(pca.pre, aes(x=PC1, y=PC2, colour=Condition, shape=Replicate)) +
  geom_point(size=4) + 
  xlab(paste("PC1: ", percentVar.pre[1], "% variance")) + 
  ylab(paste("PC2: ", percentVar.pre[2], "% variance")) + 
  theme_bw()
pre_pca

# ggpairs(pcs.pre, mapping=aes(colour=Condition),
#         columns=c("PC1", "PC2", "PC3", "PC4", "PC5"),
#         upper=list(continuous="blank"),
#         lower=list(continuous=wrap("points")),
#         columnLabels=c(paste("PC1: ", percentVar.pre[1], "% variance"),
#                        paste("PC2: ", percentVar.pre[2], "% variance"),
#                        paste("PC3: ", percentVar.pre[3], "% variance"),
#                        paste("PC4: ", percentVar.pre[4], "% variance"),
#                        paste("PC5: ", percentVar.pre[5], "% variance"))) +
#   theme_bw()

pre_pca
ggsave(pre_pca, filename="/ifs/projects/proj036/R_sessions/LNCGme01880-preB_PCA.png",
       width=6, height=6)

# GO and KEGG enrichment on 73 pre-B cell diff genes with padj <= 0.01
# need gene length data from sailfish!
pre.sailfish <- read.table("/ifs/projects/proj036/pipeline_sailfish_LNCGme01880_KO/transcript_info.tsv",
                           h=F, stringsAsFactors=F, sep="\t")
gene.lengths <- abs(pre.sailfish[, 5] - pre.sailfish[,4])
names(gene.lengths) <- pre.sailfish[, 1]

gene.lengths.up <- gene.lengths[names(gene.lengths) %in% rownames(pre.up)]
gene.lengths.down <- gene.lengths[names(gene.lengths) %in% rownames(pre.down)]

all.genes <- row.names(res.pre)
de.genes <- rownames(res.pre)[res.pre$padj <= 0.05]
lncs <- grepl(all.genes, pattern="LNCG")

up.genes <- rownames(pre.up)
down.genes <- rownames(pre.down)

gene.vector.up <- as.integer(all.genes %in% up.genes)
gene.vector.down <- as.integer(all.genes %in% down.genes)

names(gene.vector.up) <- rownames(res.pre)
names(gene.vector.down) <- rownames(res.pre)

# make sure all the genes match up propely
gene.vector.up <- gene.vector.up[names(gene.vector.up) %in% names(gene.lengths)]
gene.vector.down <- gene.vector.down[names(gene.vector.down) %in% names(gene.lengths)]

table(gene.vector.up)
table(gene.vector.down)

pwf.up <- nullp(gene.vector.up, "mm10", "ensGene", bias.data = gene.lengths)
GO.wall.up <- goseq(pwf.up, genome="mm10", id="ensGene", test.cats=c("GO:BP", "KEGG"),
                 method="Wallenius", use_genes_without_cat = T)

pwf.down <- nullp(gene.vector.down, "mm10", "ensGene", bias.data = gene.lengths)
GO.wall.down <- goseq(pwf.down, genome="mm10", id="ensGene", test.cats=c("GO:BP", "KEGG"),
                    method="Wallenius", use_genes_without_cat = T)


# adjust p-values for multiple testing first

GO.wall.up$padjust <- p.adjust(GO.wall.up$over_represented_pvalue, method="BH")
GO.wall.up$foldEnrich <- ((GO.wall.up$numDEInCat/length(up.genes))/(GO.wall.up$numInCat/length(all.genes[!lncs])))
print(dim(GO.wall.up[GO.wall.up$padjust <= 0.05,]) )

GO.wall.down$padjust <- p.adjust(GO.wall.down$over_represented_pvalue, method="BH")
GO.wall.down$foldEnrich <- ((GO.wall.down$numDEInCat/length(down.genes))/(GO.wall.down$numInCat/length(all.genes[!lncs])))
print(dim(GO.wall.down[GO.wall.down$padjust <= 0.05,]) )

# split into GO and KEGG
GO.results.up <- GO.wall.up[!is.na(GO.wall.up$ontology),]
GO.results.down <- GO.wall.down[!is.na(GO.wall.down$ontology),]

KEGG.results.up <- GO.wall.up[is.na(GO.wall.up$ontology),]
KEGG.results.down <- GO.wall.down[is.na(GO.wall.down$ontology),]

GO.sig.up <- GO.results.up[(GO.results.up$padjust <= 0.05), ]
go_cats.up <- data.frame(GO.sig.up$category, GO.sig.up$foldEnrich)

GO.sig.down <- GO.results.down[(GO.results.down$padjust <= 0.05), ]
go_cats.down <- data.frame(GO.sig.down$category, GO.sig.down$foldEnrich)

write.table(go_cats.up, "/ifs/projects/proj036/R_sessions/GO_cats4Revigo-preBcell-LNCGme01880_KO-down_reg.tsv", sep="\t", 
            quote=F, row.names = F, col.names = F)
write.table(go_cats.down, "/ifs/projects/proj036/R_sessions/GO_cats4Revigo-preBcell-LNCGme01880_KO-up_reg.tsv", sep="\t", 
            quote=F, row.names = F, col.names = F)


# only 4 enriched categories for down-regulated genes, don't need REVIGO for that!
#go.revigo.up <- read.table("/ifs/projects/proj036/R_sessions/REVIGO_preBcell-up_reg.csv", sep=",", h=T)
#go.keep.up <- go.revigo.up[go.revigo.up$eliminated == 0,]$term_ID
#GO.revigo.up <- GO.sig.up[GO.sig.up$category %in% go.keep.up,]

go.revigo.down <- read.table("/ifs/projects/proj036/R_sessions/REVIGO_preBcell-down_reg.csv", sep=",", h=T)
go.keep.down <- go.revigo.down[go.revigo.down$eliminated == 0,]$term_ID
GO.revigo.down <- GO.sig.down[GO.sig.down$category %in% go.keep.down,]

write.table(GO.revigo.down[order(GO.revigo.down$foldEnrich, decreasing = T),], 
            file = "/ifs/projects/proj036/R_sessions/Revigo_trimmed-Results-preBcell-down_reg.tsv", sep="\t",
            quote=F, row.names=F)

top_15.go.down <- GO.revigo.down[order(GO.revigo.down$foldEnrich, decreasing = T),][1:8,]

p1_go <- ggplot(top_15.go.down, aes(x=reorder(term, -top_15.go.down$padjust),
                                         y=log2(foldEnrich))) + 
  geom_bar(stat="identity", colour="white", fill="black") +
  theme_minimal() + coord_flip() + labs(x="GO Category description",
                                        y=expression(paste(log[2], " Fold Enrichment"), 
                                                                                  sep="")) + 
  geom_hline(yintercept=2, linetype="dashed", colour="white")
p1_go

p2_go <- ggplot(top_15.go.down, aes(x=reorder(term, -top_15.go.down$padjust), 
                                         y=-log10(padjust))) + 
  geom_bar(stat="identity", colour="white", fill="black") +
  theme_minimal() + coord_flip() + labs(x="GO Category description", y=expression(paste("-", log[10], " P-value"), 
                                                                                  sep="")) + 
  geom_hline(yintercept=-log10(0.01), linetype="dashed", colour="white")
p2_go

go_plots <- plot_grid(p1_go, p2_go, nrow = 2, align = "h2", labels=c("A", "B"))
go_plots

# You can save these plots whereever you want to.  I like to keep plots all in one directory,
# that way I know where they all are!
save_plot(filename="/ifs/projects/proj036/R_sessions/GO_plots-preBcell-down_reg.png", go_plots,
          nrow = 2, ncol=1, base_aspect_ratio = 1.5, base_width = 9.6)


#########################################
# KEGG description mapping and plotting #
#########################################

KEGG.sig <- KEGG.results.down[(KEGG.results.down$padjust <= 0.05), ]

# get the mapping of KEGG Id to pathway name
kegg2id <- as.list(KEGGPATHID2NAME)

# make a vector and attach to KEGG results
kegg_desc <- list()
for(x in 1:length(KEGG.sig$category)){
  desc <- kegg2id[[KEGG.sig$category[x]]]
  kegg_desc[[KEGG.sig$category[x]]] <- desc
}

kegg.vec <- unlist(kegg_desc)
kegg.cats <- names(kegg.vec)
kegg.df <- data.frame(cbind(kegg.cats, kegg.vec))
colnames(kegg.df) <- c("category", "description")

# merge back in with original results
kegg.merge <- merge(KEGG.sig, kegg.df, by="category")

# let's have a look at that plot
p1_kegg <- ggplot(kegg.merge, aes(x=reorder(description, -kegg.merge$padjust), y=log2(foldEnrich))) + geom_bar(stat="identity", fill="black",
                                                                                                               colour="white") + 
  coord_flip() + theme_minimal() + theme(axis.text.x=element_text(colour="black")) +
  labs(x="KEGG Pathway", y=expression(paste(log[2], " Fold Enrichment"), sep="")) + 
  geom_hline(yintercept=2, linetype="dashed", colour="darkgrey")

p1_kegg

p2_kegg <- ggplot(kegg.merge, aes(x=reorder(description, -kegg.merge$padjust), y=-log10(padjust))) + geom_bar(stat="identity", 
                                                                                                              colour="white", fill="black") + 
  coord_flip() + theme_minimal() + theme(axis.text.x=element_text(colour="black")) +
  labs(x="KEGG Pathway", y="Adjusted -log10(p-value)") + geom_hline(yintercept=-log10(0.01), linetype="dashed", 
                                                                    colour="darkgrey")
p2_kegg

kegg_plots <- plot_grid(p1_kegg, p2_kegg, labels=c("A", "B"), align="h2", nrow=2)
kegg_plots

# Save this KEGG plot where you want it
save_plot("/ifs/projects/proj036/R_sessions/KEGG_plots-preBcell-down_reg.png", kegg_plots, 
          nrow = 2, ncol=1, base_aspect_ratio = 1.5, base_width = 9.6)

# this is the results from the KEGG enrichment, again just change the `file=` to the relevant
# location
write.table(kegg.merge[order(kegg.merge$padjust, decreasing=F),],
            file = "KEGG-results-preBcell-down_reg.tsv", sep="\t", quote=F, row.names=F)

############################################

# Looking at expression patterns of up- and down-regulated genes and lncRNAs
down.mat <- data.frame(t(assay(rlog.pre)[rownames(pre.down),]))
down.mat$sample <- rownames(down.mat)
down.merge <- merge(down.mat, pre.design, by="sample")
down.melt <- melt(down.merge, id.vars=c("Replicate", "Condition", "Cell", "sample"))

ggplot(down.melt, aes(x=Replicate, y=value, colour=Condition)) + 
  geom_boxplot()
