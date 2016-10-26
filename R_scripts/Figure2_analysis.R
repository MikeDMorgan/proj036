##########################################
# Overlap between Fo vs GC DE genes      #
# and condition-specific DEs consistent  #
# over the first 48 hours                #
##########################################
library(reshape2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(VennDiagram)
library(goseq)
library(cowplot)
library(org.Mm.eg.db)
library(pathview)
library(KEGG.db)
library(cowplot)
library(reactome.db)
library(biomaRt)

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

# Venn diagrams of each condition intersect with Fo vs. GC DE genes and lncRNAs
png("manuscript_plots/CD40L-FovsGC-intersection-Venn.png", width=7.5, height=7.5,
    res=300, units="in")
draw.pairwise.venn(area1=length(cd40.intersect$gene_id), area2=length(fo.gc.de.genes),
                   cross.area=length(cd40.fo_gc.intersect.genes),
                   category=c("CD40L +\nIL-4", "Fo vs. GC"),
                   ext.text=F, lwd=c(0, 0), fill=c("#8b008b", "#ff8c00"),
                   alpha=0.5, cex=1.1, fontface=rep(2, 3),
                   fontfamily=rep("sans", 3), cat.fontface=rep(2, 2),
                   cat.fontfamily=rep("sans", 2), margin=0.06,
                   cat.just=list(c(0, 0), c(1, 1)))
dev.off()


# Venn diagrams of each condition intersect with Fo vs. GC DE genes and lncRNAs
png("manuscript_plots/IgM-FovsGC-intersection-Venn.png", width=7.5, height=7.5,
    res=300, units="in")
draw.pairwise.venn(area1=length(igm.intersect$gene_id), area2=length(fo.gc.de.genes),
                   cross.area=length(igm.fo_gc.intersect.genes),
                   category=c("anti-IgM", "Fo vs. GC"),
                   ext.text=F, lwd=c(0, 0), fill=c("#8b008b", "#ff8c00"),
                   alpha=0.5, cex=1.1, fontface=rep(2, 3),
                   fontfamily=rep("sans", 3), cat.fontface=rep(2, 2),
                   cat.fontfamily=rep("sans", 2), margin=0.06,
                   cat.just=list(c(0, 0), c(1, 1)))
dev.off()

# Venn diagrams of each condition intersect with Fo vs. GC DE genes and lncRNAs
png("manuscript_plots/LPS-FovsGC-intersection-Venn.png", width=7.5, height=7.5,
    res=300, units="in")
draw.pairwise.venn(area1=length(lps.intersect$gene_id), area2=length(fo.gc.de.genes),
                   cross.area=length(lps.fo_gc.intersect.genes),
                   category=c("LPS", "Fo vs. GC"),
                   ext.text=F, lwd=c(0, 0), fill=c("#8b008b", "#ff8c00"),
                   alpha=0.5, cex=1.1, fontface=rep(2, 3),
                   fontfamily=rep("sans", 3), cat.fontface=rep(2, 2),
                   cat.fontfamily=rep("sans", 2), margin=0.06,
                   cat.just=list(c(0, 0), c(1, 1)))
dev.off()

# get the genes that intersect with Fo vs. GC for all conditions, then
# intersect and plot those as venn diagrams.  Intersection should = 436 genes
fo.gc.intersect.all <- intersect(intersect(cd40.fo_gc.intersect.genes, lps.fo_gc.intersect.genes),
                                 igm.fo_gc.intersect.genes)

n12 <- intersect(cd40.fo_gc.intersect.genes, lps.fo_gc.intersect.genes)
n23 <- intersect(lps.fo_gc.intersect.genes, igm.fo_gc.intersect.genes)
n13 <- intersect(cd40.fo_gc.intersect.genes, igm.fo_gc.intersect.genes)

png("manuscript_plots/All-FovsGC-intersect-Venn.png", height=8, width=8,
    units="in", res=300)
draw.triple.venn(area1=length(cd40.fo_gc.intersect.genes), area2=length(lps.fo_gc.intersect.genes),
                 area3=length(igm.fo_gc.intersect.genes), n12=length(n12), n23=length(n23),
                 n13=length(n13), n123=length(fo.gc.intersect.all),
                 category=c("CD40L +\nIL-4", "LPS", "IgM"),
                 fill=c("#5d478b", "#ffc125", "#ff0000"), lwd=0, alpha=0.5,
                 cex=1.1, fontface=rep(2, 7), fontfamily=rep("sans", 7),
                 cat.fontfamily=rep("sans", 3), cat.fontface=rep(2, 3),
                 margin=0.06)
dev.off()

########################################################
# Gene ontology enrichment for each condition-specific #
# intersection with Fo -> GC DE genes, as well as      #
# genes that do NOT appear in any in vitro             #
# activation condition                                 #
########################################################

# CD40L + IL-4 vs. Fo -> GC DE genes
cd40.exprs <- read.table("/ifs/projects/proj036/pipeline_tpm_refcoding//tpm.dir/CD40L-000-R1.tpm", sep="\t",
                         h=T, stringsAsFactors=F)
gene.lengths <- cd40.exprs$EffectiveLength
names(gene.lengths) <- cd40.exprs$Name

cd40.gene.vector <- as.integer(names(gene.lengths) %in% cd40.fo_gc.intersect.genes)
names(cd40.gene.vector) <- names(gene.lengths)

pwf.cd40 <- nullp(cd40.gene.vector, "mm10", "ensGene", bias.data = gene.lengths)
GO.wall.cd40 <- goseq(pwf.cd40, genome="mm10", id="ensGene", test.cats=c("GO:BP", "KEGG"),
                      method="Wallenius", use_genes_without_cat = T)

GO.wall.cd40$padjust <- p.adjust(GO.wall.cd40$over_represented_pvalue, method="BH")
GO.wall.cd40$foldEnrich <- ((GO.wall.cd40$numDEInCat/length(cd40.fo_gc.intersect.genes))/(GO.wall.cd40$numInCat/length(gene.lengths)))
print(dim(GO.wall.cd40[GO.wall.cd40$padjust <= 0.05,]) )
head(GO.wall.cd40[GO.wall.cd40$padjust <= 0.05,])

# split into GO and KEGG
GO.results.cd40 <- GO.wall.cd40[!is.na(GO.wall.cd40$ontology),]
KEGG.results.cd40 <- GO.wall.cd40[is.na(GO.wall.cd40$ontology),]

GO.sig.cd40 <- GO.results.cd40[(GO.results.cd40$padjust <= 0.05), ]
go_cats.cd40 <- data.frame(GO.sig.cd40$category, GO.sig.cd40$foldEnrich)

write.table(go_cats.cd40, "/ifs/projects/proj036/R_sessions/manuscript_plots/GO_cats4Revigo-FovsGC-CD40L-intersect.tsv", sep="\t",
            quote=F, row.names = F, col.names = F)

go.revigo.cd40 <- read.table("/ifs/projects/proj036/R_sessions/manuscript_plots/REVIGO-FovsGC-CD40L-intersect.csv", sep=",", h=T)
go.keep.cd40 <- go.revigo.cd40[go.revigo.cd40$eliminated == 0,]$term_ID
GO.revigo.cd40 <- GO.sig.cd40[GO.sig.cd40$category %in% go.keep.cd40,]

top_15.go.cd40 <- GO.revigo.cd40[order(GO.revigo.cd40$padjust, decreasing = F),][1:25,]

foldCol <- colorRampPalette(brewer.pal(9, "YlOrRd"))
foldCol.lo <- foldCol(100)[1]
foldCol.hi <- foldCol(100)[length(foldCol(100))]

p1_go <- ggplot(top_15.go.cd40, aes(x=reorder(term, -top_15.go.cd40$padjust),
                                    y=-log10(padjust),
                                    fill=foldEnrich)) + 
  geom_bar(stat="identity", colour="white") +
  theme_minimal() + coord_flip() + labs(x="GO Category description",
                                        y=expression(paste(-log[10], " P-value"), 
                                                     sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="GO\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="white")
p1_go

#######################
# KEGG enriched terms #
#######################

KEGG.sig.cd40 <- KEGG.results.cd40[(KEGG.results.cd40$padjust <= 0.05), ]

# get the mapping of KEGG Id to pathway name
kegg2id <- as.list(KEGGPATHID2NAME)

# make a vector and attach to KEGG results
kegg_desc <- list()
for(x in 1:length(KEGG.sig.cd40$category)){
  desc <- kegg2id[[KEGG.sig.cd40$category[x]]]
  kegg_desc[[KEGG.sig.cd40$category[x]]] <- desc
}

kegg.vec <- unlist(kegg_desc)
kegg.cats <- names(kegg.vec)
kegg.df.cd40 <- data.frame(cbind(kegg.cats, kegg.vec))
colnames(kegg.df.cd40) <- c("category", "description")

# merge back in with original results
kegg.merge.cd40 <- merge(KEGG.sig.cd40, kegg.df.cd40, by="category")
kegg.merge.cd40 <- kegg.merge.cd40[kegg.merge.cd40$padjust <= 0.01, ]

# let's have a look at that plot
p1_kegg.cd40 <- ggplot(kegg.merge.cd40, 
                       aes(x=reorder(description, -kegg.merge.cd40$padjust),
                           y=-log10(padjust),
                           fill=foldEnrich)) + geom_bar(stat="identity",
                                                           colour="white") + 
  coord_flip() + theme_minimal() + theme(axis.text.x=element_text(colour="black")) +
  labs(x="KEGG Pathway", y=expression(paste(-log[10], " P-value"), sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="KEGG\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="darkgrey")

p1_kegg.cd40

cd40.func_plots <- plot_grid(p1_go, p1_kegg.cd40, nrow = 2, align = "h2", labels=c("GO", "KEGG"))
cd40.func_plots

save_plot(filename="/ifs/projects/proj036/R_sessions/manuscript_plots/Func-enrichment-CD40L-FovsGC-intersect.png", 
          cd40.func_plots,
          nrow = 2, ncol=1, base_aspect_ratio = 1.5, base_width = 9.6)

########################################
# anti-IgM condition vs Fo -> GC genes #
########################################
igm.exprs <- read.table("/ifs/projects/proj036/pipeline_tpm_refcoding//tpm.dir/IgM-001-R1.tpm", sep="\t",
                         h=T, stringsAsFactors=F)
gene.lengths <- igm.exprs$EffectiveLength
names(gene.lengths) <- igm.exprs$Name

igm.gene.vector <- as.integer(names(gene.lengths) %in% igm.fo_gc.intersect.genes)
names(igm.gene.vector) <- names(gene.lengths)

pwf.igm <- nullp(igm.gene.vector, "mm10", "ensGene", bias.data = gene.lengths)
GO.wall.igm <- goseq(pwf.igm, genome="mm10", id="ensGene", test.cats=c("GO:BP", "KEGG"),
                      method="Wallenius", use_genes_without_cat = T)

GO.wall.igm$padjust <- p.adjust(GO.wall.igm$over_represented_pvalue, method="BH")
GO.wall.igm$foldEnrich <- ((GO.wall.igm$numDEInCat/length(igm.fo_gc.intersect.genes))/(GO.wall.igm$numInCat/length(gene.lengths)))
print(dim(GO.wall.cd40[GO.wall.cd40$padjust <= 0.05,]) )
head(GO.wall.igm[GO.wall.igm$padjust <= 0.05,])

# split into GO and KEGG
GO.results.igm <- GO.wall.igm[!is.na(GO.wall.igm$ontology),]
KEGG.results.igm <- GO.wall.igm[is.na(GO.wall.igm$ontology),]

GO.sig.igm <- GO.results.igm[(GO.results.igm$padjust <= 0.05), ]
go_cats.igm <- data.frame(GO.sig.igm$category, GO.sig.igm$foldEnrich)

write.table(go_cats.igm, "/ifs/projects/proj036/R_sessions/manuscript_plots/GO_cats4Revigo-FovsGC-IgM-intersect.tsv", sep="\t",
            quote=F, row.names = F, col.names = F)

go.revigo.igm <- read.table("/ifs/projects/proj036/R_sessions/manuscript_plots/REVIGO-FovsGC-IgM-intersect.csv", sep=",", h=T)
go.keep.igm <- go.revigo.igm[go.revigo.igm$eliminated == 0,]$term_ID
GO.revigo.igm <- GO.sig.igm[GO.sig.igm$category %in% go.keep.igm,]

top_15.go.igm <- GO.revigo.igm[order(GO.revigo.igm$padjust, decreasing = F),][1:25,]

foldCol <- colorRampPalette(brewer.pal(9, "YlOrRd"))
foldCol.lo <- foldCol(100)[1]
foldCol.hi <- foldCol(100)[length(foldCol(100))]

p1_go.igm <- ggplot(top_15.go.igm, aes(x=reorder(term, -top_15.go.igm$padjust),
                                    y=-log10(padjust),
                                    fill=foldEnrich)) + 
  geom_bar(stat="identity", colour="white") +
  theme_minimal() + coord_flip() + labs(x="GO Category description",
                                        y=expression(paste(-log[10], " P-value"), 
                                                     sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="GO\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="white")
p1_go.igm

#######################
# KEGG enriched terms #
#######################

KEGG.sig.igm <- KEGG.results.igm[(KEGG.results.igm$padjust <= 0.05), ]

# get the mapping of KEGG Id to pathway name
kegg2id <- as.list(KEGGPATHID2NAME)

# make a vector and attach to KEGG results
kegg_desc <- list()
for(x in 1:length(KEGG.sig.igm$category)){
  desc <- kegg2id[[KEGG.sig.igm$category[x]]]
  kegg_desc[[KEGG.sig.igm$category[x]]] <- desc
}

kegg.vec <- unlist(kegg_desc)
kegg.cats <- names(kegg.vec)
kegg.df.igm <- data.frame(cbind(kegg.cats, kegg.vec))
colnames(kegg.df.igm) <- c("category", "description")

# merge back in with original results
kegg.merge.igm <- merge(KEGG.sig.igm, kegg.df.igm, by="category")
kegg.merge.igm <- kegg.merge.igm[kegg.merge.igm$padjust <= 0.01,]

# let's have a look at that plot
p1_kegg.igm <- ggplot(kegg.merge.igm,
                       aes(x=reorder(description, -kegg.merge.igm$padjust),
                           y=-log10(padjust),
                           fill=foldEnrich)) + geom_bar(stat="identity",
                                                        colour="white") + 
  coord_flip() + theme_minimal() + theme(axis.text.x=element_text(colour="black")) +
  labs(x="KEGG Pathway", y=expression(paste(-log[10], " P-value"), sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="KEGG\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="darkgrey")

p1_kegg.igm

igm.func_plots <- plot_grid(p1_go.igm, p1_kegg.igm, nrow = 2, align = "h2", labels=c("GO", "KEGG"))
igm.func_plots

save_plot(filename="/ifs/projects/proj036/R_sessions/manuscript_plots/Func-enrichment-IgM-FovsGC-intersect.png", 
          igm.func_plots,
          nrow = 2, ncol=1, base_aspect_ratio = 1.5, base_width = 9.6)

########################################
# LPS condition vs Fo -> GC genes #
########################################
lps.exprs <- read.table("/ifs/projects/proj036/pipeline_tpm_refcoding//tpm.dir/LPS-001-R1.tpm", sep="\t",
                        h=T, stringsAsFactors=F)
gene.lengths <- lps.exprs$EffectiveLength
names(gene.lengths) <- lps.exprs$Name

lps.gene.vector <- as.integer(names(gene.lengths) %in% lps.fo_gc.intersect.genes)
names(lps.gene.vector) <- names(gene.lengths)

pwf.lps <- nullp(lps.gene.vector, "mm10", "ensGene", bias.data = gene.lengths)
GO.wall.lps <- goseq(pwf.lps, genome="mm10", id="ensGene", test.cats=c("GO:BP", "KEGG"),
                     method="Wallenius", use_genes_without_cat = T)

GO.wall.lps$padjust <- p.adjust(GO.wall.lps$over_represented_pvalue, method="BH")
GO.wall.lps$foldEnrich <- ((GO.wall.lps$numDEInCat/length(lps.fo_gc.intersect.genes))/(GO.wall.lps$numInCat/length(gene.lengths)))

# split into GO and KEGG
GO.results.lps <- GO.wall.lps[!is.na(GO.wall.lps$ontology),]
KEGG.results.lps <- GO.wall.lps[is.na(GO.wall.lps$ontology),]

GO.sig.lps <- GO.results.lps[(GO.results.lps$padjust <= 0.05), ]
go_cats.lps <- data.frame(GO.sig.lps$category, GO.sig.lps$foldEnrich)

write.table(go_cats.lps, "/ifs/projects/proj036/R_sessions/manuscript_plots/GO_cats4Revigo-FovsGC-LPS-intersect.tsv", sep="\t",
            quote=F, row.names = F, col.names = F)

go.revigo.lps <- read.table("/ifs/projects/proj036/R_sessions/manuscript_plots/REVIGO-FovsGC-LPS-intersect.csv", sep=",", h=T)
go.keep.lps <- go.revigo.lps[go.revigo.lps$eliminated == 0,]$term_ID
GO.revigo.lps <- GO.sig.lps[GO.sig.lps$category %in% go.keep.lps,]

top_15.go.lps <- GO.revigo.lps[order(GO.revigo.lps$padjust, decreasing = F),][1:25,]

foldCol <- colorRampPalette(brewer.pal(9, "YlOrRd"))
foldCol.lo <- foldCol(100)[1]
foldCol.hi <- foldCol(100)[length(foldCol(100))]

p1_go.lps <- ggplot(top_15.go.lps, aes(x=reorder(term, -top_15.go.lps$padjust),
                                       y=-log10(padjust),
                                       fill=foldEnrich)) + 
  geom_bar(stat="identity", colour="white") +
  theme_minimal() + coord_flip() + labs(x="GO Category description",
                                        y=expression(paste(-log[10], " P-value"), 
                                                     sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="GO\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="white")
p1_go.lps

#######################
# KEGG enriched terms #
#######################

KEGG.sig.lps <- KEGG.results.lps[(KEGG.results.lps$padjust <= 0.05), ]

# get the mapping of KEGG Id to pathway name
kegg2id <- as.list(KEGGPATHID2NAME)

# make a vector and attach to KEGG results
kegg_desc <- list()
for(x in 1:length(KEGG.sig.lps$category)){
  desc <- kegg2id[[KEGG.sig.lps$category[x]]]
  kegg_desc[[KEGG.sig.lps$category[x]]] <- desc
}

kegg.vec <- unlist(kegg_desc)
kegg.cats <- names(kegg.vec)
kegg.df.lps <- data.frame(cbind(kegg.cats, kegg.vec))
colnames(kegg.df.lps) <- c("category", "description")

# merge back in with original results
kegg.merge.lps <- merge(KEGG.sig.lps, kegg.df.lps, by="category")
kegg.merge.lps <- kegg.merge.lps[kegg.merge.lps$padjust <= 0.01,]

# let's have a look at that plot
p1_kegg.lps <- ggplot(kegg.merge.lps,
                      aes(x=reorder(description, -kegg.merge.lps$padjust),
                          y=-log10(padjust),
                          fill=foldEnrich)) + geom_bar(stat="identity",
                                                       colour="white") + 
  coord_flip() + theme_minimal() + theme(axis.text.x=element_text(colour="black")) +
  labs(x="KEGG Pathway", y=expression(paste(-log[10], " P-value"), sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="KEGG\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="darkgrey")

p1_kegg.lps

lps.func_plots <- plot_grid(p1_go.lps, p1_kegg.lps, nrow = 2, align = "h2", labels=c("GO", "KEGG"))
lps.func_plots

save_plot(filename="/ifs/projects/proj036/R_sessions/manuscript_plots/Func-enrichment-LPS-FovsGC-intersect.png", 
          lps.func_plots,
          nrow = 2, ncol=1, base_aspect_ratio = 1.5, base_width = 9.6)





###########################################
# What pathways are enriched in the genes #
# not differentially expressed in the in  #
# vitro experiments                       #
###########################################
fo.gc.diff <- setdiff(fo.gc.de.genes, fo.gc.intersect.all)

fogc.gene.vector <- as.integer(names(gene.lengths) %in% fo.gc.diff)
names(fogc.gene.vector) <- names(gene.lengths)

pwf.fogc <- nullp(fogc.gene.vector, "mm10", "ensGene", bias.data = gene.lengths)
GO.wall.fogc <- goseq(pwf.fogc, genome="mm10", id="ensGene", test.cats=c("GO:BP", "KEGG"),
                     method="Wallenius", use_genes_without_cat = T)

GO.wall.fogc$padjust <- p.adjust(GO.wall.fogc$over_represented_pvalue, method="BH")
GO.wall.fogc$foldEnrich <- ((GO.wall.fogc$numDEInCat/length(fo.gc.diff))/(GO.wall.fogc$numInCat/length(gene.lengths)))

# split into GO and KEGG
GO.results.fogc <- GO.wall.fogc[!is.na(GO.wall.fogc$ontology),]
KEGG.results.fogc <- GO.wall.fogc[is.na(GO.wall.fogc$ontology),]

GO.sig.fogc <- GO.results.fogc[(GO.results.fogc$padjust <= 0.05), ]
go_cats.fogc <- data.frame(GO.sig.fogc$category, GO.sig.fogc$foldEnrich)

write.table(go_cats.fogc, "/ifs/projects/proj036/R_sessions/manuscript_plots/GO_cats4Revigo-FovsGC_only.tsv", sep="\t",
            quote=F, row.names = F, col.names = F)

go.revigo.fogc <- read.table("/ifs/projects/proj036/R_sessions/manuscript_plots/REVIGO-FovsGC-only.csv", sep=",", h=T)
go.keep.fogc <- go.revigo.fogc[go.revigo.fogc$eliminated == 0,]$term_ID
GO.revigo.fogc <- GO.sig.fogc[GO.sig.fogc$category %in% go.keep.fogc,]

top_15.go.fogc <- GO.revigo.fogc[order(GO.revigo.fogc$padjust, decreasing = F),][1:25,]

foldCol <- colorRampPalette(brewer.pal(9, "YlOrRd"))
foldCol.lo <- foldCol(100)[1]
foldCol.hi <- foldCol(100)[length(foldCol(100))]

p1_go.fogc <- ggplot(top_15.go.fogc, aes(x=reorder(term, -top_15.go.fogc$padjust),
                                       y=-log10(padjust),
                                       fill=foldEnrich)) + 
  geom_bar(stat="identity", colour="white") +
  theme_minimal() + coord_flip() + labs(x="GO Category description",
                                        y=expression(paste(-log[10], " P-value"), 
                                                     sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="GO\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="white")
p1_go.fogc


#######################
# KEGG enriched terms #
#######################

KEGG.sig.fogc<- KEGG.results.fogc[(KEGG.results.fogc$padjust <= 0.05), ]

# get the mapping of KEGG Id to pathway name
kegg2id <- as.list(KEGGPATHID2NAME)

# make a vector and attach to KEGG results
kegg_desc <- list()
for(x in 1:length(KEGG.sig.fogc$category)){
  desc <- kegg2id[[KEGG.sig.fogc$category[x]]]
  kegg_desc[[KEGG.sig.fogc$category[x]]] <- desc
}

kegg.vec <- unlist(kegg_desc)
kegg.cats <- names(kegg.vec)
kegg.df.fogc <- data.frame(cbind(kegg.cats, kegg.vec))
colnames(kegg.df.fogc) <- c("category", "description")

# merge back in with original results
kegg.merge.fogc <- merge(KEGG.sig.fogc, kegg.df.fogc, by="category")
kegg.merge.fogc <- kegg.merge.fogc[kegg.merge.fogc$padjust <= 0.01,]

# let's have a look at that plot
p1_kegg.fogc <- ggplot(kegg.merge.fogc,
                      aes(x=reorder(description, -kegg.merge.fogc$padjust),
                          y=-log10(padjust),
                          fill=foldEnrich)) + geom_bar(stat="identity",
                                                       colour="white") + 
  coord_flip() + theme_minimal() + theme(axis.text.x=element_text(colour="black")) +
  labs(x="KEGG Pathway", y=expression(paste(-log[10], " P-value"), sep="")) + 
  scale_fill_gradient(low=foldCol.lo, high=foldCol.hi) +
  guides(fill=guide_legend(title="KEGG\nFold Enrichment")) +
  geom_hline(yintercept=2, linetype="dashed", colour="darkgrey")

p1_kegg.fogc

fogc.func_plots <- plot_grid(p1_go.fogc, p1_kegg.fogc, nrow = 2, align = "h2", labels=c("GO", "KEGG"))
fogc.func_plots

save_plot(filename="/ifs/projects/proj036/R_sessions/manuscript_plots/Func-enrichment-FovsGC-only.png", 
          fogc.func_plots,
          nrow = 2, ncol=1, base_aspect_ratio = 2.5, base_width = 7)

