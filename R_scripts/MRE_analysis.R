########################################################
# Plot and analysis of LNCGme01880 up & down-regulated #
# genes for MREs and MRE sharing                       #
########################################################
library(ggplot2)
library(reshape2)

mre.share <- read.table("/ifs/projects/proj036/pipeline_functional_LNCGme01880_KO/mre.dir/LNCGme01880KO-mre_shared.tsv",
                        sep="\t", h=T, stringsAsFactors=F)

de.res <- read.table("/ifs/projects/proj036/R_sessions/preBcell-DE_results.tsv", h=T, 
                     sep="\t", stringsAsFactors=F)
de.res$padj[is.na(de.res$padj)] <- 1.0
de.res$GeneID <- rownames(de.res)
de.res$Group <- "NoChange"
de.res[de.res$padj <= 0.05 & de.res$log2FoldChange > 0,]$Group <- "Up"
de.res[de.res$padj <= 0.05 & de.res$log2FoldChange < 0,]$Group <- "Down"
table(de.res$Group)

# need gene length data from sailfish!
pre.sailfish <- read.table("/ifs/projects/proj036/pipeline_sailfish_LNCGme01880_KO/transcript_info.tsv",
                           h=F, stringsAsFactors=F, sep="\t")
gene.lengths <- abs(pre.sailfish[, 5] - pre.sailfish[,4])
names(gene.lengths) <- pre.sailfish[, 1]

de.mre.merge <- merge(mre.share, de.res, by.x="gene_id", by.y="GeneID")
table(de.mre.merge$Group)

pre.mre.lengths <- merge(de.mre.merge, pre.sailfish, by.x=c("gene_id"), by.y=c("V1"))
pre.mre.lengths$length <- abs((pre.mre.lengths$V5 - pre.mre.lengths$V4))/1000
pre.mre.lengths$mre_per_kb <- ((pre.mre.lengths$total_shared)/pre.mre.lengths$length)

# scale number of MREs by gene length
ggplot(pre.mre.lengths, aes(x=Group, y=total_shared, fill=Group)) + geom_boxplot()
wilcox.test(x=pre.mre.lengths$total_shared[pre.mre.lengths$Group == "Up"],
            y=pre.mre.lengths$total_shared[pre.mre.lengths$Group == "Down"])

# are the same miRNAs shared with LNCGme01880 between the up and down-regulated genes?
up.mirs <- unique(unlist(strsplit(pre.mre.lengths$shared_miRNAs[(pre.mre.lengths$total_shared >= 7) & (pre.mre.lengths$Group == "Up")], split=",", fixed=T)))
down.mirs <- unique(unlist(strsplit(pre.mre.lengths$shared_miRNAs[(pre.mre.lengths$total_shared >= 7) & (pre.mre.lengths$Group == "Down")], split=",", fixed=T)))
length(up.mirs)
length(down.mirs)
length(intersect(up.mirs, down.mirs))

ggplot(pre.mre.lengths, aes(total_shared, fill=Group)) + geom_density(alpha=0.4) +
  geom_vline(mapping=aes(xintercept=7), linetype="dashed")


