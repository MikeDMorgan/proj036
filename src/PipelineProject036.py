'''
PipelineProjec036.py
====================
 - Classes and functions used in pipeline_project036_timeseries.py
==================================================================

'''
################
# import modules
################

import sys
import glob
import gzip
import os
import itertools
import re
import math
import types
import collections
import time
import optparse, shutil
import sqlite3
import random
import tempfile
import numpy as np
import pandas as pd
from pandas.io import sql
import rpy2.rinterface as rinterface
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as numpy2ri
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineTimeseries as TS
import mygene
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import subprocess
import shlex


###########
# functions
###########


@P.cluster_runnable
def candidateCeRNAs(counts_file,
                    shared_file,
                    lncs_gtf,
                    gene_gtf,
                    annotations,
                    threshold,
                    expr_file,
                    outfile):
    '''
    Select candidate ceRNAs.  Top MRE dense and miRNA sharing
    highly correlated lncRNA: gene pairs.

    Output  heatmaps: expression profile of top lncRNAs and genes and
    normalised cross-correlation matrix
    '''

    counts_df = pd.read_table(counts_file,
                              sep="\t", header=None,
                              index_col=0)
    counts_df.columns = ['MRE_counts']
    counts_df.index.name = 'gene_id'

    shared_df = pd.read_table(shared_file,
                              sep="\t", header=0,
                              index_col=0)
    # only used when selecting a proportion of candidate ceRNAs
    # N = len(shared_df)

    if annotations.endswith("gz"):
        ann_comp = "gzip"
    else:
        ann_comp = None

    annotate_df = pd.read_table(annotations, sep="\t",
                                compression=ann_comp, index_col=0,
                                header=0)

    # take highly correlated annotated gene: lncRNA pair, r>=threshold
    high_cor = annotate_df[annotate_df['value'] >= threshold]
    high_cor['gene_id'] = high_cor.index

    # select all miRNA sharing lncRNA:gene pairs
    top_shared = shared_df.sort(columns='total_shared',
                                ascending=False)
    top_shared = top_shared[top_shared['total_shared'] > 0]
    genes = set(shared_df.index)
    lncs = set(shared_df['lncRNA_id'])

    # calculate gene and lncRNA lengths for MRE density
    gene_file = IOTools.openFile(gene_gtf, "rb")
    lnc_file = IOTools.openFile(lncs_gtf, "rb")

    gene_it = GTF.transcript_iterator(GTF.iterator(gene_file))
    lnc_it = GTF.transcript_iterator(GTF.iterator(lnc_file))
    len_dict = {}

    for git in gene_it:
        for trans in git:
            gid = trans.gene_id
            if gid in genes:
                try:
                    len_dict[gid] += trans.end - trans.start
                except KeyError:
                    len_dict[gid] = trans.end - trans.start
            else:
                pass

    for lit in lnc_it:
        for tran in lit:
            lid = tran.gene_id
            if lid in lncs:
                try:
                    len_dict[lid] += tran.end - tran.start
                except KeyError:
                    len_dict[lid] = tran.end - tran.start
            else:
                pass

    counts_df['length'] = pd.Series(len_dict, dtype=np.float64)
    counts_df['density'] = (counts_df['MRE_counts']/counts_df['length']) * 1000
    lnc_counts = counts_df.loc[lncs]

    # n = len(lnc_counts)
    top_density = lnc_counts.sort(columns='density',
                                  ascending=False)
    dense_lncs = set(top_density.index)
    shared_lncs = set(top_shared['lncRNA_id'])
    shared_genes = set(top_shared.index)
    shared_genes = [sg for sg in shared_genes]

    # select intersection and get annotations
    inter = dense_lncs.intersection(shared_lncs)

    # if none of the intersecting genes/lncRNAs are in the highly correlated
    # set then break out and give appropriate message.
    try:
        high_lnc = high_cor.loc[shared_genes]
    except KeyError:
        E.warn("None of the intersecting genes found in the highly "
               "correlated geneset. Ending here")
        P.touch(outfile)
        return 0

    # remove an NA's
    high_lnc = high_lnc.loc[np.isfinite(high_lnc['value'])]
    high_lnc.index = high_lnc['lncRNA_id']

    if not len(inter):
        P.touch(outfile)
        return 0
    else:
        pass

    # check if highly correlated lncs, MRE dense and those sharing miRNAs
    # intersect - break out if they don't
    if len(inter.intersection(high_lnc.index)):
        top_df = high_lnc.loc[inter]
        top_df = top_df.loc[np.isfinite(top_df['value'])]
        top_df.index = [x for x, y in enumerate(top_df.index)]
    else:
        E.warn("Highly correlated lncRNA set does not intersect those "
               "sharing miRNAs or are MRE dense. Ending here")
        P.touch(outfile)
        return 0        

    candidates = set(top_df['lncRNA_id'])

    # get expression from condition-vst data
    if expr_file.endswith("gz"):
        expr_comp = "gzip"
    else:
        expr_comp = None
    expr_df = pd.read_table(expr_file, sep="\t",
                            header=0, index_col=0,
                            compression=expr_comp)

    cand_expr = expr_df.loc[candidates]
    cand_merge = pd.merge(left=top_df, right=cand_expr,
                          left_on="lncRNA_id", right_index=True,
                          how='inner')

    # get gene symbols from MyGene.info
    mg = mygene.MyGeneInfo()
    try:
        mg_out = mg.querymany(cand_merge['gene_id'].tolist(),
                              scopes="ensemblgene",
                              fields="symbol",
                              species="mouse", returnall=True)['out']
    except AssertionError:
        mg_out = mg.querymany(cand_merge['gene_id'].tolist(),
                              scopes="ensemblgene",
                              fields="symbol",
                              species="mouse", returnall=True)['out']

    mg_df = pd.DataFrame(mg_out)
    symbol_df = pd.merge(left=cand_merge, right=mg_df,
                         left_on='gene_id', right_on='query')
    try:
        symbol_df.drop(['notfound', '_id', 'query'], inplace=True,
                       axis=1)
    except ValueError:
        try:
            symbol_df.drop(['notfound', 'query'], inplace=True,
                           axis=1)
        except ValueError:
            try:
                symbol_df.drop(['_id', 'query'], inplace=True,
                               axis=1)
            except ValueError:
                symbol_df.drop(['query'], inplace=True, axis=1)

    # replace ensembl ids with gene symbols for correlated genes
    E.info("matching ensembl gene ids to gene symbols")
    symbols = symbol_df[['symbol', 'gene_id']]
    symbols.index = symbols['gene_id']
    symbols.drop_duplicates(subset='gene_id', take_last=True, inplace=True)

    cand_df = pd.merge(left=cand_merge, right=symbols,
                       left_on="gene_id", right_on="gene_id",
                       how='inner')

    # remove duplicate entries
    next_df = cand_df.drop_duplicates(subset=['lncRNA_id', 'gene_id'],
                                      take_last=True)

    # select ceRNAs with highest expression, i.e. max(VST expression) >= 8.0
    # only use 0-72 hours
    E.info("filter lncRNAs on expression level")
    expr48_df = expr_df.iloc[:, 1:8]

    # try 48 hour expresion >= 8
    high_expr = expr48_df.apply(max, axis=1) > 8
    next_df.index = next_df['lncRNA_id']
    out_df = next_df.loc[high_expr]

    idxs = [iqx for iqx, iqy in enumerate(out_df.index)]
    out_df.index = idxs

    # check not all lncRNAs have been filtered out
    if len(out_df):
        out_df.to_csv(outfile, sep="\t", index_label="idx")
    else:
        E.warn("These are not the lncRNAs you are looking for."
               "No lncRNAs left in list, relax upstream parameters"
               " to allow for downstream filtering")
        P.touch(outfile)

    if len(set(out_df['lncRNA_id'])) < 2:
        return 0
    else:
        pass

    # generate normalised cross-correlation matrix of candidate ceRNAs
    E.info("Computing normalised cross-correlation matrix")
    lnc_ids = out_df['lncRNA_id'].values
    n_ids = set([lx for lx in lnc_ids])
    cor_df = pd.DataFrame(index=n_ids, columns=n_ids)
    cor_df.fillna(0.0, inplace=True)

    E.info("Calculating cross-correlations "
           "%i calculations" % (len(n_ids)*len(n_ids)))
    for pair in itertools.product(n_ids, n_ids):
        v1 = expr_df.loc[pair[0]].tolist()
        v2 = expr_df.loc[pair[1]].tolist()
        corr = TS.crossCorrelate(v1, v2, lag=0)[0]
        cor_df.loc[pair[0]][pair[1]] = corr

    E.info("Matching ceRNAs to gene symbols")
    new_sym = []
    for k in cor_df.index:
        if re.search("LNC", k):
            new_sym.append(k)
        else:
            vals = symbol_df['symbol'][symbol_df['gene_id'] == k].values[0]
            new_sym.append(vals)

    cor_df.index = new_sym
    cor_df.columns = new_sym

    # remove duplicates
    cand_merge['dups'] = cand_merge.index
    cand_merge.drop_duplicates(subset='dups', take_last=True,
                               inplace=True)
    cand_merge.drop(['dups'], axis=1, inplace=True)
    cand_lncs = cand_merge.index
    r_cor = pandas2ri.py2ri_pandasdataframe(cor_df)
    r_expr = pandas2ri.py2ri_pandasdataframe(cand_merge.iloc[:, 3:])
    r_lncs = ro.StrVector([r for r in cand_lncs])

    R.assign("cor.df", r_cor)
    R.assign("expr.df", r_expr)
    R.assign("r.lncs", r_lncs)

    cond = outfile.split("/")[-1].split("-")[0]

    R('''suppressPackageStartupMessages(library(gplots))''')
    R('''suppressPackageStartupMessages(library(RColorBrewer))''')
    R('''colnames(expr.df) <- c(0, 1, 3, 6, 12, 24, 48, 72, 96, 120)''')
    R('''rownames(expr.df) <- r.lncs''')
    R('''hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(100)''')
    R('''cor_col <- colorRampPalette(brewer.pal(9, "PuOr"))(100)''')

    E.info("Generating  heatmap of ceRNA expression")
    # heatmap of expression
    R('''png("images.dir/%s-candidates-heatmap.png", '''
      '''height=1600, width=1600)''' % cond)
    R('''heatmap.2(as.matrix(expr.df), trace="none", col=hmcol, '''
      '''density.info="none", Colv=colnames(expr.df), margins=c(6,12), '''
      '''dendrogram="none", cexRow=2, cexCol=2)''')
    R('''dev.off()''')

    E.info("Generating heatmap of cross-correlations")
    # heatmap of correlations
    R('''png("images.dir/%s-candidates-correlation_heatmap.png",'''
      '''height=1600, width=1600)''' % cond)
    R('''heatmap.2(as.matrix(cor.df), trace="none", col=cor_col, '''
      '''density.info="none", dendrogram="none", margins=c(10, 10), '''
      '''cexRow=2, cexCol=2)''')
    R('''dev.off()''')


@P.cluster_runnable
def annotateCeRNAs(ceRNA_file,
                   lnc_gtf,
                   mre_file,
                   outfile):
    '''
    Annotate putative ceRNAs from reference gtf file,
    including which MREs are shared.
    '''

    # merge input ceRNAs with shared MRE information
    if mre_file.endswith("gz"):
        mre_comp = "gzip"
    else:
        mre_comp = None

    if ceRNA_file.endswith("gz"):
        cerna_comp = "gzip"
    else:
        cerna_comp = None

    # handle bug when all ceRNAs are filtered out
    if IOTools.isEmpty(ceRNA_file):
        E.warn("No ceRNAs to parse, exiting")
        P.touch(outfile)
        return 0
    else:
        pass

    mre_df = pd.read_table(mre_file, sep="\t", header=0,
                           index_col=None, compression=mre_comp)

    cerna_df = pd.read_table(ceRNA_file, sep="\t", header=0,
                             index_col=0, compression=cerna_comp)

    matched = pd.merge(left=cerna_df, right=mre_df,
                       left_on=["gene_id", "lncRNA_id", "symbol"],
                       right_on=["gene_id", "lncRNA_id", "symbol"],
                       how='inner')

    matched = matched.drop_duplicates(subset=['gene_id', 'lncRNA_id', 'value'],
                                      take_last=True)

    # drop ceRNA pairs that do not share any miRNAs
    shared_df = matched[matched['total_shared'] != 0]
    lnc_open = IOTools.openFile(lnc_gtf, "rb")
    lnc_it = GTF.transcript_iterator(GTF.iterator(lnc_open))
    lncs = set(shared_df['lncRNA_id'])

    # get genome coordinates of ceRNAs
    lnc_dict = {}

    for lit in lnc_it:
        for lnc in lit:
            if lnc.gene_id in lncs:
                try:
                    lnc_dict[lnc.gene_id]['contig'] = lnc.contig
                    lnc_dict[lnc.gene_id]['start'].append(lnc.start)
                    lnc_dict[lnc.gene_id]['end'].append(lnc.end)
                    lnc_dict[lnc.gene_id]['strand'].add(lnc.strand)
                except KeyError:
                    lnc_dict[lnc.gene_id] = {'contig': set(lnc.contig),
                                             'start': [lnc.start],
                                             'end': [lnc.end],
                                             'strand': set(lnc.strand)}
            else:
                pass

    shared_df.index = shared_df['lncRNA_id']

    lcontig_dict = {}
    lstart_dict = {}
    lend_dict = {}
    lstrand_dict = {}

    for l_each in lnc_dict.keys():
        lcontig_dict[l_each] = lnc_dict[l_each]['contig']
        lstart_dict[l_each] = min(lnc_dict[l_each]['start'])
        lend_dict[l_each] = max(lnc_dict[l_each]['end'])
        lstrand_dict[l_each] = [ls for ls in lnc_dict[l_each]['strand']][-1]

    shared_df.loc[:, 'lncRNA_contig'] = pd.Series(lcontig_dict,
                                                  index=shared_df.index)
    shared_df.loc[:, 'lncRNA_start'] = pd.Series(lstart_dict,
                                                 index=shared_df.index)
    shared_df.loc[:, 'lncRNA_end'] = pd.Series(lend_dict,
                                               index=shared_df.index)
    shared_df.loc[:, 'lncRNA_strand'] = pd.Series(lstrand_dict,
                                                  index=shared_df.index)

    shared_df.sort(['gene_id', 'lncRNA_contig'],
                   ascending=True, inplace=True)
    shared_df.index = [x for x, y in enumerate(shared_df.index)]
    try:
        shared_df.drop(['_id'], inplace=True, axis=1)
    except ValueError:
        pass

    shared_df.to_csv(outfile, sep="\t", index_label="indx")


def netExpressionFiltering(expression,
                           direction):
    '''
    Calculate net expression over a time course, first 48-72hours
    as mean fold change between consecutive time points > 1
    Alternative: check sum of fold changes > 1
    '''

    timepoints = expression.index
    change_list = []
    for t in range(0, len(timepoints)):
        if t == 0:
            # t0 = expression.loc[timepoints[t]]
            pass
        else:
            t1 = expression.loc[timepoints[t-1]]
            t2 = expression.loc[timepoints[t]]
            fchange = t2/t1
            change_list.append(fchange)

    net_change = np.mean(change_list)

    if direction == "up":
        if net_change >= 1:
            return 1
        else:
            return 0
    elif direction == "down":
        if net_change <= 1:
            return 1
        else:
            return 0
    else:
        raise ValueError("unknown direction of expression change"
                         " on which to filter")


@P.cluster_runnable
def filterAnnotatedCeRNAs(cerna_file,
                          expression_file,
                          direction,
                          outfile):
    '''
    Filter annotated ceRNAs based on direction of expression,
    i.e. up-regulated or down-regulated
    '''

    # handle bug when all ceRNAs are filtered out
    if IOTools.isEmpty(cerna_file):
        E.warn("No ceRNAs to parse, exiting")
        P.touch(outfile)
        return 0
    else:
        pass

    if expression_file.endswith("gz"):
        expr_comp = "gzip"
    else:
        expr_comp = None

    expr_df = pd.read_table(expression_file, sep="\t", header=0,
                            index_col=0, compression=expr_comp)

    cerna_df = pd.read_table(cerna_file, sep="\t", header=0,
                             index_col=0)

    cerna_df.index = cerna_df['lncRNA_id']

    filter_dict = {}
    expr_df = expr_df.iloc[:, 1:7]
    for lnc in cerna_df.index:
        expr = expr_df.loc[lnc]
        if netExpressionFiltering(expr, direction):
            filter_dict[lnc] = cerna_df.loc[lnc]
        else:
            pass

    # may return a dictionary of dataframes, append all together
    lnc_keys = filter_dict.keys()
    if len(lnc_keys):
        df0 = filter_dict[lnc_keys[0]]
        lnc_keys.pop(0)
        for lncdf in lnc_keys:
            df1 = filter_dict[lncdf]

            df0 = df0.append(df1)
        filter_df = df0
        filter_df.index = [x for x, y in enumerate(filter_df.index)]
        filter_df.to_csv(outfile, sep="\t", index_label="indx")
    else:
        E.warn("No ceRNAs to filter on expression direction.")
        P.touch(outfile)


@P.cluster_runnable
def annotateGeneExpr(cernas, expression, outfile):
    '''
    Annotate ceRNAs with expressiond data for partner protein-coding genes
    '''

    if cernas.endswith("gz"):
        cerna_comp = "gzip"
    else:
        cerna_comp = None

    if expression.endswith("gz"):
        expr_comp = "gzip"
    else:
        expr_comp = None

    try:
        cerna_df = pd.read_table(cernas, sep="\t", header=0, index_col=0,
                                 compression=cerna_comp)
    except ValueError:
        E.warn("no ceRNA candidates to parse. Exiting")
        P.touch(outfile)
        return 0

    expr_df = pd.read_table(expression, sep="\t", header=0, index_col=0,
                            compression=expr_comp)

    genes = cerna_df['gene_id'].tolist()
    gene_expr = expr_df.loc[genes]
    cerna_lncs = cerna_df[['lncRNA_id', 'value',
                           'gene_id', 'symbol', 'shared_miRNAs']]

    merge_df = pd.merge(left=cerna_lncs, right=gene_expr, left_on='gene_id',
                        right_index=True, how='inner')
    merge_df.drop_duplicates(subset=['gene_id', 'lncRNA_id', 'symbol'],
                             take_last=True, inplace=True)
    merge_df.to_csv(outfile, index_label="indx", sep="\t")


@P.cluster_runnable
def mergeStatTable(file_list, outfile):
    '''
    merge statistical testing tables
    '''

    df = pd.read_table(file_list[0], sep="\t", header=0, index_col=0)
    file_list.pop(0)
    for fle in file_list:
        df_ = pd.read_table(fle, sep="\t", header=0, index_col=0)
        df = df.append(df_)
    df.to_csv(outfile, sep="\t", index_label="comparison")


@P.cluster_runnable
def annotateGeneList(infile,
                     lnc_gtf,
                     outfile):
    '''
    Annotate a list of correlated genes and lncRNAs with gene
    symbols and lncRNA genome co-ordinates
    '''

    if infile.endswith("gz"):
        comp = "gzip"
    else:
        comp = None

    cor_df = pd.read_table(infile, sep="\t", header=0,
                           index_col=None, compression=comp)

    gene_set = cor_df['gene_id'].tolist()
    genes = [g for g in set(gene_set)]
    mg = mygene.MyGeneInfo()
    # can throw AssertionError if too many queries at once
    # a single retry usually works
    try:
        mg_out = mg.querymany(genes, scopes="ensemblgene", fields="symbol",
                              species="mouse", returnall=True)['out']
    except AssertionError:
        mg_out = mg.querymany(genes, scopes="ensemblgene", fields="symbol",
                              species="mouse", returnall=True)['out']
    mg_df = pd.DataFrame(mg_out)
    mg_df.drop_duplicates(subset="query", take_last=True, inplace=True)
    merged = pd.merge(left=cor_df, right=mg_df,
                      how='left', left_on="gene_id",
                      right_on='query')
    try:
        merged.drop(['_id', 'notfound', 'query'], inplace=True, axis=1)
    except ValueError:
        try:
            merged.drop(['_id', 'query'], inplace=True, axis=1)
        except ValueError:
            merged.drop(['query'], inplace=True, axis=1)

    # get lncRNA co-ordinates from file
    lnc_dict = {}
    lnc_file = IOTools.openFile(lnc_gtf, "rb")
    lnc_it = GTF.transcript_iterator(GTF.iterator(lnc_file))
    for gene in lnc_it:
        start = []
        end = []
        for trans in gene:
            start.append(trans.start)
            end.append(trans.end)
            strand = trans.strand
            contig = trans.contig
            gid = trans.gene_id
            lnc_class = trans.source
            exon_class = trans.asDict()['exon_status']
        lnc_dict[gid] = {'contig': contig,
                         'start': min(start),
                         'end': max(end),
                         'strand': strand,
                         'lnc_class': lnc_class,
                         'exon_class': exon_class}
    lnc_file.close()

    lnc_df = pd.DataFrame(lnc_dict).T
    try:
        annotated = pd.merge(left=merged, right=lnc_df,
                             left_on='lncRNA', right_index=True,
                             how='left')
    except KeyError:
        annotated = pd.merge(left=merged, right=lnc_df,
                             left_on='lncRNA_id', right_index=True,
                             how='left')

    # drop cluster info if contained in dataframe
    try:
        annotated.drop(['gene_cluster', 'lncRNA_cluster'],
                       inplace=True, axis=1)
    except ValueError:
        pass
    columns = ['ensembl_id', 'lncRNA_id', 'correlation',
               'gene_symbol', 'chromosome', 'lnc_end',
               'lnc_exonic', 'lnc_class', 'lnc_start',
               'lnc_strand']
    annotated.columns = columns

    sort_cols = ['ensembl_id', 'gene_symbol',
                 'correlation', 'lncRNA_id',
                 'chromosome', 'lnc_start', 'lnc_end', 'lnc_strand',
                 'lnc_class', 'lnc_exonic']
    annotated = annotated[sort_cols]
    annotated.to_csv(outfile, sep="\t",
                     index_label="idx")


@P.cluster_runnable
def lncsPerGene(infile, outfile):
    '''
    count the number of lncRNAs correlated with each protein-coding gene.
    '''

    if infile.endswith("gz"):
        comp = "gzip"
    else:
        comp = None

    df = pd.read_table(infile, sep="\t", header=0,
                       index_col=0, compression=comp)
    genes = set(df['ensembl_id'])
    gdict = {}
    for gene in genes:
        gdf = df[df['ensembl_id'] == gene]
        cor_lncs = set(gdf['lncRNA_id'])
        gdict[gene] = len(cor_lncs)

    gser = pd.Series(gdict)
    gser.columns = ['counts']
    gser.to_csv(outfile, sep="\t", index_label="ensembl_id")


@P.cluster_runnable
def genesPerLnc(infile, outfile):
    '''
    count the number of genes correlated with each lncRNAs,
    subset by lncRNA class
    '''

    if infile.endswith("gz"):
        comp = "gzip"
    else:
        comp = None
    df = pd.read_table(infile, sep="\t", header=0,
                       index_col=0, compression=comp)

    lncs = set(df['lncRNA_id'])
    ldict = {}
    for lnc in lncs:
        ldf = df[df['lncRNA_id'] == lnc]
        lnc_class = [str(c) for c in set(ldf['lnc_class'])][0]
        cor_genes = set(ldf['ensembl_id'])
        ldict[lnc] = {'n_genes': len(cor_genes),
                      'lnc_class': lnc_class}
    sum_df = pd.DataFrame(ldict).T
    sum_df.to_csv(outfile, sep="\t", index_label="lncRNA_id")


@P.cluster_runnable
def plotGeneCounts(infile, outfile):
    '''
    Plot counts of lncRNAs per gene
    '''

    df = pd.read_table(infile, sep="\t", header=None, index_col=0)
    df.columns = ["counts"]
    r_df = pandas2ri.py2ri_pandasdataframe(df)
    R.assign("g.df", r_df)

    R('''library(ggplot2)''')
    R('''g.df$counts <- as.numeric(as.character(g.df$counts))''')
    R('''p_gene <- ggplot(g.df, aes(x=counts)) + '''
      '''geom_histogram(colour="black", fill="coral", binwidth=5) + '''
      '''labs(x="N correlated lncRNAs per gene")''')
    R('''png("%s", height=480, width=480)''' % outfile)
    R('''print(p_gene)''')
    R('''dev.off()''')


@P.cluster_runnable
def plotLncCounts(infile, outfile):
    '''
    Plot ggplot histogram of number of lncRNAs correlated per gene
    '''

    df = pd.read_table(infile, sep="\t", header=0, index_col=None)

    r_df = pandas2ri.py2ri_pandasdataframe(df)
    R.assign("r.df", r_df)
    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''r.df$lnc_class <- as.factor(r.df$lnc_class)''')
    R('''p_lnc <- ggplot(r.df, aes(x=n_genes, colour=lnc_class, '''
      '''fill=lnc_class)) + geom_density(alpha=0.5) + '''
      '''facet_grid(lnc_class ~ .) + '''
      '''labs(x="N correlated genes per lncRNA")''')
    R('''png('%s', height=480, width=480)''' % outfile)
    R('''print(p_lnc)''')
    R('''dev.off()''')


@P.cluster_runnable
def plotCorrelations(cor_file, rand_file, prox_file, outfile):
    '''
    Plot distributions of correlation coefficients across gene sets
    '''

    cor_df = pd.read_table(cor_file, sep="\t", index_col=None,
                           header=0)
    cor_df.columns = ['gene_id', 'correlation', 'lncRNA_id']

    prox_df = pd.read_table(prox_file, sep="\t", index_col=None,
                            header=0)
    rand_df = pd.read_table(rand_file, sep="\t", index_col=None,
                            header=0)
    anti_df = pd.read_table(anti_file, sep="\t", index_col=None,
                            header=0)
    cor_df['cat'] = "correlated"
    prox_df['cat'] = "proximal"
    rand_df['cat'] = "random"
    anti_df['cat'] = "anticorrelated"

    all_df = cor_df.append(prox_df)
    all_df = all_df.append(rand_df)
    all_df = all_df.append(anti_df)
    all_idx = [i for i, t in enumerate(all_df.index)]
    all_df.index = all_idx

    r_all = pandas2ri.py2ri_pandasdataframe(all_df)
    R.assign("r.df", r_all)
    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''r.df$cat <- as.factor(r.df$cat)''')
    R('''p_cor <- ggplot(r.df, aes(x=correlation, fill=cat, colour=cat)) + '''
      '''geom_density(alpha=0.2) + labs(x="Normalised cross-correlation") + '''
      '''theme_bw() + theme(text=element_text(size=18, colour="black")) + '''
      '''guides(colour=guide_legend(title="lncRNA set"), '''
      '''fill=guide_legend(title="lncRNA set"))''')
    R('''png("%s", height=480, width=480)''' % outfile)
    R('''print(p_cor)''')
    R('''dev.off()''')


@P.cluster_runnable
def plotSigMREs(cor, prox, random, outfile):
    '''
    Density plot of correlated, proximal and random lncRNA:gene
    pairs proprotion of statisticall significantly conserved MREs
    '''

    if cor.endswith("gz"):
        cor_comp = "gzip"
    else:
        cor_comp = None

    if prox.endswith("gz"):
        prox_comp = "gzip"
    else:
        prox_comp = None

    if random.endswith("gz"):
        rand_comp = "gzip"
    else:
        rand_comp = None

    cor_df = pd.read_table(cor, sep="\t", index_col=None,
                           header=0, compression=cor_comp,
                           comment='#')
    prox_df = pd.read_table(prox, sep="\t", index_col=None,
                            header=0, compression=prox_comp,
                            comment='#')
    rand_df = pd.read_table(random, sep="\t", index_col=None,
                            header=0, compression=rand_comp,
                            comment='#')

    cor_idx = cor_df.index
    lnc_tag = [x for x in cor_idx if re.search("LNC", cor_df.loc[x]['target'])]
    lnc_cor_df = cor_df.loc[lnc_tag]
    g_tag = [g for g in cor_idx if re.search("ENS", cor_df.loc[g]['target'])]
    gene_cor_df = cor_df.loc[g_tag]

    pro_idx = prox_df.index
    pro_ln = [p for p in pro_idx if re.search("LNC", prox_df.loc[p]['target'])]
    lnc_prox_df = prox_df.loc[pro_ln]
    pro_g = [w for w in pro_idx if re.search("ENS", prox_df.loc[w]['target'])]
    gene_prox_df = prox_df.loc[pro_g]

    ran_idx = rand_df.index
    ran_ln = [r for r in ran_idx if re.search("LNC", rand_df.loc[r]['target'])]
    lnc_rand_df = rand_df.loc[ran_ln]
    ran_g = [d for d in ran_idx if re.search("ENS", rand_df.loc[d]['target'])]
    gene_rand_df = rand_df.loc[ran_g]

    gprox_df = summMreCons(gene_prox_df)
    lprox_df = summMreCons(lnc_prox_df)
    grand_df = summMreCons(gene_rand_df)
    lrand_df = summMreCons(lnc_rand_df)
    gcor_df = summMreCons(gene_cor_df)
    lcor_df = summMreCons(lnc_cor_df)

    gprox_df['biotype'] = "gene"
    lprox_df['biotype'] = "lncRNA"
    grand_df['biotype'] = "gene"
    lrand_df['biotype'] = "lncRNA"
    gcor_df['biotype'] = "gene"
    lcor_df['biotype'] = "lncRNA"

    cor_sum = gcor_df.append(lcor_df)
    cor_sum['cat'] = "correlated"
    lcor_df['cat'] = "correlated"
    gcor_df['cat'] = "correlated"

    prox_sum = gprox_df.append(lprox_df)
    prox_sum['cat'] = "proximal"
    gprox_df['cat'] = "proximal"
    lprox_df['cat'] = "proximal"

    rand_sum = grand_df.append(lrand_df)
    rand_sum['cat'] = "random"
    grand_df['cat'] = "random"
    lrand_df['cat'] = "random"

    all_sum = cor_sum.append(prox_sum)
    all_sum = all_sum.append(rand_sum)
    all_idx = [ix for ix, y in enumerate(all_sum.index)]
    all_sum.index = all_idx
    r_sum = pandas2ri.py2ri_pandasdataframe(all_sum)
    R.assign("r.sum", r_sum)

    all_lncs = lcor_df.append(lprox_df)
    all_lncs = all_lncs.append(lrand_df)
    l_idx = [lx for lx, p in enumerate(all_lncs.index)]
    all_lncs.index = l_idx
    r_lncs = pandas2ri.py2ri_pandasdataframe(all_lncs)
    R.assign("r.lncs", r_lncs)

    all_genes = gcor_df.append(gprox_df)
    all_genes = all_genes.append(grand_df)
    g_idx = [gx for gx, b in enumerate(all_genes.index)]
    all_genes.index = g_idx
    # formally test differences between gene sets

    wilcoxpy = R['wilcox.test']
    test_dict = {}
    for combs in itertools.combinations(set(all_genes['cat']), r=2):
        pos1 = [x for x in combs if re.search("correlated", x)]
        pos2 = [q for q in combs if re.search("random", q)]

        if not len(pos1):
            pos1 = [p for p in combs if re.search("proximal", p)]
        elif not pos2:
            pos2 = [j for j in combs if re.search("proximal", j)]
        vec1 = all_genes['prop_sig'][all_genes['cat'] == pos1[0]]
        r_vec1 = ro.FloatVector([r for r in vec1])
        vec2 = all_genes['prop_sig'][all_genes['cat'] == pos2[0]]
        r_vec2 = ro.FloatVector([g for g in vec2])
        res = wilcoxpy(r_vec1, r_vec2, alternative="greater")
        pval = res.rx('p.value')[0][0]
        stat = res.rx('statistic')[0][0]
        test_dict[(pos1[0], pos2[0])] = {"W": stat,
                                         "p-value": pval}
    test_table = pd.DataFrame(test_dict).T
    cond = cor.split("/")[-1].split("-")[0]
    test_table['condition'] = cond
    test_table.to_csv("stats.dir/%s-sig_conserved-stats.tsv" % cond,
                      sep="\t", index_label="reference")

    r_genes = pandas2ri.py2ri_pandasdataframe(all_genes)
    R.assign("r.genes", r_genes)

    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''r.lncs$cat <- as.factor(r.lncs$cat)''')
    R('''p_sums <- ggplot(r.lncs, aes(x=prop_sig, fill=cat, colour=cat)) + '''
      '''geom_density(alpha=0.2) + '''
      '''labs(x="Proportion of significantly conserved MREs", y="density", '''
      '''title="Proportion of statistically significantly conserved\n'''
      '''MREs over lncRNAs")''')
    R('''png('%s', height=480, width=480)''' % outfile)
    R('''print(p_sums)''')
    R('''dev.off()''')


def summMreCons(dataframe):
    '''
    summarise over a dataframe of MRE conservation scores and
    p-values
    '''

    df_dict = {}
    for entry in set(dataframe['target']):
        edict = {}
        e_df = dataframe[dataframe['target'] == entry]
        mean = np.mean(e_df['score'])
        n_mres = len(e_df['score'])
        n_sig = len([s for s in e_df['pCons'] if s < 0.01])
        try:
            prop_sig = (n_sig/float(n_mres)) * 100
        except ZeroDivisionError:
            prop_sig = 0.0

        edict = {'mean_phastCons': mean,
                 'n_MREs': n_mres,
                 'n_sigP': n_sig,
                 'prop_sig': prop_sig}
        df_dict[entry] = edict

    out_df = pd.DataFrame(df_dict).T
    return out_df


@P.cluster_runnable
def annotateMreGTF(lnc_gtf, mre_gtf):
    '''
    annotated a gtf/gff file of predicted MREs with information
    from  matched lncRNA/gene gtf
    '''

    lnc_index = GTF.readAndIndex(GTF.iterator(IOTools.openFile(lnc_gtf)))
    ofile = IOTools.openFile(mre_gtf, "rb")
    for mre in GTF.iterator(ofile):
        lnc_source = lnc_index.get(mre.contig, mre.start, mre.end)
        for i in lnc_source:
            mre.source = i[2].source
            yield mre


@P.cluster_runnable
def countMREsOverLncs(mre_gtf, lnc_gtf, outfile):
    '''
    Count the number of non-redundant MREs overlapping
    lncRNA gene models in input
    count over whole genes not just transcripts
    '''

    indices = {}
    gtf_it = GTF.readAndIndex(GTF.iterator(IOTools.openFile(lnc_gtf, "rb")))
    indices['gene_id'] = gtf_it
    trans_gene_dict = {}

    counter_dict = {}

    for lnc in GTF.iterator(IOTools.openFile(lnc_gtf)):
        counter_dict[lnc.gene_id] = 0
        trans_gene_dict[lnc.transcript_id] = lnc.gene_id

    with IOTools.openFile(mre_gtf, "rb") as mre_open:
        mre_it = GTF.iterator_filtered(GTF.iterator(mre_open),
                                       feature="MRE")
        for mre in mre_it:
            overlap = mre.asDict()['target']
            gene_id = trans_gene_dict[overlap]
            counter_dict[gene_id] += 1

    with IOTools.openFile(outfile, "w") as ofile:
        for x in counter_dict.keys():
            ofile.write("%s\t%i\n" % (x,
                                      counter_dict[x]))


@P.cluster_runnable
def plotMreDensity(cor_file,
                   prox_file,
                   rand_file,
                   anti_file,
                   ref_gtf,
                   lnc_gtf,
                   outfile):
    '''
    Plot MRE density as number of MREs per nt over lncRNAs
    '''

    if cor_file.endswith("gz"):
        cor_comp = "gzip"
    else:
        cor_comp = None

    if prox_file.endswith("gz"):
        prox_comp = "gzip"
    else:
        prox_comp = None

    if rand_file.endswith("gz"):
        rand_comp = "gzip"
    else:
        rand_comp = None

    if anti_file.endswith("gz"):
        anti_comp = "gzip"
    else:
        anti_comp = None

    cor_df = pd.read_table(cor_file, sep="\t", index_col=0,
                           header=0, compression=cor_comp)
    cor_df.columns = ['MRE_counts']
    cor_df.index.name = 'gene_id'
    cor_index = cor_df.index.tolist()
    cor_lncs = set([cl for cl in cor_index if re.search("LNC", cl)])
    cor_genes = set([cg for cg in cor_index if re.search("ENS", cg)])

    prox_df = pd.read_table(prox_file, sep="\t", index_col=0,
                            header=None, compression=prox_comp)
    prox_df.index.name = 'gene_id'
    prox_df.columns = ['MRE_counts']
    prox_index = prox_df.index.tolist()
    prox_lncs = set([pl for pl in prox_index if re.search("LNC", pl)])
    prox_genes = set([pg for pg in prox_index if re.search("ENS", pg)])

    rand_df = pd.read_table(rand_file, sep="\t", index_col=0,
                            header=None, compression=rand_comp)
    rand_df.index.name = 'gene_id'
    rand_df.columns = ['MRE_counts']

    rand_index = rand_df.index.tolist()
    rand_lncs = set([rl for rl in rand_index if re.search("LNC", rl)])
    rand_genes = set([rg for rg in rand_index if re.search("ENS", rg)])

    anti_df = pd.read_table(anti_file, sep="\t", index_col=0,
                            header=0, compression=anti_comp)
    anti_df.index.name = 'gene_id'
    anti_df.columns = ['MRE_counts']

    anti_index = anti_df.index.tolist()
    anti_lncs = set([al for al in anti_index if re.search("LNC", al)])
    anti_genes = set([ag for ag in anti_index if re.search("ENS", ag)])
    
    cor_len_dict = {}
    prox_len_dict = {}
    rand_len_dict = {}
    anti_len_dict = {}

    gene_file = IOTools.openFile(ref_gtf, "rb")
    gene_iterator = GTF.transcript_iterator(GTF.iterator(gene_file))
    lnc_file = IOTools.openFile(lnc_gtf, "rb")
    lnc_iterator = GTF.transcript_iterator(GTF.iterator(lnc_file))

    # get all gene and lncRNA lengths
    for git in gene_iterator:
        for trans in git:
            gid = trans.gene_id
            if gid in cor_genes and gid in prox_genes and gid not in anti_genes:
                try:
                    cor_len_dict[gid] += (trans.end - trans.start)
                    prox_len_dict[gid] += (trans.end - trans.start)
                except KeyError:
                    cor_len_dict[gid] = trans.end - trans.start
                    prox_len_dict[gid] = trans.end - trans.start
            elif gid in cor_genes and gid not in prox_genes and gid not in anti_genes:
                try:
                    cor_len_dict[gid] += (trans.end - trans.start)
                except KeyError:
                    cor_len_dict[gid] = trans.end - trans.start              
            elif gid in rand_genes and gid in prox_genes and gid in anti_genes:
                try:
                    rand_len_dict[gid] += trans.end - trans.start
                    prox_len_dict[gid] += trans.end - trans.start
                    anti_len_dict[gid] += trans.end - trans.start
                except KeyError:
                    rand_len_dict[gid] = trans.end - trans.start
                    prox_len_dict[gid] = trans.end - trans.start
                    anti_len_dict[gid] = trans.end - trans.start
            elif gid in rand_genes and gid in anti_genes and gid not in prox_genes:
                try:
                    rand_len_dict[gid] += trans.end - trans.start
                    anti_len_dict[gid] += trans.end - trans.start
                except KeyError:
                    rand_len_dict[gid] = trans.end - trans.start
                    anti_len_dict[gid] = trans.end - trans.start
            elif gid in prox_genes and gid in anti_genes and gid not in rand_genes:
                try:
                    prox_len_dict[gid] += trans.end - trans.start
                    anti_len_dict[gid] += trans.end - trans.start
                except KeyError:
                    prox_len_dict[gid] = trans.end - trans.start
                    anti_len_dict[gid] = trans.end - trans.start
            elif gid in prox_genes and gid in rand_genes and gid not in anti_genes:
                try:
                    prox_len_dict[gid] += trans.end - trans.start
                    rand_len_dict[gid] += trans.end - trans.start
                except KeyError:
                    prox_len_dict[gid] = trans.end - trans.start
                    rand_len_dict[gid] = trans.end - trans.start
            elif gid in prox_genes and gid not in rand_genes and gid not in anti_genes:
                try:
                    prox_len_dict[gid] += trans.end - trans.start
                except KeyError:
                    prox_len_dict[gid] = trans.end - trans.start
            elif gid in rand_genes and gid not in prox_genes and gid not in anti_genes:
                try:
                    rand_len_dict[gid] += trans.end - trans.start
                except KeyError:
                    rand_len_dict[gid] = trans.end - trans.start
            elif gid in anti_genes and gid not in prox_genes and gid not in rand_genes:
                try:
                    anti_len_dict[gid] += trans.end - trans.start
                except KeyError:
                    anti_len_dict[gid] = trans.end - trans.start
            else:
                pass
    gene_file.close()

    for lit in lnc_iterator:
        for tran in lit:
            lid = tran.gene_id
            if lid in cor_lncs:
                try:
                    cor_len_dict[lid] += tran.end - tran.start
                    prox_len_dict[lid] += tran.end - tran.start
                    anti_len_dict[lid] += tran.end - tran.start
                    rand_len_dict[lid] += tran.end - tran.start
                except KeyError:
                    cor_len_dict[lid] = tran.end - tran.start
                    prox_len_dict[lid] = tran.end - tran.start
                    anti_len_dict[lid] = tran.end - tran.start
                    rand_len_dict[lid] = tran.end - tran.start
            elif lid in rand_lncs:
                try:
                    rand_len_dict[lid] += tran.end - tran.start
                    prox_len_dict[lid] += tran.end - tran.start
                    anti_len_dict[lid] += tran.end - tran.start
                except KeyError:
                    rand_len_dict[lid] = tran.end - tran.start
                    prox_len_dict[lid] = tran.end - tran.start
                    anti_len_dict[lid] = tran.end - tran.start
            elif lid in anti_lncs:
                try:
                    prox_len_dict[lid] += tran.end - tran.start
                    anti_len_dict[lid] += tran.end - tran.start
                except KeyError:
                    prox_len_dict[lid] = tran.end - tran.start
                    anti_len_dict[lid] = tran.end - tran.start
            elif lid in prox_lncs:
                try:
                    prox_len_dict[lid] += tran.end - tran.start
                except KeyError:
                    prox_len_dict[lid] = tran.end - tran.start
            else:
                pass
    lnc_file.close()

    cor_df['length'] = pd.Series(cor_len_dict, dtype=np.float64)
    prox_df['length'] = pd.Series(prox_len_dict, dtype=np.float64)
    rand_df['length'] = pd.Series(rand_len_dict, dtype=np.float64)
    anti_df['length'] = pd.Series(anti_len_dict, dtype=np.float64)

    cor_df['density'] = cor_df['MRE_counts']/(cor_df['length']/1000.0)
    prox_df['density'] = prox_df['MRE_counts']/(prox_df['length']/1000.0)
    rand_df['density'] = rand_df['MRE_counts']/(rand_df['length']/1000.0)
    anti_df['density'] = anti_df['MRE_counts']/(anti_df['length']/1000.0)

    cor_lnc_counts = cor_df.loc[cor_lncs]
    prox_lnc_counts = prox_df.loc[prox_lncs]
    rand_lnc_counts = rand_df.loc[rand_lncs]
    anti_lnc_counts = anti_df.loc[anti_lncs]

    cor_lnc_counts['cat'] = "correlated"
    prox_lnc_counts['cat'] = "proximal"
    rand_lnc_counts['cat'] = "random"
    anti_lnc_counts['cat'] = "anticorrelated"

    cor_lnc_counts['group'] = "correlated/proximal"
    prox_lnc_counts['group'] = "correlated/proximal"
    rand_lnc_counts['group'] = "random"
    anti_lnc_counts['group'] = "anticorrelated"

    all_lnc_frame = cor_lnc_counts.append(prox_lnc_counts)
    all_lnc_frame = all_lnc_frame.append(rand_lnc_counts)
    all_lnc_frame = all_lnc_frame.append(anti_lnc_counts)
    all_lnc_frame.index = [ix for ix, iy in enumerate(all_lnc_frame.index)]

    # break if all counts are zero or < 10 objects
    if max(all_lnc_frame['MRE_counts']) == 0 or len(all_lnc_frame) < 10:
        P.touch(outfile)
        return 0
    else:
        pass

    pandas2ri.activate()
    r_lnc_df = pandas2ri.py2ri_pandasdataframe(all_lnc_frame)
    # formally test differences between gene sets with wilcoxon test
    wilcoxpy = R['wilcox.test']
    test_dict = {}
    for combs in itertools.combinations(set(all_lnc_frame['group']), r=2):
        vec1 = all_lnc_frame['density'][all_lnc_frame['group'] == combs[0]].values
        r_vec1 = ro.FloatVector([f for f in vec1])
        vec2 = all_lnc_frame['density'][all_lnc_frame['group'] == combs[1]].values
        r_vec2 = ro.FloatVector([g for g in vec2])
        res = wilcoxpy(r_vec1, r_vec2, alternative="greater")
        pval = res.rx('p.value')[0][0]
        stat = res.rx('statistic')[0][0]
        test_dict[(combs[0], combs[1])] = {"W": stat,
                                           "p-value": pval}

    test_table = pd.DataFrame(test_dict).T
    cond = cor_file.split("/")[-1].split("-")[0]
    test_table.columns = ['W', 'p-value']
    test_table['condition'] = cond
    test_table.to_csv("stats.dir/%s-MRE_density-stats.tsv" % cond,
                      sep="\t", index_label="reference")
    R.assign("r.df", r_lnc_df)

    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''r.df$group <- as.factor(r.df$group)''')
    R('''p_des <- ggplot(r.df, aes(x=density, fill=group, colour=group)) + '''
      '''geom_density(alpha=0.2) + labs(x="MREs per kb", y="density") + '''
      '''theme(text=element_text(size=14, colour="black")) + theme_bw()''')
    R('''png('%s', height=480, width=480)''' % outfile)
    R('''print(p_des)''')
    R('''dev.off()''')


@P.cluster_runnable
def plotMreCounts(counts_file, outfile):
    '''
    plot output from countMREsOverLncs as histograms
    for both genes and lncRNAs
    '''
    
    mre_frame = pd.read_table(counts_file,
                              sep="\t",
                              header=None,
                              index_col=0)
    mre_frame.columns = ['MRE_counts']
    mre_frame.index.name = 'gene_id'
    mre_frame['biotype'] = ['' for cx in mre_frame.index]
    df_index = mre_frame.index.tolist()

    lncs = [cl for cl in df_index if re.search("LNC", cl)]
    genes = [cg for cg in df_index if re.search("ENS", cg)]

    lnc_counts = mre_frame.loc[lncs]
    gene_counts = mre_frame.loc[genes]
    lnc_counts['biotype'] = "lncRNA"
    gene_counts['biotype'] = "gene"

    cor_mres = gene_counts.append(lnc_counts)

    tot_val = len(cor_mres['MRE_counts'].values)
    chained = itertools.chain(lnc_counts['MRE_counts'].values,
                              gene_counts['MRE_counts'].values)
    max_val = max([s for s in chained])
    # if all values are zero, touch a sentinel file and
    # break out of function
    if max_val == 0:
        P.touch(outfile)
        return 0
    else:
        pass

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    lnc_vals = lnc_counts['MRE_counts'].values
    binwidth = int(max_val/float(tot_val/5.0))
    ax1.grid(True)

    try:
        ax1.hist(lnc_vals,
                 facecolor="blue",
                 label="lncRNA MRE counts",
                 bins=range(0, max_val + binwidth, binwidth))
    except ValueError:
        ax1.hist(lnc_vals,
                 facecolor="blue",
                 label="lncRNA MRE counts",
                 bins=range(0, max_val + binwidth))
    ax1.legend()

    ax2 = fig.add_subplot(212)
    gene_vals = gene_counts['MRE_counts'].values
    try:
        ax2.hist(gene_vals,
                 facecolor="red",
                 label="gene MRE counts",
                 bins=range(0, max_val + binwidth, binwidth))
    except ValueError:
        ax2.hist(gene_vals,
                 facecolor="red",
                 label="gene MRE counts",
                 bins=range(0, max_val + binwidth))
    ax2.grid(True)
    ax2.legend()
    fig.savefig(outfile)


def plotViolinCounts(infile, outfile):
    '''
    Generate ggplot violin plots of MRE count distributions
    for genes and lncRNAs
    '''

    # use R code for now - need to work this out in matplotlib

    mre_df = pd.read_table(infile, sep="\t", header=None, index_col=0)
    mre_df.columns = ['MRE_counts']
    mre_df.index.name = "gene_id"
    mre_df['biotype'] = ['' for px in mre_df.index]
    idx = mre_df.index.tolist()

    lncs = [pl for pl in idx if re.search("LNC", pl)]
    genes = [pg for pg in idx if re.search("ENS", pg)]

    lnc_counts = mre_df.loc[lncs]
    gene_counts = mre_df.loc[genes]
    lnc_counts['biotype'] = "lncRNA"
    gene_counts['biotype'] = "gene"

    all_counts = gene_counts.append(lnc_counts)
    all_idx = all_counts.index.tolist()

    r_df = pandas2ri.py2ri_pandasdataframe(all_counts)
    R.assign("r.df", r_df)

    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''p_g <- ggplot(r.df, aes(y=MRE_counts, x=as.factor(biotype),'''
      '''colour=as.factor(biotype), fill=as.factor(biotype))) + '''
      '''geom_violin() + labs(x="Biotype", y="MRE counts") + '''
      '''guides(colour=F, fill=F)''')
    R('''png("%(outfile)s", height=540, width=540)''' % locals())
    R('''print(p_g)''')
    R('''dev.off()''')


@P.cluster_runnable
def plotSharedCounts(cor_file,
                     prox_file,
                     random_file,
                     anti_file,
                     outfile):
    ''' density plots of shared miRNAs between gene:lncRNA pairs'''

    if cor_file.endswith("gz"):
        cor_comp = "gzip"
    else:
        cor_comp = None

    if prox_file.endswith("gz"):
        prox_comp = "gzip"
    else:
        prox_comp = None

    if random_file.endswith("gz"):
        rand_comp = "gzip"
    else:
        rand_comp = None

    if anti_file.endswith("gz"):
        anti_comp = "gzip"
    else:
        anti_comp = None

    cor_df = pd.read_table(cor_file, sep="\t", index_col=None,
                           header=0, compression=cor_comp,
                           comment='#')
    cdenom = cor_df['total_shared']/cor_df['unshared']
    cshare = cor_df['total_shared']/cdenom
    cor_df['prop_shared'] = cshare
    cor_df['ratio_shared'] = cor_df['total_shared']/cor_df['unshared']

    prox_df = pd.read_table(prox_file, sep="\t", index_col=None,
                            header=0, compression=prox_comp,
                            comment='#')
    pdenom = prox_df['total_shared']/prox_df['unshared']
    pshare = prox_df['total_shared']/pdenom
    prox_df['prop_shared'] = pshare
    prox_df['ratio_shared'] = prox_df['total_shared']/prox_df['unshared']

    rand_df = pd.read_table(random_file, sep="\t", index_col=None,
                            header=0, compression=rand_comp,
                            comment='#')
    rdenom = rand_df['total_shared']/rand_df['unshared']
    rshare = rand_df['total_shared']/rdenom
    rand_df['prop_shared'] = rshare
    rand_df['ratio_shared'] = rand_df['total_shared']/rand_df['unshared']

    anti_df = pd.read_table(anti_file, sep="\t", index_col=None,
                            header=0, compression=anti_comp,
                            comment='#')
    adenom = anti_df['total_shared']/anti_df['unshared']
    ashare = anti_df['total_shared']/adenom
    anti_df['prop_shared'] = ashare
    anti_df['ratio_shared'] = anti_df['total_shared']/anti_df['unshared']

    cor_df['cat'] = "correlated"
    prox_df['cat'] = "proximal"
    rand_df['cat'] = "random"
    anti_df['cat'] = "anticorrelated"

    all_shared = cor_df.append(rand_df)
    all_shared = all_shared.append(prox_df)
    all_shared = all_shared.append(anti_df)
    # need to re-index data frame after append to prevent duplicate indices
    new = [x for x, y in enumerate(all_shared.index)]
    all_shared.index = new

    # formally test shared miRNAs between gene sets
    wilcoxpy = R['wilcox.test']
    test_dict = {}
    for combs in itertools.combinations(set(all_shared['cat']), r=2):
        vec1 = all_shared['total_shared'][all_shared['cat'] == combs[0]]
        r_vec1 = ro.FloatVector([f for f in vec1.values])
        vec2 = all_shared['total_shared'][all_shared['cat'] == combs[1]]
        r_vec2 = ro.FloatVector([g for g in vec2.values])
        res = wilcoxpy(r_vec1, r_vec2, alternative="greater")
        pval = res.rx('p.value')[0][0]
        stat = res.rx('statistic')[0][0]
        test_dict[(combs[0], combs[1])] = {"W": stat,
                                           "p-value": pval}

    test_table = pd.DataFrame(test_dict).T
    cond = cor_file.split("/")[-1].split("-")[0]
    test_table['condition'] = cond
    test_table.to_csv("stats.dir/%s-MRE_shared-stats.tsv" % cond,
                      sep="\t", index_label="reference")

    r_share = pandas2ri.py2ri_pandasdataframe(all_shared)
    R.assign("shared.df", r_share)

    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''shared.df$cat <- as.factor(shared.df$cat)''')
    R('''p_share <- ggplot(shared.df, aes(x=total_shared, '''
      '''fill=cat, colour=cat)) + geom_density(alpha=0.2) + '''
      '''labs(x="Total number of shared miRNAs", y="density") + '''
      '''theme(text=element_text(size=14, colour="black")) + theme_bw()''')
    R('''png("%s", height=480, width=480)''' % outfile)
    R('''print(p_share)''')
    R('''dev.off()''')


def getMREs(mre_file, pairs_gtf):
    '''
    Get MREs for all lncRNAs and genes
    '''

    trans_gene_dict = {}
    catalog = {}
    # log all miRNAs for each gene and lncRNA
    with IOTools.openFile(pairs_gtf, "rb") as gfile:
        pairs_it = GTF.iterator(gfile)
        for it in pairs_it:
            trans_gene_dict[it.transcript_id] = it.gene_id
            catalog[it.gene_id] = set()

    with IOTools.openFile(mre_file, "rb") as mfile:
        mre_it = GTF.iterator_filtered(GTF.iterator(mfile))
        for mre in mre_it:
            target = mre.asDict()['target']
            mirna = mre.asDict()['miRNA']
            gene = trans_gene_dict[target]
            catalog[gene].add(mirna)

    return catalog


@P.cluster_runnable
def shareMREs(mre_file, pairs_gtf, correlations, outfile):
    '''
    Find the shared MREs between highly correlated
    lncRNA:gene pairs.
    Output:
    * list of shared miRNAs for each pair
    * number and proportion of shared MREs between pairs
    '''

    catalog = getMREs(mre_file, pairs_gtf)

    if correlations.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    cor_df = pd.read_table(correlations,
                           sep="\t",
                           compression=comp,
                           header=0)

    shared_dict = {}
    for idx in cor_df.index:
        share = {}
        gene = cor_df.loc[idx]['gene_id']
        lnc = cor_df.loc[idx]['lncRNA_id']
        gmirs = catalog[gene]
        lmirs = catalog[lnc]
        shared = gmirs.intersection(lmirs)
        not_shared = gmirs.difference(lmirs)
        try:
            lnc_prop = len(shared)/float(len(lmirs))
        except ZeroDivisionError:
            lnc_prop = 0.0
        try:
            gene_prop = len(shared)/float(len(gmirs))
        except ZeroDivisionError:
            gene_prop = 0.0

        share['gene_id'] = gene
        share['lncRNA_id'] = lnc
        share['unshared'] = len(not_shared)
        share['total_shared'] = len(shared)
        share['lncRNA_shared_proportion'] = lnc_prop
        share['gene_shared_proportion'] = gene_prop
        share['shared_miRNAs'] = ",".join([x for x in shared])

        shared_dict[idx] = share

    out_frame = pd.DataFrame(shared_dict).T
    
    # get gene symbol ids for genes
    mg = mygene.MyGeneInfo()
    try:
        q_symbol = mg.querymany(out_frame['gene_id'].tolist(),
                                scopes="ensemblgene",
                                species="mouse",
                                fields="symbol",
                                returnall=True)['out']
    except AssertionError:
        gene_set = [gsx for gsx in set(out_frame['gene_id'].tolist())]
        q_symbol = mg.querymany(gene_set,
                                scopes="ensemblgene",
                                species="mouse",
                                fields="symbol",
                                returnall=True)['out']

    q_df = pd.DataFrame(q_symbol)
    try:
        q_df.drop(['_id', 'notfound'], inplace=True, axis=1)
    except ValueError:
        pass

    outdf = pd.merge(left=out_frame, right=q_df, how='inner',
                     left_on='gene_id', right_on='query')
    outdf.drop(['query'], inplace=True, axis=1)

    outdf.to_csv(outfile, sep="\t", index=None)


@P.cluster_runnable
def countSharedMREs(mre_file,
                    pairs_gtf,
                    shared_file,
                    outfile):
    '''
    Count the number of elements in each gene and lncRNA for
    which a targetting miRNA is shared.
    '''
    if shared_file.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    shared_df = pd.read_table(shared_file,
                              sep="\t",
                              header=0,
                              index_col=None,
                              compression=comp)

    # 'shared_miRNAs' are a single string of comma-separated ids
    # need to split these to store for later in an array/list
    shared_df = shared_df.fillna("")
    shared_mirs = [h.split(",") for h in shared_df['shared_miRNAs'].values]
    shared_df['shared_miRNAs'] = shared_mirs

    # make a dictionary mapping gene_ids onto transcript ids - mre.gtf
    # only contains transcript ids
    # catalog will be all of the miRNA ids mapping to a gene/lncRNA
    trans_gene_dict = {}
    catalog = {}
    mre_dict = {}

    gene_ids = shared_df['gene_id'].values
    lnc_ids = shared_df['lncRNA_id'].values

    with IOTools.openFile(pairs_gtf, "rb") as gfile:
        pairs_it = GTF.iterator(gfile)
        for it in pairs_it:
            trans_gene_dict[it.transcript_id] = it.gene_id
            catalog[it.gene_id] = set()

    with IOTools.openFile(mre_file, "rb") as mfile:
        mre_it = GTF.iterator_filtered(GTF.iterator(mfile), "MRE")
        for mre in mre_it:
            gene = trans_gene_dict[mre.asDict()['target']]
            mre_entry = {'miRNA': mre.asDict()['miRNA'],
                         'target': gene,
                         'seed_class': mre.asDict()['seed_class'],
                         'contig': mre.contig,
                         'start': mre.start,
                         'end': mre.end,
                         'strand': mre.strand}
            mre_dict[mre.gene_id] = mre_entry

    mre_df = pd.DataFrame(mre_dict).T

    # count the number of MREs for each miRNA in each gene/lncRNA
    target_dict = {}
    targets = set(mre_df['target'].values)
    for tar in targets:
        mir_dict = {}
        tar_df = mre_df[mre_df['target'] == tar]
        mirs = set(tar_df['miRNA'].values)
        for mi in mirs:
            mir_dict[mi] = len(tar_df[tar_df['miRNA'] == mi])
            target_dict[tar] = mir_dict

    # count sites for shared miRNAs
    shared_gene_counts = {}
    for gene in gene_ids:
        mres = {}
        shared_genes = shared_df[shared_df["gene_id"] == gene]
        shared_mirs = shared_genes['shared_miRNAs'].values[0]
        for mi in shared_mirs:
            if len(mi):
                count = target_dict[gene][mi]
                mres[mi] = count
            else:
                pass
        shared_gene_counts[gene] = mres

    shared_lnc_counts = {}
    for lnc in lnc_ids:
        mres = {}
        shared_lncs = shared_df[shared_df['lncRNA_id'] == lnc]
        shared_mirs = shared_lncs['shared_miRNAs'].values
        for mi in shared_mirs[0]:
            if len(mi):
                count = target_dict[lnc][mi]
                mres[mi] = count
            else:
                pass
        shared_lnc_counts[lnc] = mres

    # generate the final table of genes, lncs and shared miRNAs with counts
    shared_mres_dict = {}
    for idx in shared_df.index:
        gene = shared_df.iloc[idx]['gene_id']
        lnc = shared_df.iloc[idx]['lncRNA_id']
        for mir in shared_df.iloc[idx]['shared_miRNAs']:
            try:
                gene_mirs = shared_gene_counts[gene][mir]
            except KeyError:
                gene_mirs = 0
            try:
                lnc_mirs = shared_lnc_counts[lnc][mir]
            except KeyError:
                lnc_mirs = 0
            mir_dict = {'miRNA': mir,
                        'gene_counts': gene_mirs,
                        'lncRNA_counts': lnc_mirs}
            shared_mres_dict[(gene, lnc)] = pd.Series(mir_dict)
    shared_mres_df = pd.DataFrame(shared_mres_dict).T

    shared_mres_df.to_csv(outfile, sep="\t")


@P.cluster_runnable
def correlateRandomPairs(pairs_file, expression, ref_gtf, outfile, seed):
    '''
    Cross-correlate random pairs of lncRNAs and protein-coding
    genes
    '''

    if pairs_file.split(".")[-1] == "gz":
        pair_comp = "gzip"
    else:
        pair_comp = None

    if expression.split(".")[-1] == "gz":
        expr_comp = "gzip"
    else:
        expr_comp = None

    expr_df = pd.read_table(expression,
                            compression=expr_comp,
                            sep="\t",
                            header=0,
                            index_col=0)

    pairs_df = pd.read_table(pairs_file,
                             compression=pair_comp,
                             sep="\t",
                             header=0,
                             index_col=0)
    lnc_ids = pairs_df['lncRNA_id']
    gene_ids = pairs_df.index

    all_lncs = [l for l in expr_df.index if re.search("LNC", l)]
    l_expr = expr_df.loc[all_lncs]
    l_expr.index.name = "lncRNA_id"

    all_genes = [g for g in expr_df.index if re.search("ENS", g)]
    g_expr = expr_df.loc[all_genes]

    # get lncRNA classifications from reference gtf
    # get transcript lengths for matching
    ofile = IOTools.openFile(ref_gtf, "rb")
    gene_it = GTF.transcript_iterator(GTF.iterator(ofile))
    class_dict = {}
    length_dict = {}
    for gene in gene_it:
        for trans in gene:
            class_dict[trans.gene_id] = {'class': trans.source,
                                         'exon': trans.asDict()['exon_status']}
            try:
                length_dict[trans.gene_id] += (trans.end - trans.start)
            except KeyError:
                length_dict[trans.gene_id] = (trans.end -  trans.start)

    ofile.close()

    # expression from pairs file
    pairs_lexpr = l_expr.loc[lnc_ids]

    # randomly sub sample genes and lncRNAs
    # match lncRNA expression to +- 0.5 of correlated lncRNAs
    # and length within 1kb
    random.seed(seed)
    lnc_idxs = set()
    while len(lnc_idxs) != len(lnc_ids):
        # randomly select a lncRNA from all expressed lncRNAs
        r_lnc = random.randint(0, len(all_lncs) - 1)
        r_lnc_name = expr_df.iloc[r_lnc].name

        # randomly select a matched lncRNA from the pairs file
        r_match = random.randint(0, len(lnc_ids) - 1)

        # check these are not the same lncRNA
        if r_lnc_name != pairs_lexpr.iloc[r_match].name:
            rclass = class_dict[r_lnc_name]['class']
            rexon = class_dict[r_lnc_name]['exon']
            lnc_len = length_dict[r_lnc_name]
            # select multi-exonic intergenic lncRNAs only
            if rclass == "intergenic" and rexon == "m":
                hi_xpr = np.mean(pairs_lexpr.iloc[r_match]) + 0.5
                lo_xpr = np.mean(pairs_lexpr.iloc[r_match]) - 0.5
                hi_len = lnc_len + 1000
                lo_len = lnc_len - 1000
                if hi_xpr < np.mean(expr_df.iloc[r_lnc]):
                    pass
                elif lo_xpr > np.mean(expr_df.iloc[r_lnc]):
                    pass
                else:
                    if lnc_len > hi_len:
                        pass
                    elif lnc_len < lo_len:
                        pass
                    else:
                        # only add random lnc if matched on expression
                        # lncRNA transcript length
                        # and lncRNA classification, but not ID
                        lnc_idxs.add(r_lnc)
            else:
                pass
        else:
            pass

    gene_idxs = set()
    while len(gene_idxs) != len(gene_ids):
        gene_idxs.add(random.randint(0, len(all_genes) - 1))

    # correlate random genes and lncRNAs
    rand_lncs = l_expr.iloc[[i for i in lnc_idxs]]
    rand_genes = g_expr.iloc[[q for q in gene_idxs]]
    r_lncs = rand_lncs.index
    r_genes = rand_genes.index

    rand_cor_df = pd.DataFrame(index=r_lncs,
                               columns=['gene_id', 'correlation'])
    for each in itertools.izip(r_lncs, r_genes):
        lval = rand_lncs.loc[each[0]].tolist()
        gval = rand_genes.loc[each[1]].tolist()
        rcor = TS.crossCorrelate(lval, gval, lag=0)
        rand_cor_df.loc[each[0]]['gene_id'] = each[1]
        rand_cor_df.loc[each[0]]['correlation'] = rcor[0]

    rand_cor_df['lncRNA_id'] = rand_cor_df.index
    rand_cor_df.index = rand_cor_df['gene_id']
    rand_cor_df.drop(['gene_id'], inplace=True, axis=1)
    rand_cor_df.to_csv(outfile, sep="\t", index_label="gene_id")


@P.cluster_runnable
def antiCorrelatePairs(pairs_file, expression, outfile, threshold):
    '''
    Get lncRNAs with paired protein-coding genes that 
    are anti-correlated in expression
    '''

    # need to restrict the number of transcripts/lncRNAs correlated.
    if pairs_file.split(".")[-1] == "gz":
        cor_comp = "gzip"
    else:
        cor_comp = None

    if expression.split(".")[-1] == "gz":
        expr_comp = "gzip"
    else:
        expr_comp = None

    pair_df = pd.read_table(pairs_file,
                            compression=cor_comp,
                            sep="\t",
                            header=0,
                            index_col=0)

    expr_df = pd.read_table(expression,
                            compression=expr_comp,
                            sep="\t",
                            header=0,
                            index_col=0)

    # select lncRNAs that are highly correlated
    # select genes that are anti-correlated with these lncRNAs
    lncs = set([l for l in pair_df['lncRNA_id'] if re.search("LNC", l)])
    genes = set([g for g in expr_df.index if re.search("ENS", g)])

    # get correlations between all lncs and genes,
    cor_frame = pd.DataFrame(index=lncs, columns=genes)
    cor_frame = cor_frame.fillna(0.0)

    pairs = itertools.product(lncs, genes)
    for each in pairs:
        lnc_val = expr_df.loc[each[0]].tolist()
        gene_val = expr_df.loc[each[1]].tolist()
        cor = TS.crossCorrelate(lnc_val, gene_val, lag=0)[0]
        cor_frame.loc[each[0], each[1]] = cor
    cor_frame = cor_frame.fillna(0.0)
    cor_frame.index.name = "lncRNA_id"
    unstack = cor_frame.unstack()
    cor_list = unstack.reset_index()
    cor_list.columns = ['gene_id', 'lncRNA_id', 'correlation']
    anticor_list = cor_list[cor_list['correlation'] <= threshold]
    anticor_list.index = anticor_list['gene_id']
    anticor_list.drop(['gene_id'], inplace=True, axis=1)

    anticor_list.to_csv(outfile, sep="\t", index_label="gene_id")
    

@P.cluster_runnable
def correlateProximalPairs(distances, pairs_file, expression, outfile):
    '''
    Get lncRNAs with most proximal protein-coding gene,
    calculate cross-correlation of expression.
    '''
    if pairs_file.split(".")[-1] == "gz":
        cor_comp = "gzip"
    else:
        cor_comp = None

    if distances.split(".")[-1] == "gz":
        dist_comp = "gzip"
    else:
        dist_comp = None

    if expression.split(".")[-1] == "gz":
        expr_comp = "gzip"
    else:
        expr_comp = None

    pair_df = pd.read_table(pairs_file,
                            compression=cor_comp,
                            sep="\t",
                            header=0,
                            index_col=0)

    dist_df = pd.read_table(distances,
                            compression=dist_comp,
                            sep="\t",
                            header=0,
                            index_col=0)

    expr_df = pd.read_table(expression,
                            compression=expr_comp,
                            sep="\t",
                            header=0,
                            index_col=0)
    
    lnc_dists = dist_df.loc[set(pair_df['lncRNA_id'].values)]
    lnc_expr = expr_df.loc[set(lnc_dists.index)]
    gene_expr = expr_df.loc[lnc_dists['closest_id']]

    gene_expr['gene_id'] = gene_expr.index
    gene_expr.drop_duplicates(subset='gene_id',
                              take_last=True,
                              inplace=True)
    gene_expr.drop(['gene_id'], inplace=True, axis=1)

    lncs = lnc_expr.index
    genes = gene_expr.index

    # get correlations between all lncs and genes,
    # regardless of proximity - subset later
    cor_frame = pd.DataFrame(index=lncs, columns=genes)
    cor_frame = cor_frame.fillna(0.0)

    pairs = itertools.product(lncs, genes)
    for each in pairs:
        lnc_val = lnc_expr.loc[each[0]].tolist()
        gene_val = gene_expr.loc[each[1]].tolist()
        cor = TS.crossCorrelate(lnc_val, gene_val, lag=0)
        cor_frame.loc[each[0]][each[1]] = cor
    cor_frame = cor_frame.fillna(0.0)
    cor_frame.index.name = "lncRNA_id"
    unstack = cor_frame.unstack()
    cor_list = unstack.reset_index()
    cor_list.columns = ['gene_id', 'lncRNA_id', 'correlation']
    cor_list.index = cor_list['lncRNA_id']
    cor_list.drop(['lncRNA_id'], inplace=True, axis=1)

    prox_cors = {}
    for idx in cor_list.index:
        cors = cor_list.loc[idx]
        cors.index = cors['gene_id']
        prox_gene = lnc_dists.loc[idx]['closest_id']
        prox_cors[idx] = {'gene_id': prox_gene,
                          'correlation': cors.loc[prox_gene]['correlation']}
    prox_cor_df = pd.DataFrame(prox_cors).T
    prox_cor_df['lncRNA_id'] = prox_cor_df.index
    prox_cor_df.index = prox_cor_df['gene_id']
    prox_cor_df.drop(['gene_id'], inplace=True, axis=1)

    prox_cor_df.to_csv(outfile, sep="\t", index_label="gene_id")


def tempCorr(x, y):
    '''
    Temporal correlation of two time series
    of the form:
    ux, uy = mean of x and y

    Corr(x, y) = sum((x-ux)(y-uy))/(sqrt(var(x)) * sqrt(var(x)))
    '''

    sum_prod = []
    sum_xsq = []
    sum_ysq = []

    for i in range(len(x) - 1):
        xi = float(x[i+1]) - float(x[i])
        yi = float(y[i+1]) - float(y[i])
        prod = xi * yi
        sum_prod.append(prod)
        sq_x = xi**2
        sq_y = yi**2

        sum_xsq.append(sq_x)
        sum_ysq.append(sq_y)

    nume = sum(sum_prod)
    denom = float(math.sqrt(sum(sum_xsq)) * math.sqrt(sum(sum_ysq)))

    if denom != 0:
        return nume/denom
    else:
        return 0


def crossCorrelate(t, s, lag=0):
    '''
    Calculate the cross-correlation of two timeseries, s and t.
    Return the normalized correlation value at lag=n.
    Uses numpy.correlate; default is to return lag=0.
    TODO: return multiple lags?
    '''

    t_mean = np.mean(t)
    s_mean = np.mean(s)
    t_std = np.std(t)
    s_std = np.std(s)
    len_t = len(t)

    t_norm = [((x - t_mean)/(t_std * len_t)) for x in t]
    s_norm = [((y - s_mean)/s_std) for y in s]

    if lag == 0:
        xcorr = np.correlate(t_norm, s_norm)
    elif lag != 0:
        xcorr = np.correlate(t_norm, s_norm, mode=2)[len_t - 1 + lag]

    return xcorr


def correlateEigengenes(ref_eigens,
                        lnc_eigens,
                        correlation,
                        lag=0):
    '''
    Correlate two sets of eigenenes.
    Specify correlation types.  Current options are temporal
    and cross-correlation.
    '''

    if ref_eigens.split(".")[-1] == "gz":
        ref_comp = "gzip"
    else:
        ref_comp = None

    if lnc_eigens.split(".")[-1] == "gz":
        lnc_comp = "gzip"
    else:
        lnc_comp = None

    ref_df = pd.read_table(ref_eigens, sep="\t",
                           header=0, index_col=0,
                           compression=ref_comp)

    lnc_df = pd.read_table(lnc_eigens, sep="\t",
                           header=0, index_col=0,
                           compression=lnc_comp)

    corr_frame = correlateLncRNAs(lnc_df,
                                  ref_df,
                                  correlation,
                                  lag)

    return corr_frame


def correlateLncRNAs(lnc_frame, gene_frame, correlation, lag=0):
    '''
    Use temporal correlation to correlate lncRNA time series
    expression profiles with input expression profiles.
    '''

    lnc_id = lnc_frame.index
    gene_id = gene_frame.index
    cor_frame = pd.DataFrame(index=lnc_id, columns=gene_id)
    cor_frame = cor_frame.fillna(0.0)

    if correlation == "temporal":
        for x in itertools.product(lnc_id, gene_id):
            lnc_vals = lnc_frame.loc[x[0]].tolist()
            gene_vals = gene_frame.loc[x[1]].tolist()
            corr = tempCorr(lnc_vals, gene_vals)
            cor_frame[x[1]][x[0]] = corr

    elif correlation == "cross-correlation":
        for x in itertools.product(lnc_id, gene_id):
            lnc_vals = lnc_frame.loc[x[0]].tolist()
            gene_vals = gene_frame.loc[x[1]].tolist()
            corr = crossCorrelate(lnc_vals, gene_vals, lag)
            cor_frame[x[1]][x[0]] = corr

    return cor_frame


def filterCorrelations(infile):
    '''
    output list of gene1:gene2:value
    '''

    if infile.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    cor_frame = pd.read_table(infile, sep="\t", header=0,
                              index_col=0, compression=comp)

    cor_list = cor_frame.unstack()
    cor_list = cor_list.reset_index()
    cor_list.columns = ['lncRNA_id', 'gene_id', 'value']
    cor_list.index = cor_list['gene_id']
    cor_list.drop(['gene_id'], inplace=True, axis=1)

    return cor_list


def correlateLncRNAsWithClusterGenes(cor_file,
                                     expr_file,
                                     set1_clusters,
                                     set2_clusters,
                                     correlation,
                                     lag=0):
    '''
    Correlate lncRNAs with genes within clusters for which
    the cluster eigengene is correlated with the lncRNA

    cor_file = correlation of set1 and set2 eigengenes

    expr_file = expression data of genes and lncRNAs in all clusters

    set1_clusters = cluster labels for expr_file gene_ids
    set2_clusters = cluster labels for lnc_expr lncRNA_ids

    correlation = correlation measure; cross-correlation or temporal
    threshold = correlation threshold (positive only)
    lag = lag to report, for cross-correlation only
    '''

    # handle gzip compressed input files

    if cor_file.split(".")[-1] == "gz":
        cor_comp = "gzip"
    else:
        cor_comp = None
    if expr_file.split(".")[-1] == "gz":
        expr_comp = "gzip"
    else:
        expr_comp = None
    if set1_clusters.split(".")[-1] == "gz":
        set1_comp = "gzip"
    else:
        set1_comp = None
    if set2_clusters.split(".")[-1] == "gz":
        set2_comp = "gzip"
    else:
        set2_comp = None

    cor_df = pd.read_table(cor_file,
                           sep="\t",
                           index_col=None,
                           header=0,
                           compression=cor_comp)
    cor_df.columns = ['set1', 'set2', 'value']

    expr_df = pd.read_table(expr_file,
                            sep="\t",
                            index_col=0,
                            header=0,
                            compression=expr_comp)

    # cluster ids for expr_file
    set1_df = pd.read_table(set1_clusters,
                            sep="\t",
                            index_col=0,
                            header=0,
                            compression=set1_comp)
    # cluster ids for lnc_expr
    set2_df = pd.read_table(set2_clusters,
                            sep="\t",
                            index_col=0,
                            header=0,
                            compression=set2_comp)

    set1_df.columns = ['gene_id', 'cluster']
    set2_df.columns = ['gene_id', 'cluster']

    set1_df.index = set1_df['gene_id']
    set2_df.index = set2_df['gene_id']

    set1_df.drop(['gene_id'], inplace=True, axis=1)
    set2_df.drop(['gene_id'], inplace=True, axis=1)

    corr_dict = {}
    for x in cor_df.index:
        set1_ids, set2_ids = cor_df.loc[x]['set1'], cor_df.loc[x]['set2']
        gene_ids = set1_df[set1_df['cluster'] == set1_ids]
        lnc_ids = set2_df[set2_df['cluster'] == set2_ids]
        lnc_vals = expr_df.loc[lnc_ids.index.tolist()]
        gene_vals = expr_df.loc[gene_ids.index.tolist()]

        # select lncRNAs and genes in correlated cluster eigengenes
        # output gene:lncRNA:correlation:gene_cluster:lncRNA_cluster

        E.info("correlations for genes in cluster %s "
               "and lncRNAs in cluster %s" % (set1_ids, set2_ids))

        cluster_cor = correlateLncRNAs(lnc_vals,
                                       gene_vals,
                                       correlation,
                                       lag)

        cluster_cor['lncRNA_id'] = cluster_cor.index
        cor_list = pd.melt(cluster_cor, id_vars='lncRNA_id')
        cor_list['gene_cluster'] = set1_ids
        cor_list['lncRNA_cluster'] = set2_ids
        cor_list.index = cor_list['gene_id']
        cor_list.drop(['gene_id'], inplace=True, axis=1)
        corr_dict[(set1_ids, set2_ids)] = cor_list

    clusters = corr_dict.keys()
    cluster1 = clusters[0]
    clusters.remove(cluster1)

    results_frame = corr_dict[cluster1]

    for clust in clusters:
        results_frame = results_frame.append(corr_dict[clust])

    return results_frame


@P.cluster_runnable
def compareFovsGC(infile, fo_gc, image_dir):
    '''
    Compare results from time point differential expression analysis
    to Fo -> GC differential analysis results.
    '''

    name_list = infile.split("/")[-1].split("_")
    p_name = name_list[0] + "_" + name_list[2]
    p_name = p_name.rstrip("-time.tsv")
    df = pd.read_table(infile,
                       sep="\t",
                       header=0,
                       index_col=0)

    # select differentially expressed genes with p <= 0.01
    # intersect gene_ids and subset these for plotting

    df = df[df['padj'] <= 0.01]
    fo_gc = fo_gc[fo_gc['padj'] <= 0.01]

    agree = []
    overlap = set(df.index).intersection(set(fo_gc.index))
    for x in overlap:
        val1 = df.loc[x]['log2FoldChange']
        val2 = fo_gc.loc[x]['log2FoldChange']
        if (val1 > 0) and (val2 > 0):
            agree.append(x)
        elif (val1 < 0) and (val2 < 0):
            agree.append(x)
        else:
            pass

    # merge dfs on gene_id, keep log2 fold changes only

    merged = pd.merge(left=df.loc[agree],
                      right=fo_gc.loc[agree],
                      how='inner',
                      left_index=True,
                      right_index=True)
    merged = merged[['log2FoldChange_x', 'log2FoldChange_y']]
    columns = ['%s_l2fc' % p_name, 'fo_gc_l2fc']
    merged.columns = columns
    merged = merged.fillna(0.0)

    ggPlotRScatter(merged, p_name, image_dir)


def ggPlotRScatter(df, p_name, image_dir):
    '''
    Generate scatter plots of Fo->GC vs time points
    for intersecting differentially expressed genes.
    Colour by |difference in log2 fold change|.
    '''

    df['diff'] = abs(df['%s_l2fc' % p_name] - df['fo_gc_l2fc'])
    df = df.fillna(0.0)

    # set up ggplot components in R
    pandas2ri.activate()
    R.assign('df', df)
    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''p_base <- ggplot(df, aes(x=df[,1], '''
      '''y=fo_gc_l2fc, colour=diff))''')
    R('''geom <- geom_point(alpha=0.5)''')
    R('''coloring <- scale_color_gradient(low='blue', high='red')''')
    R('''labels <- labs(x="%(p_name)s log2 fold change", '''
      '''y="Fo->GC log2 fold change",'''
      '''title="log2 fold change correlation between\n'''
      ''' %(p_name)s and Fo->GC differential expression '''
      '''analyses")''' % locals())
    R('''xlimits <- xlim(-15, 15)''')
    R('''ylimits <- ylim(-15, 15)''')

    # put it all together

    R('''p_tot = p_base + geom + coloring + labels + xlimits + ylimits''')

    # save to image directoy
    # need to switch on x11 plotting device for ggsave to work
    # cannot use ggsave on cluster, revert to png(), do not turn on x11
    # R.x11()
    R('''png("%(image_dir)s/%(p_name)s-vsFo_GC.png")''' % locals())
    R('''print(p_tot)''')
    R('''dev.off()''')


@P.cluster_runnable
def coreOverlapFoVsGC(infile, fo_gc):
    '''
    Take a list of differentially expressed genes
    across each condition and time point, intersect with
    FovsGC
    '''

    name_list = infile.split("/")[-1].split("_")
    out_name = name_list[0] + "_" + name_list[2]
    out_name = out_name.rstrip("-time.tsv")
    out_name = out_name + "_Fo-GC-intersect.tsv"

    df = pd.read_table(infile,
                       sep="\t",
                       header=0,
                       index_col=0)

    # select differentially expressed genes with p <= 0.01
    # intersect gene_ids and subset these for plotting

    genes = df[df['padj'] <= 0.01].index.tolist()
    fo_gc_genes = fo_gc[fo_gc['padj'] <= 0.01].index.tolist()
    agree = []
    disagree = []
    overlap = set(genes).intersection(set(fo_gc_genes))
    for x in overlap:
        val1 = df.loc[x]['log2FoldChange']
        val2 = fo_gc.loc[x]['log2FoldChange']
        if (val1 > 0) and (val2 > 0):
            agree.append(x)
        elif (val1 < 0) and (val2 < 0):
            agree.append(x)
        else:
            disagree.append(x)

    # output list to file
    with IOTools.openFile("FovsGC_compare.dir/%s" % out_name, "w") as outfile:
        for gene in agree:
            outfile.write("%s\n" % gene)


@P.cluster_runnable
def plotTimepointIntersection(infiles, outfile):
    '''
    Plot Venn diagram of intersection of gene lists
    '''

    inter_dict = {}
    for fle in infiles:
        header = fle.split("/")[-1].split("-")[0]
        in_df = pd.read_table(fle, sep="\t", header=0, index_col=0)
        genes = in_df.index.tolist()
        inter_dict[header] = genes
        biotype = fle.split("/")[-1].split("-")[-1].split("_")[-1]
        biotype = biotype.split(".")[0]
    out_dir = "/".join(outfile.split("/")[:-1])

    TS.drawVennDiagram(inter_dict, biotype, out_dir)


@P.cluster_runnable
def getTimepointIntersections(infiles, n_times, outfile):
    '''
    Take first n timepoints and intersect for each in vitro
    activation condition.
    '''

    # get gene lists from files
    file_dictionary = {}
    for infile in infiles:
        gene_list = []
        with IOTools.openFile(infile, "rb") as gene_file:
            gene_list = gene_file.read().split("\n")
        gene_list.remove('')
        file_dictionary[infile] = gene_list
    # get intersections across time points
    time_point_dict = {}
    for tme in n_times:
        tpoints = [t for t in file_dictionary.keys() if re.search(str(tme), t)]
        time_set = set(file_dictionary[tpoints[0]])
        for i in range(1, len(tpoints)):
            gene_list = file_dictionary[tpoints[i]]
            time_set = time_set.intersection(gene_list)
        time_point_dict[str(tme)] = time_set
    # intersect all time points
    core_set = set(time_point_dict[str(n_times[0])])
    for j in range(1, len(time_point_dict.keys())):
        core_set = core_set.intersection(time_point_dict[str(n_times[j])])
    core_list = list(core_set)

    core_genes = [ge for ge in list(core_list) if re.search("EN", ge)]
    core_lncs = [lc for lc in list(core_list) if re.search("LNC", lc)]
    mg = mygene.MyGeneInfo()
    out_core = mg.querymany(core_genes,
                            scopes="ensemblgene",
                            fields="symbol",
                            returnall=True)['out']
    out_df = pd.DataFrame(out_core)
    out_df.drop(['notfound'], inplace=True, axis=1)
    out_df.index = out_df['query']
    out_df.drop_duplicates(subset='query', take_last=True, inplace=True)
    out_df.drop(['query'], inplace=True, axis=1)
    out_df.to_csv(outfile,
                  sep="\t",
                  index_label="gene_id")
    condition = outfile.split("-")[0]
    lnc_out = "%s-timepoint_intersection_lncRNAs.tsv" % condition
    with IOTools.openFile(lnc_out, "w") as lnc_file:
        lnc_file.write("lncRNA_id")
        for lncrna in core_lncs:
            lnc_file.write("%s\n" % lncrna)


@P.cluster_runnable
def getCoreGenes(infiles, n_times, outfile):
    '''
    Get files of gene lists that intersect with Fo-GC gene list
    and intersect across conditions for each time point, then
    across first n_times.
    '''

    # get gene lists from files
    file_dictionary = {}
    for infile in infiles:
        gene_list = []
        with IOTools.openFile(infile, "rb") as gene_file:
            gene_list = gene_file.read().split("\n")
        gene_list.remove('')
        file_dictionary[infile] = gene_list

    # get intersections across conditions at each time point
    time_point_dict = {}
    for tme in n_times:
        tpoints = [t for t in file_dictionary.keys() if re.search(str(tme), t)]
        time_set = set(file_dictionary[tpoints[0]])
        for i in range(1, len(tpoints)):
            gene_list = file_dictionary[tpoints[i]]
            time_set = time_set.intersection(gene_list)
        time_point_dict[str(tme)] = time_set

    # intersect all time points
    core_set = set(time_point_dict[str(n_times[0])])
    for j in range(1, len(time_point_dict.keys())):
        core_set = core_set.intersection(time_point_dict[str(n_times[j])])
    core_list = list(core_set)

    core_genes = [ge for ge in list(core_list) if re.search("EN", ge)]
    core_lncs = [lc for lc in list(core_list) if re.search("LNC", lc)]
    mg = mygene.MyGeneInfo()
    out_core = mg.querymany(core_genes,
                            scopes="ensemblgene",
                            fields="symbol",
                            returnall=True)['out']
    out_df = pd.DataFrame(out_core)
    out_df.drop(['notfound'], inplace=True, axis=1)
    out_df.index = out_df['query']
    out_df.drop_duplicates(subset='query', take_last=True, inplace=True)
    out_df.drop(['query'], inplace=True, axis=1)
    out_df.to_csv(outfile,
                  sep="\t",
                  index_label="gene_id")

    lnc_out = "%s-%s-core_lncRNAs.tsv" % (outfile.split("-")[0],
                                          outfile.split("-")[1])
    with IOTools.openFile(lnc_out, "w") as lnc_file:
        for lncrna in core_lncs:
            lnc_file.write("%s\n" % lncrna)


@P.cluster_runnable
def getConditionGenes(list_of_files, reference, outfile):
    '''
    Get time point intersection of genes and lncRNAs
    specific to each condition for n time points.
    '''

    inter_dict = {}
    # parse gene set per condition from files
    for fle in list_of_files:
        header = fle.split("/")[-1].split("-")[0]
        in_df = pd.read_table(fle, sep="\t", header=0, index_col=0)
        genes = in_df.index.tolist()
        inter_dict[header] = set(genes)

    # use difference_update iteratively to get difference between
    # reference and all other conditions
    spec_set = set(inter_dict[reference])

    for pair in itertools.product(inter_dict.keys(), inter_dict.keys()):
        # only get condition genes for specific condition
        if pair[0] == reference:
            if pair[0] == pair[1]:
                pass
            else:
                spec_set.difference_update(inter_dict[pair[1]])
        else:
            pass

    # detect genes - if so output gene symbols, otherwise skip and
    # write outfile
    if len([g for g in spec_set if re.search("ENS", g)]):
        geneset = [q for q in spec_set]
        mg = mygene.MyGeneInfo()
        out_core = mg.querymany(geneset,
                                scopes="ensemblgene",
                                fields="symbol",
                                returnall=True)['out']
        out_df = pd.DataFrame(out_core)
        try:
            out_df.drop(['notfound'], inplace=True, axis=1)
        except ValueError:
            pass
        out_df.index = out_df['query']
        out_df.drop_duplicates(subset='query', take_last=True, inplace=True)
        out_df.drop(['query'], inplace=True, axis=1)
        out_df.to_csv(outfile,
                      sep="\t",
                      index_label="gene_id")
    else:
        with IOTools.openFile(outfile, "w") as ofile:
            for obj in spec_set:
                ofile.write("%s\n" % obj)


@P.cluster_runnable
def correlateGeneLncRNA(gene_file, lnc_file, expression_file, outfile):
    '''
    Correlate gene and lncRNA expression for a given condition
    using temporal correlation.
    '''

    if gene_file.endswith("gz"):
        gene_comp = "gzip"
    else:
        gene_comp = None

    if lnc_file.endswith("gz"):
        lnc_comp = "gzip"
    else:
        lnc_comp = None

    if expression_file.endswith("gz"):
        expr_comp = "gzip"
    else:
        expr_comp = None

    gene_list = pd.read_table(gene_file, sep="\t", index_col=0, header=0,
                              compression=gene_comp)
    lnc_list = pd.read_table(lnc_file, sep="\t", index_col=0, header=None,
                             compression=lnc_comp)
    expr_df = pd.read_table(expression_file, sep="\t", header=0,
                            index_col=0, compression=expr_comp)

    gene_ids = gene_list.index.tolist()
    lnc_ids = lnc_list.index.tolist()

    corr_df = pd.DataFrame(index=gene_ids,
                           columns=lnc_ids)
    corr_df = corr_df.fillna(0.0)
    for x in itertools.product(gene_ids, lnc_ids):
        gene_val = expr_df.loc[x[0]].tolist()
        # difference in lncRNA expression thresholds may result in some lncRNA
        # not being included in final expression tables
        try:
            lnc_val = expr_df.loc[x[1]].tolist()
            corr = crossCorrelate(gene_val, lnc_val, lag=0)
            corr_df.loc[x[0], x[1]] = corr
        except KeyError:
            corr_df.loc[x[0], x[1]] = 0.0

    corr_df.to_csv(outfile, sep="\t",
                   index_label=None)


@P.cluster_runnable
def plotConditionHeatmap(infile, outfile):
    '''
    Plot heatmap of condition speicific cross-correlation
    between genes and lncRNAs
    '''

    if infile.endswith("gz"):
        comp = "gzip"
    else:
        comp = None
    cor_df = pd.read_table(infile, sep="\t",
                           header=0, index_col=0,
                           compression=comp)
    # remove duplicate entries
    genes = cor_df.index
    cor_df['genes'] = genes
    cor_df.drop_duplicates(['genes'], inplace=True,
                           take_last=True)
    cor_df.drop(['genes'], inplace=True, axis=1)
    plotHeatmap(cor_df, outfile)


def plotHeatmap(dataframe_object, outfile):
    '''
    plot heatmap with R::heatmap.2
    '''

    r_matrix = pandas2ri.py2ri_pandasdataframe(dataframe_object)
    R.assign("mat.vals", r_matrix)
    R('''suppressPackageStartupMessages(library(RColorBrewer))''')
    R('''suppressPackageStartupMessages(library(gplots))''')
    R('''hmcol <- colorRampPalette(brewer.pal(9, "PuOr"))(100)''')
    R('''png("%(outfile)s", units="in", res=300, '''
      '''height=5.3, width=5.3)''' % locals())
    R('''heatmap.2(as.matrix(mat.vals), trace="none", dendrogram="none", '''
      '''col=hmcol, labRow=F, labCol=F, margins=c(6,6))''')
    R('''dev.off''')


def correlate_lncsEigen(lncFile, eigenFile, correlation, lag=0):
    '''
    correlate lncRNA expression against module eigengene expression
    '''

    # average lncRNA expression across replicates and correlate with eigengene
    # expression

    if lncFile.split(".")[-1] == "gz":
        compression = "gzip"
    else:
        compression = None

    lncs_data = pd.read_table(lncFile, sep="\t",
                              header=0, index_col=0,
                              compression=compression)

    eigen_data = pd.read_table(eigenFile,
                               index_col=0,
                               header=0,
                               sep="\t",
                               compression=compression)

    corr_df = correlateLncRNAs(lncs_data,
                               eigen_data,
                               correlation,
                               lag)
    return corr_df


@P.cluster_runnable
def generateBackground(infiles, outfile):
    '''
    Output gene list of background genes for GO enrichment
    from a list of files.
    Requires gene id be in ensembl format
    '''

    # just select protein-coding genes for background, exclude
    # other genes/transcripts that may skew or distort
    # final enrichments

    gene_set = set()
    if type(infiles) != tuple:
        if infiles.endswith("gz"):
            comp = "gzip"
        else:
            comp = None
        df = pd.read_table(infiles, sep="\t", index_col=0,
                           header=0, compression=comp)
        genes = df.index.tolist()
        genes = [x for x in genes if re.search("ENS", x)]
        gene_set.update(genes)

    else:
        for fle in infiles:
            if fle.endswith("gz"):
                comp = "gzip"
            else:
                comp = None
            df = pd.read_table(fle, sep="\t", index_col=0,
                               header=0, compression=comp)
            genes = df.index.tolist()
            genes = [x for x in genes if re.search("ENS", x)]
            gene_set.update(genes)

    with IOTools.openFile(outfile, "w") as output:
        for x in gene_set:
            output.write("%s\n" % x)


@P.cluster_runnable
def goEnrichment(gene_set, bg_set, genome, db_ids, outfile, database):
    '''
    Perform GO enrichment on a single gene_set using the R/Bioconductor
    package GOseq.  Adjust for gene length only, not expression (later?)
    database choices = GO:BP, GO:MF, GO:CC, KEGG
    '''

    gene_df = pd.read_table(gene_set, sep="\t", header=0, index_col=0)
    gene_df['gene_id'] = gene_df.index
    pc_genes = [x for x in gene_df['gene_id'] if re.search("ENS", x)]
    gene_df = gene_df.loc[pc_genes]
    gene_df.drop_duplicates(subset='gene_id', inplace=True, take_last=True)

    bg_df = pd.read_table(bg_set, sep="\t", header=0, index_col=0)
    bg_df['gene_id'] = bg_df.index
    bg_df.drop_duplicates(subset='gene_id', inplace=True, take_last=True)

    geneset_r = ro.StrVector([x for x in gene_df.index.tolist()])
    bgset_r = ro.StrVector([l for l in bg_df.index.tolist()])

    # make a vector of integers for degs
    R('''core.de <- c(%s)''' % geneset_r.r_repr())
    R('''bg <- c(%s)''' % bgset_r.r_repr())
    R('''core.vector <- as.integer(bg%in%core.de)''')
    R('''names(core.vector) <- bg''')

    # generate pwf and perform GO enrichment
    R('''sink(file="sink_file.txt")''')
    R('''suppressPackageStartupMessages(library(goseq))''')
    R('''source("/ifs/projects/proj036/go_enrichment/GO2term.R")''')
    R('''pwf.core <- nullp(core.vector, "%(genome)s", '''
      ''' "%(db_ids)s", plot.fit=F)''' % locals())
    R('''go.core <- goseq(pwf.core, "%(genome)s", '''
      '''"%(db_ids)s", use_genes_without_cat=T, '''
      '''test.cats=c("%(database)s"))''' % locals())
    R('''GO.res <- GO2Term(go.core)''')
    R('''sink(file=NULL)''')

    go_df = pandas2ri.ri2pandas("GO.res")

    go_df.to_csv(outfile, sep="\t", index_label="category")


@P.cluster_runnable
def clusterGOEnrichment(cluster_file, genome, db_ids, label, out_dir):
    '''
    Perform GO enrichment on genes within each cluster - uses custom R script
    and Bioconductor package goseq
    database choices = GO:BP, GO:CC, GO:MF, KEGG
    '''

    R('''sink(file="sink_file.txt")''')
    R('''source("/ifs/projects/proj036/go_enrichment/GOseq_analysis.R")''')
    R('''go_enrichment(infile="%(cluster_file)s", species="%(genome)s", '''
      '''gene_id_db="%(db_ids)s", outfile_header="%(label)s", '''
      '''out_dir="%(out_dir)s")''' % locals())
    R('''sink(file=NULL)''')


@P.cluster_runnable
def topGO(go_file, expression_file, cluster_file, outfile):
    '''
    Calculate foldenrichments for each cluster gene ontolgoy
    enrichment file
    '''

    go_name = go_file.split("/")[-1].split("_")[1].split("GO")[0]
    go_frame = pd.read_table(go_file, sep="\t", index_col=None, header=0)
    go_frame.index = go_frame['category']

    express = pd.read_table(expression_file, sep="\t",
                            index_col=0, header=0)
    expr_genes = express.index.tolist()
    N = len(expr_genes)

    if cluster_file.endswith("gz"):
        comp = "gzip"
    else:
        comp = None
    clust_df = pd.read_table(cluster_file, sep="\t",
                             index_col=0, header=None,
                             compression=comp)
    clust_df.columns = ['gene_id', 'cluster']
    clusters = set(clust_df['cluster'])
    clust_dict = {}
    for x in clusters:
        clust_dict[x] = 0

    for gene in clust_df.index.values:
        clust_dict[clust_df.loc[gene]['cluster']] += 1

    # calculate fold enrichments for all GO terms
    n = int(clust_dict[go_name])
    fun = lambda x: round((x['numDEInCat']/float(n))/(x['numInCat']/float(N)),
                          3)
    enrich = go_frame.apply(fun, axis=1)
    go_frame['foldEnrichment'] = enrich

    go_frame.to_csv(outfile, sep="\t", index_label="category")


@P.cluster_runnable
def summariseGO(list_of_go, expression_file, cluster_file, outfile):
    '''
    Summarise gene ontology enrichments over a list of cluster
    '''

    go_dict = {}
    for go in list_of_go:
        go_name = go.split("/")[-1].split("_")[1].split("GO")[0]
        _df = pd.read_table(go, sep="\t", index_col=None, header=0)
        _df.index = _df['category']
        go_dict[go_name] = _df

    express = pd.read_table(expression_file, sep="\t",
                            index_col=0, header=0)
    expr_genes = express.index.tolist()
    N = len(expr_genes)

    if cluster_file.endswith("gz"):
        comp = "gzip"
    else:
        comp = None
    clust_df = pd.read_table(cluster_file, sep="\t",
                             index_col=0, header=None,
                             compression=comp)
    clust_df.columns = ['gene_id', 'cluster']
    clusters = set(clust_df['cluster'])
    clust_dict = {}
    for x in clusters:
        clust_dict[x] = 0

    for gene in clust_df.index.values:
        clust_dict[clust_df.loc[gene]['cluster']] += 1
    cluster_series = pd.Series(clust_dict)

    # calculate fold enrichments for all GO terms
    for clust in go_dict.keys():
        n = int(clust_dict[clust])
        func = lambda x: round((x['numDEInCat']/float(n))/(x['numInCat']/float(N)),
                               3)
        df_ = go_dict[clust]
        enrich = df_.apply(func, axis=1)
        df_['foldEnrichment'] = enrich
        go_dict[clust] = df_

    # summarise over all GO enrichments in clusters
    # take top ten enriched terms from each cluster
    top_dict = {}
    for each in go_dict.keys():
        go_df = go_dict[each]
        go_df.sort(columns='padjust', inplace=True, ascending=True)
        top_ = go_df.loc[go_df['ont'] == "BP"][0:10]
        top_['padjust'] = [float("%0.3g" % x) for x in top_['padjust'].values]
        top_dict[each] = top_

    # select top enrichment from each cluster for summary table
    one_dict = {}
    for name in top_dict.keys():
        one_df = top_dict[name]
        one_series = pd.Series(one_df.iloc[0])
        one_dict[name] = one_series

    GO_df = pd.DataFrame(one_dict).T
    GO_df['cluster'] = GO_df.index.tolist()
    GO_df.sort(inplace=True, columns='padjust', ascending=True)
    GO_df['padjust'] = [float("%0.3g" % x) for x in GO_df['padjust'].values]
    GO_df = pd.merge(GO_df,
                     pd.concat([GO_df['cluster'], cluster_series], axis=1),
                     left_index=True, right_index=True)
    GO_df = GO_df[['category', 'term', 'padjust', 'foldEnrichment', 'cluster']]
    GO_df.sort(columns='padjust', inplace=True, ascending=True)
    GO_df.drop(['cluster'], inplace=True, axis=1)
    GO_df.to_csv(outfile, sep="\t", index_label="cluster")


@P.cluster_runnable
def classifyLncRNA(lnc_list, lnc_gtf, lnc_class, direction, threshold):
    '''
    Classify a lncRNA set based on their correlation with protein-coding
    genes and positional classification.
    e.g. negatively correlated anti-sense lncRNAs or positively
    correlated intergenic lncRNAs
    Only select multi-exon lncRNAs
    '''

    if lnc_list.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    cor_df = pd.read_table(lnc_list, sep="\t",
                           header=0, index_col=None,
                           compression=comp)

    if direction == "positive":
        dir_cor = cor_df[cor_df['value'] >= 0]
    elif direction == "negative":
        dir_cor = cor_df[cor_df['value'] <= 0]
    else:
        raise AttributeError("Unrecognised correlation direction"
                             "Please supply positive or negative")
    try:
        cor_df.index = cor_df['set2']
    except KeyError:
        cor_df.index = cor_df['gene_id']

    lnc_iter = GTF.transcript_iterator(GTF.iterator(IOTools.openFile(lnc_gtf)))
    try:
        lnc_cors = set(dir_cor['set1'].tolist())
    except KeyError:
        try:
            lnc_cors = set(dir_cor['lncRNA'].tolist())
        except KeyError:
            lnc_cors = set(dir_cor['lncRNA_id'].tolist())
    try:
        gene_ids = set(cor_df['set2'].tolist())
    except KeyError:
        gene_ids = set(cor_df['gene_id'].tolist())

    lnc_frame = pd.DataFrame(index=gene_ids, columns=['value'])

    for trans in lnc_iter:
        lnc_id = list(set([x.gene_id for x in trans]))[0]
        lnc_source = list(set([x.source for x in trans]))[0]
        exon_status = set([x.asDict()['exon_status_locus'] for x in trans])
        lnc_status = list(exon_status)[0]

        # lncRNAs with exon_status_locs == m are multi-exon lncRNAs
        
        if lnc_source == lnc_class and lnc_status == 'm':
            if lnc_id in lnc_cors:
                try:
                    temp_frame = cor_df[cor_df['set1'] == lnc_id]
                except KeyError:
                    try:
                        temp_frame = cor_df[cor_df['lncRNA'] == lnc_id]
                    except KeyError:
                        temp_frame = cor_df[cor_df['lncRNA_id'] == lnc_id]
                lnc_frame = lnc_frame.append(temp_frame)
            else:
                pass
    
    # need to remove duplicate entries and NA's across genes and lncRNAs
    not_na = lnc_frame[np.isfinite(lnc_frame['value'])]
    try:
        not_na = not_na.drop_duplicates(subset=['set2', 'set1', 'value'],
                                        take_last=True)
    except KeyError:
        try:
            not_na = not_na.drop_duplicates(subset=['gene_id', 'lncRNA_id',
                                                    'value'],
                                            take_last=True)
        except KeyError:
            not_na = not_na.drop_duplicates(subset=['gene_id', 'lncRNA',
                                                    'value'],
                                            take_last=True)
    # drop cluster information if present
    try:
        not_na.drop(['gene_cluster', 'lncRNA_cluster'],
                    inplace=True, axis=1)
    except ValueError:
        pass

    # catch bug in number of columns due to column name differences
    if len(not_na.columns) > 3:
        not_na.drop(['gene_id'], inplace=True, axis=1)
    else:
        pass
    not_na.columns = ['gene_id', 'lncRNA_id', 'value']
    not_na.index = not_na['gene_id']
    not_na = not_na.drop(['gene_id'], axis=1)

    return not_na


def filter_correlation(infile,
                       threshold):
    '''
    filter lncRNAs with correlation below threshold
    '''
    
    cor_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    lncs_dict = {}

    modules = cor_frame.index.values.tolist()
    E.info("filtering correlations at threshold: %0.1f" % threshold)

    for lncRNA in cor_frame.columns.values:
        mod_dict = {}
        for mod in modules:
            if abs(cor_frame[lncRNA][mod]) >= threshold:
                mod_dict[mod] = cor_frame[lncRNA][mod]
            else:
                pass
        if len(mod_dict) > 0:
            lncs_dict[lncRNA] = mod_dict
        else:
            pass

    output_frame = pd.DataFrame(lncs_dict)
    output_frame['gene_id'] = output_frame.index
    cor_list = pd.melt(output_frame, id_vars='gene_id')
    flat_cor = cor_list[cor_list['value'] >= threshold]
    flat_cor.index = flat_cor['gene_id']
    flat_cor.columns = ['gene_id', 'lncRNA', 'value']
    flat_cor.drop(['gene_id'], inplace=True, axis=1)

    return flat_cor


def classify_lncrna(infile,
                    lnc_gtf,
                    summary_file,
                    out_gtf):
    '''
    classify lncRNAs based on the direction of their correlation with
    module eigengenes
    '''

    mod_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    modules = mod_frame.index.values.tolist()
    
    lnc_file = IOTools.openFile(lnc_gtf, "r")

    with open(out_gtf, "w") as gtf_file:
        for line in GTF.readFromFile(lnc_file):
            entry = GTF.Entry()
            if line.gene_id in mod_frame:
                entry.copy(line)
                for mod in modules:
                    mod_val = mod_frame[entry.gene_id][mod]
                    entry.addAttribute(key=mod,
                                       value=str(mod_val))
                    gtf_file.write("%s\n" % entry)
    
    sources = set()
    classification_dict = {}

    with open(out_gtf, "r") as new_gtf:
        for entry in GTF.readFromFile(new_gtf):
            sources.add(entry.source)
    
    for source in sources:
        classification_dict[source] = 0

    with open(out_gtf, "r") as openFile:
        for entry in GTF.readFromFile(openFile):
            classification_dict[entry.source] += 1
    (pd.Series(classification_dict)).to_csv(summary_file,
                                            sep="\t")
    
    pos_lncs = {}
    neg_lncs = {}
    
    with open(out_gtf, "r") as newfile:
        for entry in GTF.readFromFile(newfile):
            for mod in modules:
                try:
                    assert entry.asDict()[mod]
                    if float(entry.asDict()[mod]) < 0:
                        neg_lncs[entry.gene_id] = {'source': entry.source,
                                                   mod: 'negative'}

                    elif float(entry.asDict()[mod]) > 0:
                        pos_lncs[entry.gene_id] = {'source': entry.source,
                                                   mod: 'positive'}
                except(KeyError):
                    pass

    pos_series = pd.Series(pos_lncs)
    neg_series = pd.Series(neg_lncs)
    
    all_series = pos_series.append(neg_series)

    return all_series


def classCorrLncRNA(cluster_file,
                    gene_express,
                    lncRNA_express,
                    lnc_gtf,
                    threshold,
                    lncrna_class,
                    corr_direction):
    '''
    Classify lncRNAs based on correlation direction with protein
    coding gene expression
    '''
    
    def correlateDataFrames(df1, df2):
        '''
        Correlate 2 different dataframes, assuming matching column IDs
        but different indexes.  Default = Pearson correlation.
        '''

        # setup indices for dfs
        
        idx1 = df1.index.tolist()
        idx2 = df2.index.tolist()

        indexes = itertools.product(idx1, idx2)

        # create empty correlation dataframe
        correlate_frame = pd.DataFrame(index=idx1,
                                       columns=idx2)
        correlate_frame.fillna(0.0)

        for index in indexes:
            df1_series = df1.loc[index[0]].values[:-1]
            df2_series = df2.loc[index[1]].values

            # calculate Pearson correlation using numpy.corrcoef
            # np.corrcoef returns correlation matrix - need to index for
            # non-identity value

            corr = np.corrcoef(df1_series, df2_series)[0][1]
            correlate_frame.loc[index[0], index[1]] = corr

        return correlate_frame

    # for each cluster correlate lncRNAs against protein-coding genes
    # use Pearson correlation for now, but what about cross-correlation
    # or some sort of temporal correlation measure?

    clusters = pd.read_table(cluster_file,
                             sep="\t",
                             header=0,
                             index_col=0)
    clusters.columns = ['gene', 'cluster']

    cluster_cols = set(clusters['cluster'].tolist())

    gene_express = pd.read_table(gene_express,
                                 sep="\t",
                                 index_col=0,
                                 header=0)

    lnc_express = pd.read_table(lncRNA_express,
                                sep="\t",
                                index_col=0,
                                header=0)

    lnc_express['lncRNA'] = lnc_express.index.tolist()

    lnc_dict = {}
    lnc_index = GTF.iterator(IOTools.openFile(lnc_gtf))
    for lnc in lnc_index:
        lnc_dict[lnc.gene_id] = lnc.source
    
    class_frame = pd.DataFrame({'class': lnc_dict})
    intergenic_lncs = lnc_express[class_frame['class'] == lncrna_class]

    file_prefix = cluster_file.split("/")[1].split("-")[0]

    # just pull out the intergenic lncRNAs

    correlation_dict = {}
    all_correlations = {}

    for col in cluster_cols:

        # setup dataframe for cluster

        col_data = clusters[clusters['cluster'] == col]
        col_data.index = col_data['gene']
        col_data.drop(['gene'], inplace=True, axis=1)
        col_data_genes = col_data.index.tolist()
        
        # setup cluster-specific gene expression dataframe

        col_gene_express = gene_express.loc[col_data_genes]
        
        cor_frame = correlateDataFrames(intergenic_lncs,
                                        col_gene_express)

        # select lncRNAs on correlation

        correlated_lncs = []

        if corr_direction == "positive":
            for lncrna in cor_frame.index:
                if any([True for x in cor_frame.loc[lncrna] if x > threshold]):
                    correlated_lncs.append(lncrna)
                else:
                    pass
        elif corr_direction == "negative":
            for lncrna in cor_frame.index:
                lnc_correlations = cor_frame.loc[lncrna]
                logic_list = [True for x in lnc_correlations if x < -threshold]
                if any(logic_list):
                    correlated_lncs.append(lncrna)
                else:
                    pass

        lncs_cor_frame = cor_frame.loc[correlated_lncs]
        correlation_dict[col] = correlated_lncs

        # write out each correlation matrix to a separate file with cluster ID
        # write out list of correlated lncRNA IDs to file

        class_dir = "lncRNA_classification.dir"
        correlation_out = "%s/%s-%s-%s-correlations.tsv" % (class_dir,
                                                            file_prefix,
                                                            col,
                                                            corr_direction)
        lncs_cor_frame.to_csv(correlation_out,
                              sep="\t")
        correlated_lncs_out = "%s/%s-%s-%s-lncRNAs.tsv" % (class_dir,
                                                           file_prefix,
                                                           col,
                                                           corr_direction)
        lnc_out = IOTools.openFile(correlated_lncs_out, "w")
        for lnc in correlated_lncs:
            lnc_out.write("%s\n" % lnc)
        lnc_out.close()
        
        all_correlations[col] = cor_frame

    # iteratively merge each correlation frame onto the previous one
    # use lncRNAs IDs as index/keys

    total_frame = pd.concat(all_correlations.values(),
                            axis=1)

    return total_frame


@P.cluster_runnable
def list2GTF(list_of_ids, gtf_file, out_gtf):
    '''
    Turn a list of gene/lncRNA ids into a gtf file
    '''

    gtf_it = GTF.transcript_iterator(GTF.iterator(IOTools.openFile(gtf_file)))
    with IOTools.openFile(out_gtf, "w") as outfile:
        for trans in gtf_it:
            for exon in trans:
                if exon.gene_id in list_of_ids:
                    entry = GTF.Entry()
                    entry = entry.copy(exon)
                    outfile.write("%s\n" % entry)
                else:
                    pass


def correlationPairs(infile, threshold):
    '''
    Take a list of gene:lncRNA pairs with correlation coefficients.
    Output gene:lncRNA pairs with correlation >= threshold.
    Input table is in format gene_id:gene_cluster:lncRNA:lncRNA_cluster:value
    '''

    if infile.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    cors_df = pd.read_table(infile, sep="\t", header=0,
                            compression=comp, index_col=0)
    try:
        lnc_set = set(cors_df['lncRNA_id'])
    except KeyError:
        lnc_set = set(cors_df['lncRNA'])
    lncs_nn = {}
    idx = 0
    for lnc in lnc_set:
        try:
            l_df = cors_df[cors_df['lncRNA_id'] == lnc]
        except KeyError:
            l_df = cors_df[cors_df['lncRNA'] == lnc]
        if threshold >= 0:
            mgene = l_df[l_df['value'] >= threshold]
        elif threshold <= 0:
            mgene = l_df[l_df['value'] <= threshold]
        for xlnc in mgene.index:
            lncs_nn[str(idx)] = {'lncRNA_id': lnc,
                                 'gene_id': xlnc,
                                 'value': mgene.loc[xlnc]['value']}
            idx += 1

    out_frame = pd.DataFrame(lncs_nn).T
    out_frame.index = out_frame['gene_id']
    out_frame.drop(['gene_id'], axis=1, inplace=True)
    return out_frame


def maxCorrelationPairs(infile):
    '''
    Take a list of gene:lncRNA pairs with correlation coefficients.
    Output gene:lncRNA pairs with maximal cross-correlation.
    Input table is in format gene_id:gene_cluster:lncRNA:lncRNA_cluster:value
    '''

    if infile.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    cors_df = pd.read_table(infile, sep="\t", header=0,
                            compression=comp, index_col=0)
    try:
        lnc_set = set(cors_df['lncRNA_id'])
    except KeyError:
        lnc_set = set(cors_df['lncRNA'])
    lncs_nn = {}
    idx = 0
    for lnc in lnc_set:
        try:
            l_df = cors_df[cors_df['lncRNA_id'] == lnc]
        except KeyError:
            l_df = cors_df[cors_df['lncRNA'] == lnc]
        max_cor = max(l_df['value'])
        mgene = l_df[l_df['value'] == max_cor].index[0]
        lncs_nn[str(idx)] = {'lncRNA_id': lnc,
                             'gene_id': mgene,
                             'value': max_cor}
        idx += 1

    out_frame = pd.DataFrame(lncs_nn).T
    out_frame.index = out_frame['gene_id']
    out_frame.drop(['gene_id'], axis=1, inplace=True)
    return out_frame


def testGOCatOverlap(eigen_ids,
                     correlation_frame,
                     threshold,
                     go_terms_dict,
                     all_terms):
    '''
    Test significance of overlap for GO enrichment categories
    '''

    eigen_combs = itertools.combinations(eigen_ids, 2)
    fisherpy = R['fisher.test']
    q_py = R['p.adjust']
    sig_dict = {}
    for x in eigen_combs:
        if x[0] != x[1] and correlation_frame.loc[x[0]][x[1]] >= threshold:
            contingency = np.zeros((2, 2))
            g1 = set(go_terms_dict[x[0]])
            g2 = set(go_terms_dict[x[1]])

            intersect = g1.intersection(g2)
            g1_diff = g1.difference(g2)
            g2_diff = g2.difference(g1)
            all_diff = all_terms.difference(g1).difference(g2)

            contingency[0, 0] = len(intersect)
            contingency[1, 0] = len(g1_diff)
            contingency[0, 1] = len(g2_diff)
            contingency[1, 1] = len(all_diff)
            f_p = fisherpy(numpy2ri(contingency))
            f_py = [list(k) for k in np.array(f_p)]
            pvalue = f_py[0][0]
            odds = f_py[2][0]
            l_ci = f_py[1][0]
            u_ci = f_py[1][1]

            sig_dict[x] = {'OR': odds,
                           'pvalue': pvalue,
                           'lower_ci': l_ci,
                           'upper_ci': u_ci}

    sig_df = pd.DataFrame(sig_dict).T
    sig_df['qvalue'] = q_py(sig_df['pvalue'], "BH")

    # select clusters for merging with adjusted p < 0.01

    to_merge = sig_df[sig_df['qvalue'] < 0.01]
    return to_merge


def mergeClusters(eigen_file, consensus_file, go_dir, threshold):
    '''
    Merge clusters based on their functional enrichment and the
    correlation of their expression profiles.
    '''

    name = eigen_file.split("/")[1].split("-")[0]
    eigen_df = pd.read_table(eigen_file,
                             sep="\t",
                             header=0,
                             index_col=0)

    go_files = [x for x in os.listdir(go_dir) if re.search(name + "-enrich",
                                                           x)]

    # calculate all eigengene correlations
    eigen_cor = pd.DataFrame(index=eigen_df.index, columns=eigen_df.index)
    eigen_cor = eigen_cor.fillna(0.0)
    eigen_ids = itertools.combinations_with_replacement(eigen_df.index, 2)
    for each in eigen_ids:
        val1 = eigen_df.loc[each[0]].tolist()
        val2 = eigen_df.loc[each[1]].tolist()
        corr = tempCorr(val1, val2)
        eigen_cor.loc[each[0]][each[1]] = corr
        eigen_cor.loc[each[1]][each[0]] = corr

    go_terms_dict = {}
    all_terms = set()

    for fle in go_files:
        name = fle.split("_")[1].split(".")[0]
        _go_df = pd.read_table(go_dir + fle, sep="\t", header=0, index_col=1)
        _go_df.sort(columns='foldEnrichment', ascending=False)
        go_terms_dict[name] = _go_df[_go_df['padjust'] < 0.01].index.tolist()
        all_terms.update(_go_df.index.tolist())

    # test for statistically significant overlap in GO categories
    # between cluster GO enrichments with Fisher's exact test
    # only test overlap for clusters with correlation >= threshold

    sig_df = testGOCatOverlap(eigen_ids=eigen_df.index,
                              correlation_frame=eigen_cor,
                              threshold=threshold,
                              go_terms_dict=go_terms_dict,
                              all_terms=all_terms)

    return sig_df


def makeSplicedFasta(infile, outfile):
    '''
    Merge fasta sequences together into a single
    spliced transcript sequence
    '''

    fasta_dict = {}
    with IOTools.openFile(infile, "rb") as fafile:
        for line in fafile.readlines():
            if line[0] == '>':
                header = line.rstrip("\n")
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.rstrip("\n")

    with IOTools.openFile(outfile, "w") as ofile:
        for key, value in fasta_dict.items():
            ofile.write("%s\n%s\n" % (key, value))


def targetScanWrapper(miRNA_file, target_file, outfile):
    '''
    Python wrapper for MRE prediction by targetScan
    '''

    # target scan must be in the current working directoy
    assert os.path.exists("targetscan_60.pl")

    job_options = "-l mem_free=4G"

    statement = '''
    perl targetscan_60.pl %(miRNA_file)s  %(target_file)s %(outfile)s
    > %(outfile)s.log'''

    P.run()


def clusterSummary(list_of_files, outfile):
    '''
    Generate a summary table from consensus clustering
    '''

    # condition: n clusters: median objects per cluster: median length of objs
    file_dict = {}
    for fle in list_of_files:
        fname = fle.split("/")[-1]
        condition = fname.split("-")[0]
        ref = fname.split("-")[1] + "gtf.gz"
        df_ = pd.read_table(fle, sep="\t", header=0, index_col=0)
        df_.columns = ["gene_id", "cluster"]
        clust_dict = {}

        for idx in df_.index:
            cluster = df_.loc[idx]['cluster']
            try:
                clust_dict[cluster] += 1
            except KeyError:
                clust_dict[cluster] = 1

        med_size = np.median(clust_dict.value())
        file_dict[fname] = {'condition': condition,
                            'reference': ref,
                            'median_cluster_size': med_size}

    outframe = pd.DataFrame(file_dict).T
    outframe.to_csv(outfile, sep="\t", index_label="input_file")


@P.cluster_runnable
def plotClusterHeatmaps(eigengenes, expression, clusters, image_dir):
    '''
    Generate a plot of expression for each cluster with matching
    eigengene expression
    '''

    if expression.endswith("gz"):
        expr_comp = "gzip"
    else:
        expr_comp = None

    if eigengenes.endswith("gz"):
        eigen_comp = "gzip"
    else:
        eigen_comp = None

    if clusters.endswith("gz"):
        clust_comp = "gzip"
    else:
        clust_comp = None

    expr = pd.read_table(expression, sep="\t",
                         header=0, index_col=0,
                         compression=expr_comp)

    clust_df = pd.read_table(clusters, sep="\t",
                             header=None, index_col=0,
                             compression=clust_comp)
    clust_df.columns = ['gene_id', 'cluster']

    eigens = pd.read_table(eigengenes, sep="\t",
                           header=0, index_col=0,
                           compression=eigen_comp)

    mg = mygene.MyGeneInfo()
    condition = eigengenes.split("/")[-1].split("-")[0]
    reference = eigengenes.split("/")[-1].split("-")[2]

    all_clusts = set(clust_df['cluster'])
    for clust in all_clusts:
        genes = clust_df[clust_df['cluster'] == clust]['gene_id'].tolist()
        gene_expr = expr.loc[genes]
        clust_eigen = eigens.loc[clust]

        if reference == "refcoding":
            # get gene symbols - if missing replace with ensembl ID
            mg_out = mg.querymany(genes, scopes="ensemblgene", fields="symbol",
                                  species="mouse", returnall=True)['out']
            sym_df = pd.DataFrame(mg_out)
            sym_df.index = sym_df['query']
            c_df = pd.merge(left=gene_expr, right=sym_df, how='inner',
                            left_index=True, right_index=True)
            # get notfound IDs and replace with ensembl
            try:
                nf_df = c_df.loc[c_df['notfound'] == True]
                nf_df['symbol'] = nf_df['query']
                c_df.loc[nf_df.index] = nf_df
                c_df.drop(['_id', 'notfound', 'query'], inplace=True, axis=1)
            except KeyError:
                c_df.drop(['_id', 'query'], inplace=True, axis=1)

            # drop extraneous columns and remove duplicate entries based on
            # gene symbol
            c_df.index = c_df['symbol']
            c_df.drop_duplicates(subset=['symbol'],
                                 take_last=True, inplace=True)
            c_df.drop(['symbol'], inplace=True, axis=1)
            c_ids = c_df.index

        else:
            c_df = gene_expr
            c_ids = gene_expr.index

        # push objects into R and plot heatmaps
        r_ids = ro.StrVector([rs for rs in c_ids])
        r_clust = pandas2ri.py2ri_pandasdataframe(c_df)
        r_eigen = ro.FloatVector([fe for fe in clust_eigen.values])

        R.assign("gnames", r_ids)
        R.assign("gexprs", r_clust)
        R.assign("geigen", r_eigen)

        # plot heatmaps
        R('''suppressPackageStartupMessages(library(gplots))''')
        R('''suppressPackageStartupMessages(library(RColorBrewer))''')

        R('''colnames(gexprs) <- c(0,1,3,6,12,24,48,72,96,120)''')
        R('''rownames(gexprs) <- gnames''')

        # create color vector proportional to eigengene expression
        R('''eigen_func <- colorRampPalette(brewer.pal(9, "BuPu"))''')
        R('''eigen_col <- eigen_func(length(unique(geigen'''
          ''')))[as.factor(-1*geigen)]''')
        R('''hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(100)''')

        outfile = "-".join([condition, reference,
                            clust, "expression_heatmap.png"])
        outfile = image_dir + "/" + outfile
        # output to png device
        R('''png('%s', height=480, width=480)''' % outfile)
        R('''heatmap.2(as.matrix(gexprs), trace="none", col=hmcol,'''
          '''dendrogram="none",density.info="none", ColSideColors=eigen_col,'''
          '''margins=c(6,12), cexCol=2.0, labRow=F, Colv=colnames(gexprs))''')
        R('''dev.off()''')

        if reference == "refcoding":
            txt_file = outfile.rstrip("_heatmap.png")
            txt_file = txt_file + "_gene_symbols.tsv"
            # output file with gene symbol and ensembl IDs
            out_df = sym_df['symbol']
            out_df.columns = ['ensembl']
            out_df.to_csv(txt_file, sep="\t",
                          index_label="gene_symbol")
        else:
            txt_file = outfile.rstrip("_heatmap.png")
            txt_file = txt_file + "_lncRNA_ids.tsv"
            with open(txt_file, "w") as ofile:
                for lnc in c_df.index:
                    ofile.write("%s\n" % lnc)


@P.cluster_runnable
def plotEigenHeatmaps(eigengenes, image_dir):
    '''
    Plot a heatmap of eigengene correlations
    '''

    if eigengenes.endswith("gz"):
        comp = "gzip"
    else:
        comp = None
    cor_df = pd.read_table(eigengenes, sep="\t",
                           index_col=0, header=0,
                           compression=comp)
    cols = cor_df.columns
    rows = cor_df.index

    # push into R environment for plotting
    r_cols = ro.StrVector([rc for rc in cols])
    r_rows = ro.StrVector([rr for rr in rows])
    r_df = pandas2ri.py2ri_pandasdataframe(cor_df)

    R.assign("r.cols", r_cols)
    R.assign("r.rows", r_rows)
    R.assign("cor.mat", r_df)

    R('''suppressPackageStartupMessages(library(gplots))''')
    R('''suppressPackageStartupMessages(library(RColorBrewer))''')

    cond = eigengenes.split("/")[-1].split("-")[0]

    ofile = "-".join([cond, "eigengene-correlation-heatmap.png"])
    ofile = "/".join([image_dir, ofile])

    R('''rownames(cor.mat) <- r.rows''')
    R('''colnames(cor.mat) <- r.cols''')
    R('''hmcol <- colorRampPalette(brewer.pal(9, "PuOr"))(100)''')
    R('''png('%s', height=480, width=480)''' % ofile)
    R('''heatmap.2(as.matrix(cor.mat), trace="none", col=hmcol,'''
      '''dendrogram="none", density.info="none")''')
    R('''dev.off()''')


@P.cluster_runnable
def mirEnrichment(cerna_file, 
                  mre_file,
                  pairs_gtf,
                  mirna_file):
    '''
    Test for enrichment of specific miRNAs amongst ceRNAs and partner
    gene 3' UTRs

    Requirements:
        * .gtf file of MREs
        * gtf
    
    '''

    if cerna_file.endswith("gz"):
        cerna_comp = "gzip"
    else:
        cerna_comp = None

    cerna_df = pd.read_table(cerna_file, sep="\t", header=0,
                             index_col=0,
                             compression=cerna_comp)

    catalog = getMREs(mre_file, pairs_gtf)

    mirnas = set()
    with IOTools.openFile(mirna_file, "rb") as ofile:
        for line in ofile.readlines():
            mirnas.add(line.rstrip("\n"))

    fisherpy = R["fisher.test"]
    padjpy = R["p.adjust"]
   
    results = []
    pvalues = []

    lmirs = {}
    gmirs = {}
    lnc_seen = set()
    gene_seen = set()

    # create dicts of ceRNA miRs and partner miRs
    for idx in cerna_df.index:
        gene = cerna_df.loc[idx]['gene_id']
        lnc = cerna_df.loc[idx]['lncRNA_id']
        if gene not in gene_seen:
            gmirs[gene] = [g for g in catalog[gene]]
        else:
            pass

        if lnc not in lnc_seen:
            lmirs[lnc] = [l for l in catalog[lnc]]
        else:
            pass

        lnc_seen.update(lnc)
        gene_seen.update(gene)

    # generate contingency tables and test enrichment
    # of each miRNA
    for mir in mirnas:
        contingency = np.zeros((2, 2))
        for lnc in lmirs.keys():
            if mir in lmirs[lnc]:
                contingency[0, 0] += 1
            elif mir not in lmirs[lnc]:
                contingency[1, 0] += 1

        for gene in gmirs.keys():
            if mir in gmirs[gene]:
                contingency[0, 1] += 1
            elif mir not in gmirs[gene]:
                contingency[1, 1] += 1
            
        # run Fisher's exact test in R
        f = fisherpy(numpy2ri.numpy2ri(contingency), alternative="greater")

        # get overlap numbers
        ncerna = contingency[0, 0]
        npartners = contingency[0, 1]
        tcerna = contingency[0, 0] + contingency[1, 0]
        tpartners = contingency[0, 1] + contingency[1, 1]
        
        # convert fishers back to python
        fx = [list(x) for x in np.array(f)]

        # fisher.test returns pval, CIs, OR
        pvalue, ci_low, ci_hi, OR = fx[0][0], fx[1][0], fx[1][1], fx[2][0]
        pvalues.append(pvalue)

        # set default OR to 1
        if (ncerna + npartners == 0 or
            (ncerna == tcerna and npartners == tpartners)):
            OR = 1

        results.append([mir, OR, ci_low, ci_hi, pvalue, ncerna,
                        npartners, tcerna, tpartners])

    qvalues = padjpy(pvalues)
    for i in range(len(results)):
        yield results[i], [qvalues[i]]


def runSailfishIndex(fasta_file, outdir, threads,
                     kmer):
    '''
    Wrapper for sailfish index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fasta_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = '''
    sailfish index --transcripts %s --out %s --threads %i --kmerSize %i
    ''' % (fasta_file, outdir, threads, kmer)

    os.system(command)


def runSailfishQuant(fasta_index, fastq_files, output_dir,
                     paired=False, library="ISF", threads=4,
                     gene_gtf=None):
    '''
    Wrapper for sailfish quant command
    '''

    decompress = False
    if len(fastq_files) > 1:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass
    else:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " sailfish quant --index %s -l %s  -o %s " % (fasta_index,
                                                            library,
                                                            output_dir)

    states.append(command)

    if threads:
        states.append(" --threads %i " % threads)
    else:
        pass

    if gene_gtf:
        states.append(" --geneMap %s " % gene_gtf)
    else:
        pass

    # sailfish does not handle compress files natively,
    # need to decompress on the fly with advanced
    # bash syntax
    if decompress and paired:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1.gz", fq)])
        fstr_format = " ".join(["%s" for hq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

        second_mates = tuple([sq for sq in fastq_files if re.search("fastq.2.gz", sq)])
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 <(zcat %s)" % sdecomp_format

        states.append(decomp_second)

    elif decompress and not paired:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1.gz", fq)])
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -r <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

    elif paired and not decompress:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1", fq)])
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 %s " % fdecomp_format

        states.append(decomp_first)

        second_mates = tuple([sq for sq in fastq_files if re.search("fastq.2", sq)])
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 %s " % sdecomp_format

        states.append(decomp_second)

    statement = " ".join(states)

    # subprocess cannot handle process substitution
    # therefore needs to be wrapped in /bin/bash -c '...'
    # for bash to interpret the substitution correctly

    process = subprocess.Popen(statement, shell=True,
                               executable="/bin/bash")

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError(
            "-------------------------------------------\n"
            "Child was terminated by signal %i: \n"
            "The stderr was \n%s\n%s\n"
            "-------------------------------------------" %
            (-process.returncode, stderr, statement))


def runKallistoIndex(fasta_file, outfile, kmer=31):
    '''
    Wrapper for kallisto index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fast_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = "kallisto index --index=%s  %s" % (outfile,
                                                 fasta_file)

    os.system(command)


def runKallistoQuant(fasta_index, fastq_files, output_dir,
                     bias=False, bootstrap=None,
                     seed=1245, threads=None, plaintext=False):
    '''
    Wrapper for kallisto quant command
    '''

    if len(fastq_files) > 1:
        fastqs = " ".join(fastq_files)
    else:
        fastqs = fastq_files

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " kallisto quant --index=%s --output-dir=%s" % (fasta_index,
                                                              output_dir)
    states.append(command)

    if bias:
        states.append(" --use-bias ")
    else:
        pass

    if bootstrap:
        states.append(" --bootstrap=%i --seed=%i " % (bootstrap,
                                                      seed))
    else:
        pass

    if plaintext:
        states.append(" --plaintext ")
    else:
        pass

    if threads:
        states.append(" --threads=%i " % threads)
    else:
        pass

    states.append(" %s " % fastqs)

    statement = " ".join(states)

    # need to rename output files to conform to input/output
    # pattern as required.  Default name is abundance*.txt
    # when using plaintext output
    # kaliisto requires an output directory - create many small
    # directories, one for each file.
    # then extract the abundance.txt file and rename using the
    # input/output pattern

    os.system(statement)
