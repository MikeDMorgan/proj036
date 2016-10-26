'''
timeseries_analysis.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Functions and script for time-series analysis of gene expression data;
focus is on RNA-seq data.

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import glob
import gzip
import re
import math
import types
import collections
import time
import shutil
import re
import optparse
import itertools
import pandas as pd
import numpy as np
import pandas.rpy.common as com
import rpy2.rinterface as rinterface
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
import CGAT.GTF as GTF
import mygene


def DESeqNormalize(infile,
                   time_points,
                   reps,
                   conditions=None):
    '''
    variance stabilizing transformation of timeseries RNA-seq data
    '''
    
    pandas2ri.activate()
    reps = reps

    # load library
    R('''suppressMessages(library(DESeq))''')

    # generates a lists for the design data frame
    # of the proper length
    # these need to be rpy2 objects to be parsed
    # properly in the string formatting
    
    E.info("converting to pandas dataframe object")

    data_frame = pd.read_table(infile, index_col=0, header=0, sep="\t")
    rdf = com.convert_to_r_dataframe(data_frame)

    if not conditions:
        time_rep_comb = [x for x in itertools.product(time_points, reps)]
        time_cond = ro.StrVector([x[0] for x in time_rep_comb])
        rep_cond = ro.StrVector([x[1] for x in time_rep_comb])

        R.assign('countsTable', rdf)
        R('''design <- data.frame(row.names=colnames(countsTable),'''
          '''times=%s, replicates=%s)''' % (time_cond.r_repr(),
                                            rep_cond.r_repr()))
    elif conditions:
        design_dict = {}
        for x in data_frame.columns.values:
            sample_dict = {}
            sample_dict['condition'] = str(x).split(".")[0]
            sample_dict['times'] = int(str(x).split(".")[1])
            sample_dict['replicates'] = str(x).split(".")[2]
            design_dict[x] = sample_dict
        design_frame = pd.DataFrame(design_dict)
        design_frame = design_frame.T

        des_cond = design_frame['condition'].values.tolist()
        des_time = design_frame['times'].values.tolist()
        des_reps = design_frame['replicates'].values.tolist()

        cond_cond = ro.StrVector([x for x in des_cond])
        time_cond = ro.StrVector([x for x in des_time])
        rep_cond = ro.StrVector([x for x in des_reps])

        R.assign('countsTable', rdf)
        R.assign('design', design_frame)

    # create the count data set and normalize to library size
    # transform with variance stabilizing transformation
    # only select genes with an average of ten reads mapping

    E.info("calculating size factors and dispersion")

    R('''cds <- newCountDataSet(countsTable, design)''')
    R('''cds_size <- estimateSizeFactors(cds)''')
    R('''cds_disp <- estimateDispersions(cds_size, method="blind")''')
    R('''notZero <- (rowMeans(counts(cds))>10)''')

    E.info("applying variance stabilizing transformation")

    R('''vst <- varianceStabilizingTransformation(cds_disp)''')
    R('''vst <- vst[notZero, ]''')

    # format data set to long format with condition and replicate labels
    # convert to a numpy array

    R('''replicates <- c(%s)''' % rep_cond.r_repr())
    R('''times <- c(%s)''' % time_cond.r_repr())
    if conditions:
        R('''conditions <- c(%s)''' % cond_cond.r_repr())
        R('''trans_vst = data.frame(t(exprs(vst)), '''
          '''times, replicates, conditions)''')
    else:
        R('''trans_vst = data.frame(t(exprs(vst)), times, replicates)''')

    data_file = com.load_data('trans_vst')

    return data_file


def maSigPro(infile,
             order_terms=1,
             fdr=0.01,
             adjust="BH",
             stepwise="backward",
             include_p=0.01,
             rsq=0.2,
             var_group="all"):
    '''
    Generate differentially expressed genes for each experimental
    condition across a time series.  Uses the bioconductor
    package maSigPro to derive a set of genes of interest.
    '''

    ref_gtf = str(infile).split("-")[1]
    data_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    design_dict = {}

    for x in data_frame.index.values:
        sample_dict = {}
        condition = str(x).split(".")[0]
        sample_dict[condition] = 1
        sample_dict['times'] = int(str(x).split(".")[1])
        sample_dict['replicates'] = str(x).split(".")[2]
        design_dict[x] = sample_dict

    design_frame = pd.DataFrame(design_dict)
    design_frame = design_frame.T
    cols = ['times', 'replicates', condition]
    design_frame = design_frame[cols]
    design_file = "deseq.dir/%s-%s-design.tsv" % (condition, ref_gtf)
    design_frame.to_csv(design_file, sep="\t")
    data_file = "deseq.dir/%s-%s-data.tsv" % (condition, ref_gtf)
    results_file = "deseq.dir/%s-%s-maSigPro.tsv" % (condition, ref_gtf)

    # data frame columns must be in the order time-replicate-condition
    # for maSigPro
    # define the numnber of higher-order terms included in the models

    masigpro_out = "deseq.dir/maSigPro.out"
       
    R('''suppressMessages(library(maSigPro))''')
    R('''input_data <- read.table('%(infile)s', sep="\t", '''
      '''h=T, row.names=1)''' % locals())
    R('''input_data <- t(input_data[0:(length(input_data)-2)])''')

    E.info("constructing experimental design matrix")

    R('''input_design <- data.matrix(read.table('%(design_file)s', '''
      '''sep="\t", h=T, row.names=1))''' % locals())
    R('''%(condition)s_mat <- make.design.matrix(input_design, '''
      '''degree = %(order_terms)i )''' % locals())
    R('''sink(file = '%(masigpro_out)s')''' % locals())

    E.info("fitting linear model for each gene with "
           "%i polynomial terms" % order_terms)

    R('''%(condition)s_fit <- p.vector(input_data, %(condition)s_mat, '''
      '''Q = %(fdr)f, MT.adjust = '%(adjust)s')''' % locals())

    # fit a linear model to each of the genes called as
    # differentially expressed
    # report genes with model R-squared > threshold
    # maSigPro gives an un-suppressable output to stdout
    # therefore sink is used to shunt this to a temporary file 'maSigPro.out'
    
    R('''%(condition)s_step <- T.fit(%(condition)s_fit, '''
      '''step.method='%(stepwise)s', alfa=%(include_p)f)''' % locals())

    E.info("selecting significantly differentially "
           "expressed genes at FDR=%0.3f" % fdr)

    R('''sink(file=NULL)''')
    R('''%(condition)s_sigs <- get.siggenes(%(condition)s_step, '''
      '''rsq=%(rsq)f, vars='%(var_group)s')''' % locals())
    R('''write.table(%(condition)s_sigs$sig.genes$%(condition)s$group.coeffs'''
      ''',file="deseq.dir/%(condition)s-%(ref_gtf)s-coefficients.tsv", '''
      '''sep="\t")''' % locals())
    R('''write.table(%(condition)s_sigs$sig.genes$%(condition)s$sig.pvalues,'''
      '''file="deseq.dir/%(condition)s-%(ref_gtf)s-pvalues.tsv",'''
      ''' sep="\t")''' % locals())
    R('''write.table(%(condition)s_sigs$summary, '''
      '''file='deseq.dir/%(condition)s-%(ref_gtf)s-geneids.tsv', '''
      '''sep="\t")''' % locals())
    # merge the p-value and coefficient results into a single file
    p_file = "deseq.dir/%(condition)s-%(ref_gtf)s-pvalues.tsv" % locals()
    coef_file = "deseq.dir/%s-%s-coefficients.tsv" % (condition,
                                                      ref_gtf)
    p_frame = pd.read_table(p_file, sep="\t")
    coef_frame = pd.read_table(coef_file, sep="\t")
    results_frame = pd.merge(coef_frame, p_frame,
                             how='right',
                             left_index=True,
                             right_index=True)

    results_frame.to_csv(results_file, sep="\t")

    R('''diff_genes <- data.frame(%(condition)s_fit$SELEC)''' % locals())
    diff_genes = com.load_data('diff_genes')

    return diff_genes


def combined_maSigPro(infile,
                      order_terms=1,
                      fdr=0.05,
                      padjust="BH",
                      stepwise="backward",
                      pinclude=0.01,
                      rsquared=0.2,
                      var_group="groups"):

    ref_gtf = str(infile).split("-")[1]
    data_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    design_dict = {}

    # generate the experimental design data frame
    for x in data_frame.index.values:
        sample_dict = {}
        for cond in set(data_frame['conditions']):
            if cond == data_frame.loc[x][-1]:
                sample_dict[cond] = 1
            else:
                sample_dict[cond] = 0

        sample_dict['times'] = int(str(x).split(".")[1])
        sample_dict['replicates'] = str(x).split(".")[2]
        design_dict[x] = sample_dict

    # set the columns in order of time-replicate-conditions

    design_frame = pd.DataFrame(design_dict)
    design_frame = design_frame.T
    cols = list(set(data_frame['conditions']))
    cols.sort()
    cols.insert(0, 'times')
    cols.insert(1, 'replicates')

    design_frame = design_frame.T
    design_frame = design_frame[cols]

    # sort so that LPS is the control group

    design_frame = design_frame.sort(columns=['LPS', 'times'],
                                     ascending=[0, 1])

    results_file = "combined_analysis.dir/combined-%s-maSigPro.tsv" % ref_gtf

    # define the numnber of higher-order terms included in the models

    masigpro_out = "combined_analysisdeseq.dir/%s-maSigPro.out" % ref_gtf

    R.assign('input_data', data_frame)
    R.assign('input_design', design_frame)
    
    R('''suppressMessages(library(maSigPro))''')
    R('''input_data <- t(input_data[0:(length(input_data)-3)])''')

    E.info("constructing experimental design matrix")

    R('''combined_mat <- make.design.matrix(input_design,'''
      ''' degree = %(order_terms)i )''')
    R('''sink(file = '%(masigpro_out)s')''' % locals())

    E.info("fitting linear model for each gene with %i"
           " polynomial terms" % order_terms)

    R('''combined_fit <- p.vector(input_data, combined_mat,'''
      ''' Q = %(fdr)f, MT.adjust = '%(adjust)s')''' % locals())

    # fit a linear model to each of the genes called as
    # differentially expressed
    # report genes with model R-squared > threshold
    # maSigPro gives an un-suppressable output to stdout
    # therefore sink is used to shunt this to a temporary file 'maSigPro.out'
    
    R('''combined_step <- T.fit(combined_fit, step.method='%(stepwise)s', '''
      '''alfa=%(include_p)f)''' % locals())

    E.info("selecting significantly differentially "
           "expressed genes at FDR=%0.3f" % fdr)

    R('''sink(file=NULL)''')
    R('''combined_sigs <- get.siggenes(combined_step, '''
      '''rsq=%(rsq)f, vars='%(var_group)s')''' % locals())

    R('''write.table(combined_sigs$sig.genes$%(condition)s$group.coeffs,'''
      '''file="deseq.dir/%(condition)s-%(ref_gtf)s-coefficients.tsv",'''
      '''sep="\t")''' % locals())
    R('''write.table(%(condition)s_sigs$sig.genes$%(condition)s$sig.pvalues,'''
      '''file="deseq.dir/%(condition)s-%(ref_gtf)s-pvalues.tsv",'''
      ''' sep="\t")''' % locals())
    R('''write.table(%(condition)s_sigs$summary,'''
      '''file='deseq.dir/%(condition)s-%(ref_gtf)s-geneids.tsv','''
      ''' sep="\t")''' % locals())
    # merge the p-value and coefficient results into a single file
    p_file = "deseq.dir/%(condition)s-%(ref_gtf)s-pvalues.tsv" % locals()
    coef_file = "deseq.dir/%s-%s-coefficients.tsv" % (condition,
                                                      ref_gtf)

    p_frame = pd.read_table(p_file, sep="\t")
    coef_frame = pd.read_table(coef_file, sep="\t")
    results_frame = pd.merge(coef_frame, p_frame,
                             how='right',
                             left_index=True,
                             right_index=True)
    results_frame.to_csv(results_file, sep="\t")

    R('''diff_genes <- data.frame(%(condition)s_fit$SELEC)''' % locals())
    diff_genes = com.load_data('diff_genes')

    return diff_genes


def CovarFilter(infile,
                time_points,
                replicates,
                quantile):
    '''
    Filter gene list based on the distribution of the
    sums of the covariance of each gene.  This is highly
    recommended to reduce the total number of genes used
    in the dynamic time warping clustering to reduce the
    computational time.  The threshold is placed at the
    intersection of the expected and observed value
    for the given quantile.
    '''

    time_points.sort()
    time_rep_comb = [x for x in itertools.product(time_points, replicates)]
    time_cond = ro.StrVector([x[0] for x in time_rep_comb])
    rep_cond = ro.StrVector([x[1] for x in time_rep_comb])

    E.info("loading data frame")

    R('''diff_data <- read.table('%(infile)s', h=T, row.names=1)''' % locals())
    
    # need to be careful about column headers and transposing data frames

    R('''trans_data <- data.frame(t(diff_data))''')
    R('''times <- c(%s)''' % time_cond.r_repr())
    R('''replicates <- c(%s)''' % rep_cond.r_repr())
    
    # calculate the covariance matrix for all genes
    # sum each gene's covariance vector

    E.info("calculating sum of covariance of expression")

    R('''covar.mat <- abs(cov(trans_data))''')
    R('''sum.covar <- rowSums(covar.mat)''')
    R('''exp.covar <- abs(qnorm(ppoints(sum.covar),'''
      '''mean=mean(sum.covar), sd=sd(sum.covar)))''')
    R('''sum.covar.quant <- quantile(sum.covar)''')
    R('''exp.covar.quant <- quantile(exp.covar)''')

    E.info("filter on quantile")

    R('''filtered_genes <- names(sum.covar[sum.covar > '''
      '''sum.covar.quant[%(quantile)i]'''
      ''' & sum.covar > exp.covar.quant[%(quantile)i]])''' % locals())
    R('''filtered_data <- data.frame(diff_data[filtered_genes, ])''')
    R('''filtered_frame <- data.frame(t(diff_data[filtered_genes, ]),'''
      '''times, replicates)''')

    filtered_frame = com.load_data('filtered_frame').T

    return filtered_frame


def treeCutting(infile,
                expression_file,
                cluster_file,
                cluster_algorithm):
    '''
    Use dynamic tree cutting to derive clusters for each
    resampled distance matrix
    '''
    wgcna_out = "tmp.dir/WGCNA.out"

    R('''sink(file='%(wgcna_out)s')''' % locals())
    R('''suppressMessages(library(WGCNA, flashClust))''')

    E.info("loading distance matrix")

    R('''distance_data <- data.matrix(read.table('%(infile)s', '''
      '''h=T, row.names=1))''' % locals())

    E.info("clustering data by %s linkage" % cluster_algorithm)

    R('''clustering <- flashClust(as.dist(distance_data),'''
      ''' method='%(cluster_algorithm)s')''' % locals())
    R('''cluster_cut <- cutreeDynamic(dendro=clustering, '''
      '''minClusterSize=50, deepSplit=2)''')
    R('''color_cut <- labels2colors(cluster_cut)''')
    R('''write.table(color_cut, file = '%(cluster_file)s','''
      '''sep="\t")''' % locals())
    R('''cluster_matched <- data.frame(cbind(rownames(distance_data),'''
      '''color_cut))''')
    R('''colnames(cluster_matched) = c("gene_id", "cluster")''')
    R('''cluster_matched <- data.frame(cluster_matched$gene_id,'''
      '''cluster_matched$cluster)''')
    R('''sink(file=NULL)''')

    cluster_frame = com.load_data('cluster_matched')
    
    return cluster_frame


def consensusClustering(infile,
                        cutHeight,
                        cluster_algorithm):
    '''
    hierachichal clustering based on gene-cluster correlation across
    resampled datasets.  cut tree based on arbitrary threshold, i.e.
    the proportion of resampled datasets where two genes fall into
    the same cluster
    '''
    condition = infile.split("/")[1].split("-")[0]
    wgcna_out = "tmp.dir/consensus-WGCNA.out"

    R('''sink(file='%(wgcna_out)s')''' % locals())
    R('''suppressMessages(library(WGCNA, flashClust))''')

    E.info("loading distance matrix")

    R('''distance_data <- data.matrix(read.table('%(infile)s','''
      '''h=T, row.names=1))''' % locals())

    E.info("clustering data by %s linkage" % cluster_algorithm)

    R('''clustering <- flashClust(as.dist(1-distance_data),'''
      '''method='%(cluster_algorithm)s')''' % locals())
    if cutHeight > float(0.01):
        R('''cluster_cut <- cutreeDynamic(dendro=clustering, '''
          '''distM=as.matrix(1-distance_data), minClusterSize=50,'''
          '''cutHeight=%(cutHeight)s)''' % locals())
    else:
        R('''cluster_cut <- cutreeDynamic(dendro=clustering, '''
          '''distM=as.matrix(1-distance_data), '''
          '''minClusterSize=50)''' % locals())
    R('''color_cut <- labels2colors(cluster_cut)''')
    R('''cluster_matched <- data.frame(cbind(rownames(distance_data),'''
      '''color_cut))''')
    R('''colnames(cluster_matched) = c("gene_id", "cluster")''')
    R('''cluster_matched <- data.frame(cluster_matched$gene_id,'''
      '''cluster_matched$cluster)''')
    
    # plot and save dendrogram of clustering

    R('''png("plots.dir/%(condition)s-dendrogram-consensus_clustering.png")'''
      % locals())
    R('''plotDendroAndColors(dendro=clustering, colors=color_cut,'''
      '''groupLabels="Dynamic tree cut",'''
      '''dendroLabels=F, addGuide=T, guideHang=0.05, '''
      '''hang=0.03, main="%(condition)s")''' % locals())
    R('''dev.off()''')
    R('''sink(file=NULL)''')

    cluster_frame = com.load_data('cluster_matched')
    
    return cluster_frame


def eigenGenes(infile,
               cluster_file):
    '''
    Use WGCNA to derive eigengenes expression profiles
    from clustering and tree cutting
    '''
    
    header = cluster_file.split("/")[1].split("-")[0]

    # reshape data and derive module eigengenes average expression
    R('''sink(file='sink_file.txt')''')
    R('''suppressMessages(library(reshape2))''')
    R('''suppressMessages(library(WGCNA))''')
    R('''source("/ifs/devel/michaelm/Time-series/summarySE.R")''')
    R('''cluster_match <- read.table('%(cluster_file)s', h=T)''' % locals())
    R('''express_data <- read.table('%(infile)s', '''
      '''h=T, row.names=1)''' % locals())
    R('''data_colour = cluster_match[,2]''')
    R('''data_melt <- melt(express_data, id.vars=c("times", "replicates"))''')
    R('''data_sum <- summarySE(data_melt, measurevar="value",'''
      '''groupvars=c("times", "variable"))''')
    R('''data_mod <- data.frame(data_sum$times,'''
      ''' data_sum$variable, data_sum$value)''')
    R('''colnames(data_mod) <- c("times", "gene", "value")''')
    R('''data_wide <- dcast(data_mod, gene ~ times, value.var="value")''')
    R('''rownames(data_wide) <- data_wide$gene''')
    R('''times <- unique(express_data$times)''')
    R('''data_wide <- data.frame(data_wide[,-1])''')
    R('''sink(file=NULL)''')
    
    # derive module eigengenes based on consensus clustering
    R('''MElist <- moduleEigengenes(t(data_wide), colors=data_colour)''')
    R('''eigen_express <- data.frame(MElist$averageExpr, times)''')
    R('''varExp <- data.frame(MElist$varExplained)''')
        
    var_explained = com.load_data('varExp')
    var_file = "eigengenes.dir/%s-eigengene_varExplained.tsv" % header
    var_explained.to_csv(var_file,
                         sep="\t")
    eigen_expr = com.load_data('eigen_express')

    return eigen_expr


def clusterPCA(infile,
               cluster_file,
               image_dir):
    '''
    PCA for each module within an experimental condition across
    the time series.
    Take PC1 as the module eigengene and return the loadings and proportion
    of variance explained for the eigengene.
    The eigengene expression across the time series is taken to be the
    value for PC1 at each timepoint as a vector.
    This is basically what WGCNA moduleEigengenes does but it
    does not recover the PC loadings.
    '''
    
    header = cluster_file.split("-")[0]

    # reshape data
    R('''sink(file='sink_file.txt')''')
    R('''suppressMessages(library(reshape2))''')
    R('''suppressMessages(library(WGCNA))''')
    R('''source("/ifs/devel/michaelm/Time-series/summarySE.R")''')
    R('''source("Rscripts/clusterEigengenes.R)''')
    R('''cluster_match <- read.table('%(cluster_file)s', h=T)''' % locals())
    R('''express_data <- read.table('%(infile)s', '''
      '''h=T, row.names=1)''' % locals())
    R('''sink(file=NULL)''')
    R('''print(head(express_data))''')
    R('''colnames(cluster_match) <- c("genes", "cluster")''')
    R('''data_melt <- melt(express_data, id.vars=c("times", "replicates"))''')
    R('''data_sum <- summarySE(data_melt, measurevar="value", '''
      '''groupvars=c("times", "variable"))''')
    R('''data_mod <- data.frame(data_sum$times,'''
      ''' data_sum$variable, data_sum$value)''')
    R('''colnames(data_mod) <- c("times", "gene", "value")''')
    R('''data_wide <- dcast(data_mod, gene ~ times, value.var="value")''')
    R('''rownames(data_wide) <- data_wide$gene''')
    R('''times <- as.numeric(unique(express_data$times))''')
    R('''data_wide <- data.frame(data_wide[,-1])''')
    
    # derive module eigengenes - return a dataframe of eigengene expression

    R('''eigen_frame <- eigenExpress(clusterPCA(cluster_frame=cluster_match,'''
      '''expression_frame=data_wide, n=times), n=times)''')
    
    # generate loadings plot for each eigengene

    R('''eigenLoad(clusterPCA(cluster_frame=cluster_match,'''
      '''expression_frame=data_wide, n=times), image.dir=%(image_dir)s,'''
      '''condition=%(header)s)''' % locals())

    # generate expression profile plots for all eigengenes

    R('''eigenPlot(eigen_frame, image.dir=%(image_dir)s,'''
      '''condition=%(header)s)''' % locals())
    
    eigen_frame = com.load_data("eigen_frame")

    return eigen_frame


def pairwise_timecourse(infile,
                        fdr="BH"):
    '''
    Use the empirical Bayes framework in the R timecourse package
    to calculate Hotellings T2 for differentially expressed genes
    between two conditions.
    Calculate p-values based on the equivalent F-distribution
    with appropriate FDR adjustment.
    '''

    def F_stats(data_frame,
                hotellings_frame,
                header):

        '''
        Returns the equivalent F-value from Hotellings T2 and
        the exact p-value given F(d1, d2).
        FDR adjustment according to Benjamani and Hochberg
        '''

        stat_dict = {}
        nx = len(data_frame['conditions'])/2.0
        ny = len(data_frame['conditions'])/2.0
        
        p = len(set(data_frame['times']))
        f_dict = {}
        p_dict = {}

        pf = ro.r.pf
        for gene in hotellings_frame.index.values:
            T2 = float(hotellings_frame.loc[gene]['HotellingsT2_%s' % header])
            fstat = ((nx+ny-p-1.0)/float((nx+ny-2.0) * p)) * float(T2)
            df2 = p
            df1 = nx + ny - 1.0 - p
            p_value = 1.0 - float(pf(float(fstat), df1, df2)[0])
            f_dict[gene] = fstat
            p_dict[gene] = p_value

        adjust = ro.r('p.adjust')
        p_series = pd.Series(p_dict)
        p_copy = p_series.copy()
        p_vector = ro.FloatVector(p_copy.values.tolist())
        p_adjust = adjust(p_vector,
                          method="bonferroni")
        p_adjust = pd.Series(p_adjust,
                             index=p_series.index.values)
        f_series = pd.Series(f_dict)

        stat_dict['fstat_%s' % header] = f_series
        stat_dict['pvalue_%s' % header] = p_series
        stat_dict['padjust_%s' % header] = p_adjust
        
        stat_frame = pd.DataFrame(stat_dict)
        out_frame = pd.merge(hotellings_frame,
                             stat_frame,
                             left_index=True,
                             right_index=True)

        return out_frame

    combined_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    analysis_group = infile.split("-")[0].split("/")[1]

    # generate all pairwise combinations of experimental conditions
    comb_cond = set(combined_frame['conditions'])
    pairwise = [x for x in itertools.combinations(comb_cond, 2)]
    condition_set = set(combined_frame['conditions'])

    # create one combined design frame to sample from

    design_dict = {}
    for x in combined_frame.index.values:
        sample_dict = {}
        for cond in condition_set:
            if cond == combined_frame.loc[x][-1]:
                sample_dict[cond] = 1
            else:
                sample_dict[cond] = 0
        sample_dict['times'] = int(str(x).split(".")[1])
        sample_dict['replicates'] = str(x).split(".")[2]
        design_dict[x] = sample_dict

        design_frame = pd.DataFrame(design_dict).T

    # create a dictionary of design frames for each pairwise combination

    sub_design_frame = {}

    for pair in pairwise:
        design1 = design_frame[design_frame[pair[0]] == 1]
        design2 = design_frame[design_frame[pair[1]] == 1]
        designs = [design1, design2]
        
        concated = pd.concat(designs)
        for cond in condition_set:
            if cond not in pair:
                del concated[cond]
        sub_design_frame[pair] = concated

    # create a dictionary of expression data frames for
    # each pairwise combination

    groups = combined_frame.groupby('conditions')

    sub_data_frame = {}
    for pair in pairwise:
        frame1 = combined_frame.loc[groups.groups[pair[0]]]
        frame2 = combined_frame.loc[groups.groups[pair[1]]]

        frames = [frame1, frame2]
        concated = pd.concat(frames)
        sub_data_frame[pair] = concated

    # iterate over all pair-wise combinations of experimental conditions
    # and calculate Hotellings T2, F-statistics and p-values for
    # differentially expressed genes
    # report all results

    pair_dict = {}
    R('''suppressPackageStartupMessages(library(timecourse))''')
    for pair in pairwise:

        input_data = pd.DataFrame(sub_data_frame[pair])
        timecourse_design  = sub_design_frame[pair]
        timecourse_design = timecourse_design.sort([pair[0],
                                                    'replicates'],
                                                   ascending=[0,1])
        
        timecourse_data = input_data.T[timecourse_design.index.values.tolist()]

        reps = ro.StrVector([str(x) for x in timecourse_data.loc['replicates'].values])
        times = ro.IntVector([int(x) for x in timecourse_data.loc['times'].values])

        timecourse_data = timecourse_data.drop(['times'])
        timecourse_data = timecourse_data.drop(['replicates'])
        timecourse_data = timecourse_data.drop(['conditions'])
        
        input_file = "combined_analysis.dir/%s_%svs%s-timecourse-data.tsv" % (analysis_group,
                                                                              pair[1],
                                                                              pair[0])
        con_head = "%svs%s" % (pair[1],pair[0])

        timecourse_data.to_csv(input_file, sep="\t")

        # calculate Hotellings T2 for differentially expressed genes between conditions

        R('''sink(file='sink_file.txt')''')
        R('''input_data <- read.table('%(input_file)s', h=T, row.names=1)''' % locals())
        R('''trt <- rep(c('%s', '%s'), each = 30)''' % (pair[1], pair[0]))
        R('''reps <- rep(rep(c('%(reps)s'), each=10), each=2)''' % locals())
        R('''times <- rep(rep(as.numeric('%(times)s'),3),2)''' % locals())
        R('''sizes <- matrix(3, nrow=dim(input_data)[1], ncol=2)''')
        R('''MB.paired <- mb.long(input_data, method="2D", times = 10, reps=sizes,'''
          '''condition.grp=trt)''')
        R('''hotellings <- data.frame(cbind(rownames(MB.paired$M),'''
          '''MB.paired$HotellingT2, MB.paired$pos.HotellingT2))''')
        R('''rownames(hotellings) <- hotellings$X1''')
        R('''hotellings <- hotellings[,-1]''')
        R('''colnames(hotellings) <- c("HotellingsT2_%(con_head)s", "rank_%(con_head)s")''' % locals())
        R('''sink(file=NULL)''')
        hotellings_frame = com.load_data('hotellings')

        f_frame = F_stats(data_frame=input_data,
                          hotellings_frame=hotellings_frame,
                          header=con_head)

        f_frame.to_csv("combined_analysis.dir/%s_%svs%s-hotellings.tsv" % (analysis_group,
                                                                           pair[1],
                                                                           pair[0]), sep="\t")
        pair_dict[pair] = f_frame

    return pair_dict


def clusterAgreement(infile):
    '''
    calculate co-occurence of genes within resampled clusters
    '''
    
    # read in aggregated cluster assignment file, genes as rows, iterations as columns

    df = pd.read_table(infile, sep="\t", header=None, index_col=0)
    genes = df.index.values

    # instantiate an empy matrix to count the number of times each gene appears with others
    # at each iteration

    dmat = pd.DataFrame(index=genes,
                        columns=genes)
    dmat = dmat.fillna(0)
    
    # generate all pair-wise combinations of genes, index into dmat using these
    reps = df.columns.values

    #combinations = itertools.combinations_with_replacement(genes, 2)

    '''
    for x in combinations:
        for k in reps:
            if df.loc[x[0]][k] == df.loc[x[1]][k]:
                if x[0] == x[1]:
                    dmat[x[0]][x[1]] += 1
                else:
                    dmat[x[0]][x[1]] += 1
                    dmat[x[1]][x[0]] += 1
    
    '''

    # alternative, faster code for consensus clustering - might be able to improve
    # it by using less nested for loops
    
    # generate a set for each cluster that contains the genes belonging to each cluster
    # count the number of times two genes occur in a cluster and add to the dataframe
    # repeat for each resampling iteration
    # cluster sets are generated at every iteration - time improvement over directly
    # accessing and storing as two dataframes is a factor of 7-8

    for i in reps:
        # from the input dataframe generate a list of sets, one set for each cluster
        # in that resampling iteration.
        # repeat for every resampling iteration.
 
        clusters = set(df[i].values.tolist())
        combinations = itertools.combinations_with_replacement(genes, 2)
        cluster_dict = {}
        for col in clusters:
            cluster_dict[col] = []
        for gene in genes:
            cluster_dict[df[i][gene]].append(gene)
        rep_list = []
        
        # for each cluster add all the genes with that cluster ID to the set
        # add all of the cluster sets to a list container

        for col in clusters:
            col_set  = set()
            gene_members = itertools.combinations_with_replacement(cluster_dict[col], 2)
            col_set.add(gene_members)
            rep_list.append(col_set)
    
        # count if two genes occur in the same cluster, adding to both sides
        # of the symmetrical dataframe

        for cluster_set in rep_list:
            for combs in cluster_set:
                for x in combs:
                    if x[0] == x[1]:
                        dmat[x[0]][x[1]] += 1
                    else:
                        dmat[x[0]][x[1]] += 1
                        dmat[x[1]][x[0]] += 1

    # calculate the proportion of co-occurences

    prob = lambda x: x/float(len(reps))

    probs_df = dmat.applymap(prob)

    return probs_df


def correlate_lncsEigen(lncFile, eigenFile):
    '''
    correlate lncRNA expression against module eigengene expression
    '''

    # average lncRNA expression across replicates and correlate with eigengene
    # expression
    
    sinkfile = P.getTempFilename(".")
    
    R('''sink(file='%(sinkfile)s')''' % locals())
    R('''suppressMessages(library(reshape2))''')
    R('''source("/ifs/devel/michaelm/Time-series/summarySE.R")''')
    R('''lncs_data <- read.table('%(lncFile)s', sep="\t", h=T, row.names=1)''' % locals())
    R('''sink(file=NULL)''')
    R('''eigen_data <- read.table('%(eigenFile)s', sep="\t", h=T, row.names=1)''' % locals())
    R('''lncs_melt <- melt(lncs_data, id.vars=c("times", "replicates"))''')
    R('''lncs_sum <- summarySE(lncs_melt, measurevar="value", groupvars=c("times", "variable"))''')
    R('''lncs_mod <- data.frame(lncs_sum$times, lncs_sum$variable, lncs_sum$value)''')
    R('''colnames(lncs_mod) <- c("times", "variable", "value")''')
    R('''lncs_wide <- dcast(lncs_mod, variable ~ times, value.var="value")''')
    R('''rownames(lncs_wide) <- lncs_wide$variable''')
    R('''lncs_wide <- lncs_wide[,-1]''')
    
    R('''lncs_eigen_cor <- cor(eigen_data[, -(length(eigen_data))], t(lncs_wide))''')
    
    cor_frame = com.load_data('lncs_eigen_cor')
    
    return cor_frame


def filter_correlation(infile,
                       threshold):
    '''
    filter lncRNAs with correlation below threshold
    '''
    
    cor_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    
    lncs_dict = {}

    modules = cor_frame.index.values.tolist()

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
    output_frame = output_frame.replace("NaN", 0.0)

    return output_frame


def classify_lncrna(infile,
                    lnc_gtf,
                    summary_file,
                    out_gtf):
    '''
    classify lncRNAs based on the direction of their correlation with
    module eigengenes
    '''

    lnc_cor_dict = {}
    mod_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    modules = mod_frame.index.values.tolist()
    
    lnc_file = IOTools.openFile(lnc_gtf, "r")

    with open(out_gtf, "w") as gtf_file:
        for line in GTF.readFromFile(lnc_file):
            entry = GTF.Entry()
            if line.gene_id in mod_frame:
                entry.copy(line)
                for mod in modules:
                    entry.addAttribute(key=mod,
                                       value=str(mod_frame[entry.gene_id][mod]))
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
                        neg_lncs[entry.gene_id] = {'source' : entry.source,
                                                   mod : 'negative'}

                    elif float(entry.asDict()[mod]) > 0:
                        pos_lncs[entry.gene_id] = {'source' : entry.source,
                                                   mod : 'positive'}
                except(KeyError):
                    pass

    pos_series = pd.Series(pos_lncs)
    neg_series = pd.Series(neg_lncs)
    
    all_series = pos_series.append(neg_series)

    return all_series


def avTimeExpression(infile):
    '''
    Calculate average expression over replicates at each time point
    '''

    # import into R and reformat with reshape, average using
    # custom R function

    #R('''sink(file='sink_file.txt')''')
    R('''source("/ifs/devel/michaelm/Time-series/summarySE.R")''')
    R('''library(reshape2)''')
    R('''indata <- read.table('%(infile)s', sep="\t", '''
      '''row.names=1, h=T)''' % locals())

    # make sure times is numeric rather than factor
    R('''times <- sapply(colnames(indata), '''
      '''function(x) unlist(strsplit(x, "[.]"))[2])''')
    R('''replicates <- sapply(colnames(indata), '''
      '''function(x) unlist(strsplit(x, "[.]"))[3])''')
    R('''tdata = data.frame(t(indata))''')
    R('''tdata$times <- as.numeric(as.character(times))''')
    R('''tdata$replicates <- replicates''')

    # melt data frame and calculate average over replicates

    R('''data_melted <- melt(tdata, id.vars=c("times", "replicates"))''')
    R('''data_sum <- summarySE(data_melted, measurevar="value", groupvars=c("times", "variable"))''')
    R('''data_sum <- data.frame(data_sum$times, data_sum$variable, data_sum$value)''')
    R('''colnames(data_sum) <- c("times", "genes", "value")''')
    R('''data_av <- dcast(data_sum, genes ~ times, value.var="value")''')
    R('''rownames(data_av) <- data_av$genes''')
    R('''data_av <- data_av[,-1]''')
    R('''sink(file=NULL)''')
    
    data_frame = com.load_data('data_av')
    return data_frame

    
def classCorrLncRNA(cluster_file,
                    gene_express,
                    lncRNA_express,
                    lnc_gtf,
                    threshold,
                    lncrna_class,
                    corr_direction):
    '''
    Classify :lncrna_class: lncRNAs based on :corr_direction: correlation with protein
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
    lnc_class_frame = pd.merge(class_frame,
                               lnc_express,
                               how="right",
                               left_index=True,
                               right_index=True)

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
                if any([True for x in cor_frame.loc[lncrna] if x < -threshold]):
                    correlated_lncs.append(lncrna)
                else:
                    pass

        lncs_cor_frame = cor_frame.loc[correlated_lncs]
        correlation_dict[col] = correlated_lncs

        # write out each correlation matrix to a separate file with cluster ID
        # write out list of correlated lncRNA IDs to file

        correlation_out = "lncRNA_classification.dir/%s-%s-%s-correlations.tsv" % (file_prefix,
                                                                                   col,
                                                                                   corr_direction)
        lncs_cor_frame.to_csv(correlation_out,
                              sep="\t")
        correlated_lncs_out = "lncRNA_classification.dir/%s-%s-%s-lncRNAs.tsv" % (file_prefix,
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


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--time", dest="timepoints", type="string",
                      help="a comma-separated list of time points measured")

    parser.add_option("--replicates", dest="reps", type="string",
                      help="a comma-separated list of replicate IDs")

    parser.add_option("--conditions", dest="conditions", type="string",
                      help="a comma-separated list of experimental conditions")
    
    parser.add_option("--orders", dest="orders", type="int",
                      help="order of polynomial terms to include in"
                      "maSigPro linear model")

    parser.add_option("--fdr", dest="fdr", type="string",
                      help="FDR for calling DEGs")

    parser.add_option("--padjust", dest="padjust", type="string",
                      help="multiple testing correction to apply to"
                      "control FDR")

    parser.add_option("--stepwise", dest="stepwise", type="string",
                      help="stepwise regression to use")

    parser.add_option("--pinclude", dest="pinclude", type="string",
                      help="p-value for inclusion in stepwise regression")

    parser.add_option("--rsquared", dest="rsquared", type="string",
                      help="rsquared cut-off for DEG reporting")

    parser.add_option("--var-group", dest="vargroup", type="string",
                      help="variable group reporting. each, all or"
                      "group")

    parser.add_option("--task", dest="task", type="string",
                      help="analysis task to be executed")

    parser.add_option("--infile", dest="infile", type="string",
                      help="input file path")
    
    parser.add_option("--quantile", dest="quantile", type="int",
                      help="see pipeline.ini for explanation")

    parser.add_option("--cluster-algorithm", dest="cluster", type="string",
                      help="hierarchical clustering algorithm")

    parser.add_option("--expression-file", dest="express", type="string",
                      help="matching expression data from input distance matrix")

    parser.add_option("--cluster-file", dest="clustfile", type="string",
                      help="file to output cluster labels to")

    parser.add_option("--output-file", dest="outfile", type="string",
                      help="output file to write to")

    parser.add_option("--cut-height", dest="cutHeight", type="string",
                      help="threshold at which to define consensus clusters"
                      "as valid")
    
    parser.add_option("--cor-threshold", dest="cor_threshold",type="string",
                      help="threshold at which to filter lncRNA expression"
                      "correlation with eigengene expression")

    parser.add_option("--summary-file", dest="summary_file", type="string",
                      help="summary file")

    parser.add_option("--output-gtf", dest="output_gtf", type="string",
                      help="output gtf file to write attributes to")

    parser.add_option("--image-dir", dest="images_dir", type="string",
                      help="directory to write plots/figures to")

    parser.add_option("--cor-direction", dest="corr_direction", type="choice",
                      choices=("positive", "negative"), help="direction of "
                      "correlation to classify lncRNAs on")

    parser.add_option("--lncRNA-class", dest="lncrna_class", type="choice",
                      choices=("intergenic",
                               "antisense",
                               "sense_downstream",
                               "sense_upstream",
                               "antisense_upstream",
                               "antisense_downstream"),
                      help="classification of lncRNA to use for correlation analysis")
                               

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    parser.set_defaults(cutHeight=0,
                        conditions=None)
        
    if options.task == "deseq":
        timepoints = [int(x) for x in options.timepoints.split(",")]
        timepoints.sort()
        reps = [x for x in options.reps.split(",")]
        if options.conditions == None:
            conditions = None
        else:
            conditions = [x for x in options.conditions.split(",")]
        infile = options.stdin
        data_frame = DESeqNormalize(infile=infile,
                                    time_points=timepoints,
                                    reps=reps,
                                    conditions=conditions)

    elif options.task == "masigpro":
        data_frame = maSigPro(infile=infile,
                              order_terms=int(options.orders),
                              fdr=float(options.fdr),
                              adjust=options.padjust,
                              stepwise=options.stepwise,
                              include_p=float(options.pinclude),
                              rsq=float(options.rsquared),
                              var_group=options.vargroup)

    elif options.task == "sumcovar":
        timepoints = [int(x) for x in options.timepoints.split(",")]
        reps = [x for x in options.reps.split(",")]
        data_frame = CovarFilter(infile=infile,
                                 time_points=timepoints,
                                 replicates=reps,
                                 quantile=int(options.quantile))

    elif options.task == "cluster":
        data_frame = treeCutting(infile=infile,
                                 expression_file=options.express,
                                 cluster_file=options.clustfile,
                                 cluster_algorithm=options.cluster)

    elif options.task == "clustagree":
        data_frame = clusterAgreement(infile)

    elif options.task == "consensus-cluster":
        data_frame = consensusClustering(infile=infile,
                                         cutHeight=float(options.cutHeight),
                                         cluster_algorithm=options.cluster)
    elif options.task == "eigengene":
        files = infile.split(",")
        infile = files[1]
        cluster_file = files[0]
        data_frame = eigenGenes(infile=infile,
                                cluster_file=cluster_file)

    elif options.task == "timecourse":
        data_dict = pairwise_timecourse(infile,
                                        fdr=float(options.fdr))
        frame_list = []
        for key in data_dict.keys():
            frame_list.append(pd.DataFrame(data_dict[key]))
        frame1 = pd.merge(frame_list[0],
                          frame_list[1],
                          how="right",
                          left_index=True,
                          right_index=True)
        frame2 = pd.merge(frame1,
                          frame_list[2],
                          how="right",
                          left_index=True,
                          right_index=True)
        data_frame = frame2

    elif options.task == "correlate-lncs":
        files = infile.split(",")
        eigenFile = files[0]
        lncFile = files[1]

        data_frame = correlate_lncsEigen(lncFile=lncFile,
                                         eigenFile=eigenFile)
    elif options.task == "filter-correlation":
        data_frame = filter_correlation(infile=infile,
                                        threshold=float(options.cor_threshold))

    elif options.task == "classify-lncrna":
        files = infile.split(",")
        infile = files[0]
        lnc_gtf = files[1]
        summary_file = options.summary_file
        out_gtf = options.output_gtf

        data_frame = classify_lncrna(infile=infile,
                                      lnc_gtf=lnc_gtf,
                                      summary_file=summary_file,
                                      out_gtf=out_gtf)
    elif options.task == "pca":
        files = infile.split(",")
        infile = files[1]
        cluster_file = files[0]
        data_frame = clusterPCA(infile=infile,
                                cluster_file=cluster_file,
                                image_dir=options.images_dir)

    elif options.task == "average_expression":
        data_frame = avTimeExpression(infile)

    elif options.task == "correlate_direction":

        files = infile.split(",")
        cluster_file = files[0]
        gene_express = files[1]
        lncRNA_express = files[2]
        lnc_gtf = files[3]

        threshold = float(options.cor_threshold)
        lncrna_class = options.lncrna_class
        corr_direction = options.corr_direction

        data_frame = classCorrLncRNA(cluster_file,
                                     gene_express,
                                     lncRNA_express,
                                     lnc_gtf,
                                     threshold,
                                     lncrna_class,
                                     corr_direction)
    else:
        pass

    data_frame.to_csv(options.stdout, sep="\t",
                      header=True, index_label="gene_id")

    # Write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
