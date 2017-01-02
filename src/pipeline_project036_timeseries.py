################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
Pipeline project036_timeseries
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline template.

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_project036_timeseries.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_project036_timeseries.tgz
   tar -xvzf pipeline_project036_timeseries.tgz
   cd pipeline_project036_timeseries
   python <srcdir>/pipeline_project036_timeseries.py make full

.. note::
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *
from ruffus.combinatorics import *
import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3
import random, itertools
import tempfile
import numpy as np
import pandas as pd
from pandas.io import sql
import rpy2.rinterface as rinterface
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.GTF as GTF
import PipelineProject036

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# sample = PipelineTracks.AutoSample

# assign tracks
GENESETS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("*.gtf.gz"),
    "(\S+).gtf.gz")

TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(glob.glob("*.bam"),
                                                                         "(\S+).bam")

REPLICATE = PipelineTracks.Aggregate(TRACKS, labels=("replicate", ))
TIME = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))

EXPRESSION = PARAMS['externals_expression_tables']
DIFF_DIR = PARAMS['externals_diff_directory']
COUNTS = PARAMS['externals_counts']
CLUSTERS = PARAMS['externals_clustering']
CONSENSUS = PARAMS['externals_consensus']
EIGENS = PARAMS['externals_eigengenes']
COMBINED = PARAMS['externals_combined']
###################################################################
###################################################################
###################################################################


def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################


@follows(connect)
@transform("/ifs/projects/proj036/jethro_bamfiles/deseq.dir/*_Jethro.tsv",
           regex("/ifs/projects/proj036/jethro_bamfiles/deseq.dir/FovsGC_DESeq2_Jethro.tsv"),
           r"/ifs/projects/proj036/jethro_bamfiles/deseq.dir/FovsGC_DESeq2_Jethro.load")
def loadFovsGC(infile, outfile):
    P.load(infile, outfiles)


@follows(loadFovsGC,
         mkdir("consensus_cluster.dir"))
@transform("%s/*-consensus.tsv" % CONSENSUS,
           regex("%s/(.+)-(.+)-consensus.tsv"),
           r"consensus_cluster.dir/\1-\2-consensus.load")
def loadConsenusClustering(infile, outfile):
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@follows(loadFovsGC,
         mkdir("report_summary.dir"))
@collate("%s/*-consensus.tsv" % CONSENSUS,
         regex("%s/(.+)-(.+)-consensus.tsv"),
         r"report_summary.dir/clustering_summary")
def clusteringSummary(infiles, outfile):
    '''
    Generate a summary table of the consensus clustering
    across all conditions
    '''

    PipelineProject036.clusterSummary(infiles, outfile)
         
###################################################################
###################################################################
###################################################################


@follows(mkdir("expression.dir"))
@transform("%s/*-merged.tsv.gz" % COMBINED,
           regex("%s/(.+)-merged.tsv.gz" % COMBINED),
           r"expression.dir/\1-vst.tsv.gz")
def normaliseExpressionTable(infile, outfile):
    '''
    Library size normalise and perform variance stabilising
    transformation a la DESeq
    '''

    time_agg = TIME.__dict__['track2groups'].keys()
    time_points = [int(str(x).split("-")[1]) for x in time_agg]
    time_points = set(time_points)
    time_points = list(time_points)
    time_points.sort()
    time_points = [str(x) for x in time_points]
    rep_agg = REPLICATE.__dict__['track2groups'].keys()
    replicates = [str(x).split("-")[2] for x in rep_agg]

    time_points = ",".join(time_points)
    replicates = ",".join(replicates)

    statement = '''
    python %(scriptsdir)s/expression2expression.py
    --task=deseq
    --log=%(outfile)s.log
    --replicates=%(replicates)s
    --time=%(time_points)s
    %(infile)s
    | gzip > %(outfile)s'''

    P.run()


@follows(normaliseExpressionTable)
@transform(normaliseExpressionTable,
           suffix("-vst.tsv.gz"),
           "-average_expression.tsv.gz")
def averageMergedExpression(infile, outfile):
    '''
    Average expression over genes and lncRNAs for each
    condition.

    '''

    statement = '''
    python %(scriptsdir)s/expression2expression.py
    --task=average_expression
    --log=%(outfile)s.log
    %(infile)s
    | gzip > %(outfile)s'''
    P.run()


##################################################################
###################################################################
###################################################################


@follows(mkdir("diff_timepoints.dir"),
         mkdir("diff_condition.dir"),
         mkdir("images.dir"))
@transform("%s/*.tsv" % DIFF_DIR,
           regex(r"%s/(.+)_(.+)_(.+)-time.tsv" % DIFF_DIR),
           r"images.dir/\1_\3-vsFo_GC.png")
def compareFovsGC(infile, outfile):
    '''
    Compare results from time point differential expression
    analysis to Fo -> GC differential expression analysis.
    Generate scatter plots for intersection of differentially
    expressed genes.
    '''
    # pass a data frame to cluster_runnable decorated
    # function rather than creating multiple sql
    # queries to avoid locking issues with sql database

    dbh = connect()
    fo_gc = sql.read_sql("SELECT * FROM FovsGC_deseq;",
                         dbh,
                         index_col="gene_id")

    image_dir = outfile.split("/")[0]
    PipelineProject036.compareFovsGC(infile,
                                     fo_gc,
                                     image_dir,
                                     submit=True)

##################################################################
###################################################################
###################################################################


@follows(compareFovsGC,
         mkdir("FovsGC_compare.dir"))
@transform("%s/*.tsv" % DIFF_DIR,
           regex("%s/(.+)_(.+)_(.+)-time.tsv" % DIFF_DIR),
           r"FovsGC_compare.dir/\1_\3_Fo-GC-intersect.tsv")
def coreOverlapFovsGC(infile, outfile):
    '''
    Take all intersecting genes over all conditions for
    time points 0-48 hours, intersect with FovsGC genes
    '''

    dbh = connect()
    fo_gc = sql.read_sql('''SELECT * FROM FovsGC_deseq;''',
                         dbh,
                         index_col="gene_id")

    PipelineProject036.coreOverlapFoVsGC(infile=infile,
                                         fo_gc=fo_gc,
                                         submit=True)

###################################################################
###################################################################
###################################################################


@follows(coreOverlapFovsGC)
@collate(coreOverlapFovsGC,
         regex("FovsGC_compare.dir/(.+)_(.+)_Fo-GC-intersect.tsv"),
         r"FovsGC_compare.dir/\1-timepoint_intersection_genes.tsv")
def defineTimepointIntersection(infiles, outfile):
    '''
    Intersect all time-point Fo-GC intersections across
    all time points and conditions.
    Output core genes only.
    '''

    n_times = [1, 3, 6, 12, 24, 48]
    PipelineProject036.getTimepointIntersections(infiles=infiles,
                                                 n_times=n_times,
                                                 outfile=outfile,
                                                 submit=True)


@follows(defineTimepointIntersection)
@collate("FovsGC_compare.dir/*.tsv",
         regex("FovsGC_compare.dir/(.+)-timepoint_intersection_(.+).tsv"),
         r"images.dir/Fo_GC-timepoint_intersection_\2.png")
def plotTimepointIntersection(infiles, outfile):
    '''
    Plot intersection of all core genes for all conditions
    for genes and lncRNAs separately
    '''

    PipelineProject036.plotTimepointIntersection(infiles,
                                                 outfile,
                                                 submit=True)


@transform("time_DE_union.dir/*-summarised-expression.tsv",
           regex("time_DE_union.dir/(.+)-summarised-expression.tsv"),
           r"time_DE_union.dir/\1-CORE-cross_correlation.tsv")
def timeUnionCrossCorrelation(infile, outfile):
    '''
    Intersect the Fo -> GC DE genes with timepoint union
    DE genes, then calculate their cross-correlation
    '''

    job_memory = "8G"

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/CrossCorrelation.py
    --log=%(outfile)s.log
    %(infile)s
    > %(outfile)s'''

    P.run()


@follows(timeUnionCrossCorrelation)
@transform("time_DE_union.dir/*-timepoint_union_genes.tsv",
           regex("time_DE_union.dir/(.+)-timepoint_union_genes.tsv"),
           r"time_DE_union.dir/\1-fo_gc-intersection.tsv")
def defineUnionIntersectGenes(infile, outfile):
    '''
    Intersect the union of temporally differentially expressed
    genes with the Fo -> GC DE genes for each condition
    '''

    dbh = connect()
    fo_gc = sql.read_sql('''SELECT * FROM FovsGC_deseq;''',
                         dbh,
                         index_col="gene_id")

    PipelineProject036.unionOverlapFoVsGC(infile=infile,
                                          fo_gc=fo_gc,
                                          outfile=outfile,
                                          submit=True)



@follows(defineUnionIntersectGenes)
@collate(defineUnionIntersectGenes,
         regex("time_DE_union.dir/(.+)-fo_gc-intersection.tsv"),
         r"time_DE_union.dir/Union-core_genes.tsv")
def defineUnionCoreGenes(infiles, outfile):
    '''
    Intersect each of the union of temporally DE
    gene sets that intersect Fo -> GC into a
    single set of core genes
    '''

    PipelineProject036.getCoreUnionGenes(infiles=infiles,
                                         outfile=outfile,
                                         submit=True)


@follows(timeUnionCrossCorrelation,
         defineUnionCoreGenes)
@transform("tpm.dir/*.tpm",
           regex("tpm.dir/(.+).tpm"),
           add_inputs("time_DE_union.dir/Union-core_genes.tsv"),
           r"time_DE_union.dir/\1-lncRNA-Union_core-correlation.tsv")
def correlateUnionCoreGenesLncRNAs(infiles, outfile):
    '''
    Cross-correlation of core-defined protein-coding genes
    with all lncRNAs
    '''

    job_memory = "8G"

    infiles = ",".join(infiles)
    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/lnc2correlate.py
    --log=%(outfile)s.log
    --method=cross-correlation
    --task=correlate-lncs-genes
    %(infiles)s
    > %(outfile)s
    '''

    P.run()


@follows(coreOverlapFovsGC)
@collate(coreOverlapFovsGC,
         regex("(.+)_(.+)-intersect.tsv"),
         r"FovsGC_compare.dir/\2-core_genes.tsv")
def defineCoreGenes(infiles, outfile):
    '''
    Intersect all time-point Fo-GC intersections across
    all time points and conditions.
    Output core genes only.
    '''

    n_times = [1, 3, 6, 12, 24, 48]
    PipelineProject036.getCoreGenes(infiles=infiles,
                                    n_times=n_times,
                                    outfile=outfile,
                                    submit=True)


@follows(defineCoreGenes)
@transform(defineCoreGenes,
           suffix(".tsv"),
           ".load")
def loadCoreGenes(infile, outfile):
    P.load(infile, outfile)


@follows(loadCoreGenes)
@transform("FovsGC_compare.dir/Fo-GC-core_lncRNAs.tsv",
           suffix(".tsv"),
           ".load")
def loadCorelncRNAs(infile, outfile):
    P.load(infile, outfile)
###################################################################
###################################################################
###################################################################


@follows(defineTimepointIntersection, defineCoreGenes)
@collate("FovsGC_compare.dir/*.tsv",
         regex("FovsGC_compare.dir/(.+)-timepoint_intersection_(.+).tsv"),
         add_inputs(r"FovsGC_compare.dir/*-timepoint_intersection_\2.tsv"),
         r"FovsGC_compare.dir/\1-specific_\2.tsv")
def defineConditionGenes(infiles, outfile):
    '''
    Intersect all time points for each condition.
    '''

    header = infiles[0][0].split("/")[-1].split("-")[0]
    list_of_files = [f for f in infiles[0] if not re.search(header, f)]
    list_of_files.append(infiles[0][0])
    PipelineProject036.getConditionGenes(list_of_files=list_of_files,
                                         reference=header,
                                         outfile=outfile,
                                         submit=True)


@follows(defineConditionGenes)
@transform(defineConditionGenes,
           suffix(".tsv"),
           ".load")
def loadConditionGenes(infile, outfile):
    P.load(infile, outfile)
###################################################################
###################################################################
###################################################################


@follows(defineConditionGenes,
         averageMergedExpression)
@transform(averageMergedExpression,
           regex("expression.dir/(.+)-average_expression.tsv.gz"),
           add_inputs([r"FovsGC_compare.dir/Fo-GC-core_genes.tsv",
                       r"FovsGC_compare.dir/Fo-GC-core_lncRNAs.tsv"]),
           r"FovsGC_compare.dir/\1-core_correlations.tsv")
def correlateCoreGenesLncRNAs(infiles, outfile):
    '''
    Correlate expression of core genes and lncRNAs for each condition
    '''

    expr_file = infiles[0]
    gene_file = infiles[1][0]
    lnc_file = infiles[1][1]
    PipelineProject036.correlateGeneLncRNA(gene_file=gene_file,
                                           lnc_file=lnc_file,
                                           expression_file=expr_file,
                                           outfile=outfile,
                                           submit=True)
###################################################################
###################################################################
###################################################################


@follows(defineConditionGenes,
         averageMergedExpression,
         correlateCoreGenesLncRNAs)
@transform(averageMergedExpression,
           regex("expression.dir/(.+)-average_expression.tsv.gz"),
           add_inputs([r"FovsGC_compare.dir/\1-specific_genes.tsv",
                       r"FovsGC_compare.dir/\1-specific_lncRNAs.tsv"]),
           r"FovsGC_compare.dir/\1-specific_correlations.tsv")
def conditionCorrelationGenesLncRNAs(infiles, outfile):
    '''
    Condition specific expression correlation of genes and lncRNAs
    '''

    expr_file = infiles[0]
    gene_file = infiles[1][0]
    lnc_file = infiles[1][1]
    PipelineProject036.correlateGeneLncRNA(gene_file=gene_file,
                                           lnc_file=lnc_file,
                                           expression_file=expr_file,
                                           outfile=outfile,
                                           submit=True)


@follows(conditionCorrelationGenesLncRNAs)
@transform(conditionCorrelationGenesLncRNAs,
           regex("FovsGC_compare.dir/(.+)-specific_correlations.tsv"),
           r"images.dir/\1-specific_correlations.png")
def heatmapConditionCorrelationGenesLncRNAs(infile, outfile):
    '''
    plot heatmap of cross-correlation between condition specific
    genes and lncRNAs
    '''

    PipelineProject036.plotConditionHeatmap(infile,
                                            outfile,
                                            submit=True)

###################################################################
###################################################################
###################################################################


@follows(defineConditionGenes,
         defineCoreGenes,
         mkdir("go_enrichment.dir/"))
@collate(averageMergedExpression,
         regex("expression.dir/(.+)-average_expression.tsv.gz"),
         r"go_enrichment.dir/Core_background_go.tsv")
def generateCoreBackground(infile, outfile):
    '''
    Generate list of background genes for GO enrichment of core genes
    Background genes are all differentially expressed genes between
    follicular and germinal centre B cells
    '''

    dbh = connect()
    fo_gc = sql.read_sql("SELECT * FROM FovsGC_deseq;",
                         dbh,
                         index_col="gene_id")
    temp = P.getTempFilename(shared=True)
    fo_gc.to_csv(temp, sep="\t")
    PipelineProject036.generateBackground(infiles=temp,
                                          outfile=outfile,
                                          submit=True)
###################################################################
###################################################################
###################################################################


@follows(defineCoreGenes,
         generateCoreBackground)
@transform(defineCoreGenes,
           regex("FovsGC_compare.dir/Fo-GC-core_genes.tsv"),
           add_inputs(r"go_enrichment.dir/Core_background_go.tsv"),
           r"go_enrichment.dir/Fo-GC-core_genes-go.tsv")
def geneOntologyCoreGenes(infiles, outfile):
    '''
    GO enrichment for gene_list of core overlap genes
    '''

    gene_file = infiles[0]
    bg_file = infiles[1]
    PipelineProject036.goEnrichment(gene_set=gene_file,
                                    bg_set=bg_file,
                                    genome=PARAMS['genome'],
                                    db_ids="ensGene",
                                    outfile=outfile,
                                    database=PARAMS['go_database'],
                                    submit=True)

###################################################################
###################################################################
###################################################################


@follows(generateCoreBackground,
         defineConditionGenes)
@transform(averageMergedExpression,
           regex("expression.dir/(.+)-average_expression.tsv.gz"),
           r"go_enrichment.dir/\1-specific_background_go.tsv")
def generateSpecificBackground(infile, outfile):
    '''
    Generate background geneset for condition-specific GO enrichment
    '''
    dbh = connect()
    fo_gc = sql.read_sql("SELECT * FROM FovsGC_deseq;",
                         dbh,
                         index_col="gene_id")
    temp = P.getTempFilename(shared=True)
    fo_gc.to_csv(temp, sep="\t")

    PipelineProject036.generateBackground(infiles=temp,
                                          outfile=outfile,
                                          submit=True)
###################################################################
###################################################################
###################################################################


@follows(defineConditionGenes,
         generateSpecificBackground)
@transform(defineConditionGenes,
           regex("FovsGC_compare.dir/(.+)-specific_genes.tsv"),
           add_inputs(r"go_enrichment.dir/\1-specific_background_go.tsv"),
           r"go_enrichment.dir/\1-specfic_genes_go.tsv")
def geneOntologyConditionGenes(infiles, outfile):
    '''
    GO enrichment on condition specific genes overlapping
    with Fo -> GC analysis
    '''

    gene_file = infiles[0]
    bg_file = infiles[1]

    PipelineProject036.goEnrichment(gene_set=gene_file,
                                    bg_set=bg_file,
                                    genome=PARAMS['genome'],
                                    db_ids="ensGene",
                                    outfile=outfile,
                                    database=PARAMS['go_database'],
                                    submit=True)
###################################################################
###################################################################
###################################################################


@follows(mkdir("go_enrichment.dir/"))
@subdivide("%s/*-consensus.tsv.gz" % CONSENSUS,
           regex("%s/(.+)-refcoding-consensus.tsv.gz" % CONSENSUS),
           r"go_enrichment.dir/\1_*")
def clusterGeneOntology(infile, outfiles):
    '''
    GO enrichment on each cluster from the consensus clustering
    '''

    header = infile.split("/")[-1].split("-")[0]

    PipelineProject036.clusterGOEnrichment(cluster_file=infile,
                                           genome=PARAMS['genome'],
                                           db_ids="ensGene",
                                           label=header,
                                           out_dir="go_enrichment.dir",
                                           submit=True)


@follows(clusterGeneOntology)
@transform("go_enrichment.dir/*GO.tsv",
           regex("go_enrichment.dir/(.+)_(.+)GO.tsv"),
           r"go_enrichment.dir/\1-\2GO.load")
def loadClusterGeneOntology(infile, outfile):
    P.load(infile, outfile)


@follows(loadClusterGeneOntology)
@transform("go_enrichment.dir/*GO.tsv",
           regex("go_enrichment.dir/(.+)_(.+)GO.tsv"),
           r"go_enrichment.dir/\1-\2-enrich.tsv")
def topClusterGeneOntology(infile, outfile):
    '''
    Get fold enrichments and output top GO
    terms
    '''

    condition = infile.split("/")[-1].split("_")[0]
    expr_dir = PARAMS['externals_expression_tables']
    expression_file = "%s/%s-refcoding-filtered-vst.tsv" % (expr_dir,
                                                            condition)
    con_dir = PARAMS['externals_consensus']
    cluster_file = "%s/%s-refcoding-consensus.tsv.gz" % (con_dir,
                                                         condition)

    PipelineProject036.topGO(go_file=infile,
                             expression_file=expression_file,
                             cluster_file=cluster_file,
                             outfile=outfile,
                             submit=True)


@follows(topClusterGeneOntology)
@transform(topClusterGeneOntology,
           suffix(".tsv"),
           ".load")
def loadTopClusterGO(infile, outfile):
    P.load(infile, outfile)


@follows(loadClusterGeneOntology)
@collate("go_enrichment.dir/*GO.tsv",
         regex("go_enrichment.dir/(.+)_(.+)GO.tsv"),
         r"go_enrichment.dir/\1-GO_summary.tsv")
def summariseClusterGO(infiles, outfile):
    '''
    Summarise GO enrichments over all clusters.
    '''

    condition = infiles[0].split("/")[-1].split("_")[0]
    expr_dir = PARAMS['externals_expression_tables']
    expression_file = "%s/%s-refcoding-filtered-vst.tsv" % (expr_dir,
                                                            condition)
    con_dir = PARAMS['externals_consensus']
    cluster_file = "%s/%s-refcoding-consensus.tsv.gz" % (con_dir,
                                                         condition)

    PipelineProject036.summariseGO(list_of_go=infiles,
                                   expression_file=expression_file,
                                   cluster_file=cluster_file,
                                   outfile=outfile,
                                   submit=True)


@follows(summariseClusterGO)
@transform(summariseClusterGO,
           suffix(".tsv"),
           ".load")
def loadSummariseClusterGO(infile, outfile):
    P.load(infile, outfile)
           
###################################################################
###################################################################
###################################################################


@follows(mkdir("eigengenes.dir"))
@transform("%s/*-cluster_eigengenes.tsv.gz" % EIGENS,
           regex("%s/(.+)-PCA-refcoding-cluster_eigengenes.tsv.gz" % EIGENS),
           add_inputs(r"%s/\1-PCA-lncRNA_merged-cluster_eigengenes.tsv.gz" % EIGENS),
           r"eigengenes.dir/\1-correlated-eigengenes.tsv.gz")
def correlateEigengenes(infiles, outfile):
    '''
    Correlate eigengene profiles of protein-coding genes and lncRNAs
    '''

    infiles = ",".join(infiles)
    
    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=correlate-eigens
    %(infiles)s | gzip > %(outfile)s'''

    P.run()


@follows(correlateEigengenes)
@transform("%s/*-cluster_eigengenes.tsv.gz" % EIGENS,
           regex("%s/(.+)-PCA-(.+)-cluster_eigengenes.tsv.gz" % EIGENS),
           add_inputs([r"expression.dir/\1-average_expression.tsv.gz",
                       r"consensus_cluster.dir/\1-\2-consensus.tsv.gz"]),
           r"images.dir/\1-\2_eigengene-heatmaps.tsv")
def plotClusterHeatMaps(infiles, outfile):
    '''
    Generate heatmap of cluster expression with eigengene expression.
    Use sentinel file for pipeline.
    '''

    eigen_file = infiles[0]
    express = infiles[1][0]
    consensus = infiles[1][1]

    PipelineProject036.plotClusterHeatmaps(eigengenes=eigen_file,
                                           expression=express,
                                           clusters=consensus,
                                           image_dir="images.dir",
                                           submit=True)

    P.touch(outfile)


@follows(correlateEigengenes)
@transform(correlateEigengenes,
           regex("eigengenes.dir/(.+)-correlated-eigengenes.tsv.gz"),
           r"images.dir/\1-eigengene-correlation-heatmap.png")
def plotCorrelateEigengenes(infile, outfile):
    '''
    Generate a heatmap of eigengene correlations
    '''
    
    PipelineProject036.plotEigenHeatmaps(eigengenes=infile,
                                         image_dir="images.dir",
                                         submit=True)
###################################################################
###################################################################
###################################################################


@follows(correlateEigengenes)
@transform(correlateEigengenes,
           regex("eigengenes.dir/(.+)-correlated-eigengenes.tsv.gz"),
           r"eigengenes.dir/\1-eigengene-correlation.tsv.gz")
def matchEigengenes(infile, outfile):
    '''
    Correlate lncRNA eigengenes and protein-coding eigengenes.
    '''

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=filter-correlations
    --cor-threshold=%(lncrna_correlation_threshold)s
    %(infile)s | gzip > %(outfile)s '''

    P.run()

###################################################################
###################################################################
###################################################################


@follows(matchEigengenes,
         averageMergedExpression,
         mkdir("lncRNA_classification.dir"))
@transform(matchEigengenes,
           regex("eigengenes.dir/(.+)-eigengene-correlation.tsv.gz"),
           add_inputs([r"expression.dir/\1-average_expression.tsv.gz",
                       r"%s/\1-refcoding-consensus.tsv" % CONSENSUS,
                       r"%s/\1-lncRNA_merged-consensus.tsv" % CONSENSUS]),
           r"lncRNA_classification.dir/\1-lncRNA_merged-gene_correlations.tsv.gz")
def correlateLncRNAsWithClusterGenes(infiles, outfile):
    '''
    Correlate expression profiles of lncRNAs with genes
    within clusters also having correlation with cluster
    eigengene.
    Output is a flat file of gene:lncRNA:correlation:cluster
    '''

    input_files = [infiles[0],
                   infiles[1][0],
                   infiles[1][1],
                   infiles[1][2]]

    input_files = ",".join(input_files)

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=correlate-cluster-genes
    %(input_files)s | gzip > %(outfile)s'''

    P.run()


@follows(correlateUnionCoreGenesLncRNAs)
@transform(correlateUnionCoreGenesLncRNAs,
           regex("time_DE_union.dir/(.+)-lncRNA-Union_core-correlation.tsv"),
           r"lncRNA_classification.dir/\1-CORE_Union_lncRNA-gene_correlations.tsv.gz")
def flattenCoreUnionCorrelations(infile, outfile):
    '''
    Select highly correlated lncRNA:protein-coding gene pairs
    and flatten, one row per correlation
    '''

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=filter-correlations
    %(infile)s | gzip > %(outfile)s
    '''

    P.run()


@follows(defineConditionGenes,
         defineCoreGenes,
         correlateLncRNAsWithClusterGenes,
         conditionCorrelationGenesLncRNAs)
@transform(conditionCorrelationGenesLncRNAs,
           regex("FovsGC_compare.dir/(.+)-specific_correlations.tsv"),
           r"lncRNA_classification.dir/\1-specific_lncRNA_merged-gene_correlations.tsv.gz")
def flattenConditionCorrelations(infile, outfile):
    '''
    flatten correlations between condition-specific genes and lncRNAs
    '''

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=filter-correlations
    %(infile)s | gzip > %(outfile)s'''

    P.run()


@follows(defineConditionGenes,
         defineCoreGenes,
         correlateLncRNAsWithClusterGenes,
         conditionCorrelationGenesLncRNAs,
         flattenConditionCorrelations)
@transform(correlateCoreGenesLncRNAs,
           regex("FovsGC_compare.dir/(.+)-core_correlations.tsv"),
           r"lncRNA_classification.dir/\1-core_lncRNA_merged-gene_correlations.tsv.gz")
def flattenCoreCorrelations(infile, outfile):
    '''
    flatten correlations between core genes and lncRNAs
    '''

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=filter-correlations
    %(infile)s | gzip > %(outfile)s'''

    P.run()


@jobs_limit(2)
@follows(correlateLncRNAsWithClusterGenes,
         flattenCoreCorrelations,
         flattenConditionCorrelations,
         flattenCoreUnionCorrelations)
@transform([correlateLncRNAsWithClusterGenes,
            flattenCoreCorrelations,
            flattenConditionCorrelations,
            flattenCoreUnionCorrelations],
           suffix("-gene_correlations.tsv.gz"),
           add_inputs("lncRNA_merged.gtf.gz"),
           "-annotated.tsv")
def annotateCorrelatedLncRNAs(infiles, outfile):
    '''
    Add gene and lncRNA annotations to correlation
    list - lncRNA co-ordinates and gene symbols
    '''

    lnc_gtf = infiles[1]
    infile = infiles[0]

    PipelineProject036.annotateGeneList(infile=infile,
                                        lnc_gtf=lnc_gtf,
                                        outfile=outfile,
                                        submit=True,
                                        job_options="-l mem_free=5G")


@follows(annotateCorrelatedLncRNAs)
@transform(annotateCorrelatedLncRNAs,
           suffix(".tsv"),
           ".load")
def loadAnnotatedCorrelations(infile, outfile):
    P.load(infile, outfile)


@follows(loadAnnotatedCorrelations)
@transform(annotateCorrelatedLncRNAs,
           suffix("-annotated.tsv"),
           "-gene_counts.tsv")
def countLncsPerGene(infile, outfile):
    '''
    Count the number of lncRNAs correlated with
    each protein-coding gene
    '''

    PipelineProject036.lncsPerGene(infile,
                                   outfile,
                                   submit=True)


@follows(annotateCorrelatedLncRNAs,
         countLncsPerGene)
@transform(annotateCorrelatedLncRNAs,
           suffix("-annotated.tsv"),
           "-lncRNA_counts.tsv")
def countGenesPerLnc(infile, outfile):
    '''
    Count the number of correlated genes per lncRNA,
    subset by lncRNA class
    '''

    PipelineProject036.genesPerLnc(infile,
                                   outfile,
                                   submit=True,
                                   job_options="-l mem_free=7G")


@follows(countLncsPerGene)
@transform(countLncsPerGene,
           suffix(".tsv"),
           ".load")
def loadCountsPerGene(infile, outfile):
    P.load(infile, outfile)


@follows(countGenesPerLnc,
         loadCountsPerGene)
@transform(countGenesPerLnc,
           suffix(".tsv"),
           ".load")
def loadCountsPerLnc(infile, outfile):
    P.load(infile, outfile)


@follows(loadCountsPerGene,
         loadCountsPerLnc)
@transform(countLncsPerGene,
           regex("lncRNA_classification.dir/(.+)-(.+)-gene_counts.tsv"),
           r"images.dir/\1-\2-histogram-gene_counts.png")
def plotLncsPerGene(infile, outfile):
    '''
    Plot histogram of lncRNAs per gene
    '''

    PipelineProject036.plotGeneCounts(infile,
                                      outfile,
                                      submit=True,
                                      job_options="-l mem_free=4G")


@follows(loadCountsPerLnc,
         plotLncsPerGene)
@transform(countGenesPerLnc,
           regex("lncRNA_classification.dir/(.+)-(.+)-lncRNA_counts.tsv"),
           r"images.dir/\1_\2-density-lncRNA_counts.png")
def plotCountsPerLnc(infile, outfile):
    '''
    Plot histogram of genes per lncRNA, subset by lncRNA class
    '''

    PipelineProject036.plotLncCounts(infile,
                                     outfile,
                                     submit=True,
                                     job_options="-l mem_free=8G")

###################################################################
###################################################################
###################################################################


@follows(correlateLncRNAsWithClusterGenes,
         flattenConditionCorrelations,
         flattenCoreCorrelations,
         plotCountsPerLnc)
@transform([correlateLncRNAsWithClusterGenes,
            flattenCoreCorrelations,
            flattenConditionCorrelations,
            flattenCoreUnionCorrelations],
           regex("lncRNA_classification.dir/(.+)-gene_correlations.tsv.gz"),
           add_inputs(r"lncRNA_merged.gtf.gz"),
           r"lncRNA_classification.dir/\1-positive_lncRNAs.tsv.gz")
def classifyPosLncRNAs(infiles, outfile):
    '''
    Classify multi-exonic intergenic lncRNAs whose expression is positively
    correlated with protein-coding gene expression above a threshold.
    '''

    input_files = ",".join(infiles)
    job_options = "-l mem_free=6G"

    statement = '''python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=classify-lncrna
    --lncRNA-class=intergenic
    --cor-direction=positive
    --cor-threshold=%(lncrna_correlation_threshold)s
    %(input_files)s
    | gzip > %(outfile)s'''

    P.run()


@follows(correlateLncRNAsWithClusterGenes,
         flattenConditionCorrelations,
         flattenCoreCorrelations,
         flattenCoreUnionCorrelations,
         plotCountsPerLnc)
@transform([correlateLncRNAsWithClusterGenes,
            flattenCoreCorrelations,
            flattenConditionCorrelations,
            flattenCoreUnionCorrelations],
           regex("lncRNA_classification.dir/(.+)-gene_correlations.tsv.gz"),
           add_inputs(r"lncRNA_merged.gtf.gz"),
           r"lncRNA_classification.dir/\1-negative_lncRNAs.tsv.gz")
def classifyNegLncRNAs(infiles, outfile):
    '''
    Classify multi-exonic intergenic lncRNAs whose expression is negatively
    correlated with protein-coding gene expression above a threshold.
    '''

    input_files = ",".join(infiles)

    job_options = "-l mem_free=6G"

    statement = '''python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/lnc2correlate.py
    --log=%(outfile)s.log
    --task=classify-lncrna
    --lncRNA-class=intergenic
    --cor-direction=negative
    --cor-threshold=%(lncrna_correlation_threshold)s
    %(input_files)s
    | gzip > %(outfile)s'''

    P.run()


####################################################################
####################################################################
####################################################################


@follows(classifyPosLncRNAs,
         classifyNegLncRNAs,
         mkdir("correlations.dir"))
@collate("*.gtf.gz",
         regex("refcoding.gtf.gz"),
         add_inputs(r"lncRNA_merged.gtf.gz"),
         r"correlations.dir/gene_distances.tsv.gz")
def getGeneDistances(infiles, outfile):
    '''
    Returns a list of lncRNAs:gene distances
    '''

    genes = infiles[0][0]
    lncs = infiles[0][1]

    job_options = "-l mem_free=2G"

    statement = '''zcat %(lncs)s |
    cgat %(scriptsdir)s/gtf2table.py
    --gff-file=%(genes)s
    --filename-format=gtf
    --log=neighbours.log
    --counter=distance-genes
    | gzip > %(outfile)s'''

    P.run()
####################################################################
####################################################################
####################################################################


@follows(classifyPosLncRNAs,
         getGeneDistances)
@transform(classifyPosLncRNAs,
           regex("lncRNA_classification.dir/(.+)-positive_lncRNAs.tsv.gz"),
           r"correlations.dir/\1-correlation_pairs.tsv")
def getCorrelatePairs(infile, outfile):
    '''
    Select correlated gene pairs >= threshold for each lncRNA
    '''

    job_options = "-l mem_free=4G"

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/lnc2correlate.py
    --task=filter-pairs
    --cor-threshold=%(lncrna_correlation_threshold)s
    --log=%(outfile)s.log
    %(infile)s
    > %(outfile)s'''

    P.run()


@follows(classifyPosLncRNAs,
         getCorrelatePairs,
         getGeneDistances)
@transform(classifyPosLncRNAs,
           regex("lncRNA_classification.dir/(.+)-positive_lncRNAs.tsv.gz"),
           r"correlations.dir/\1-maxCor_pairs.tsv")
def getMaxCorrelatePairs(infile, outfile):
    '''
    Select most higly correlated gene pair for each lncRNA
    '''

    job_options = "-l mem_free=4G"

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/lnc2correlate.py
    --task=max-correlation
    --log=%(outfile)s.log
    %(infile)s
    > %(outfile)s'''

    P.run()

####################################################################
####################################################################
####################################################################


@follows(getGeneDistances,
         getCorrelatePairs)
@transform(getCorrelatePairs,
           regex("correlations.dir/(.+)-(.+)-correlation_pairs.tsv"),
           add_inputs([r"correlations.dir/gene_distances.tsv.gz",
                       r"expression.dir/\1-average_expression.tsv.gz"]),
           r"correlations.dir/\1-\2-proximal_pairs.tsv")
def getProximalPairs(infiles, outfile):
    '''
    Get a list of lncRNAs with closet protein-coding gene
    and their correlations
    '''

    cor_pairs = infiles[0]
    distance_file = infiles[1][0]
    expr_file = infiles[1][1]

    PipelineProject036.correlateProximalPairs(distances=distance_file,
                                              pairs_file=cor_pairs,
                                              expression=expr_file,
                                              outfile=outfile,
                                              submit=True)

@follows(getProximalPairs)
@transform(getCorrelatePairs,
           regex("correlations.dir/(.+)-(.+)-correlation_pairs.tsv"),
           add_inputs(r"lncRNA_classification.dir/\1-\2-negative_lncRNAs.tsv.gz"),
           r"correlations.dir/\1-\2-anticorrelation_pairs.tsv")
def getAntiCorrelations(infiles, outfile):
    '''
    Select the most strongly anti-correlated protein-coding
    gene for each lncRNA
    '''

    cor_pairs = infiles[0]
    expr_file = infiles[1]

    PipelineProject036.antiCorrelatePairs(pairs_file=cor_pairs,
                                          expression=expr_file,
                                          outfile=outfile,
                                          threshold=-float(PARAMS['lncrna_correlation_threshold']),
                                          submit=True,
                                          job_options="-l mem_free=32G")

####################################################################
####################################################################
####################################################################


@follows(getProximalPairs,
         getCorrelatePairs)
@transform(getProximalPairs,
           regex("correlations.dir/(.+)-(.+)-proximal_pairs.tsv"),
           add_inputs([r"expression.dir/\1-average_expression.tsv.gz",
                       "lncRNA_merged.gtf.gz"]),
           r"correlations.dir/\1-\2-random_pairs.tsv")
def getRandomPairs(infiles, outfile):
    '''
    Randomly select lncRNA:gene pairs to create a
    null distribution of correlations
    '''

    cor_pairs = infiles[0]
    expr_file = infiles[1][0]
    ref_gtf = infiles[1][1]

    PipelineProject036.correlateRandomPairs(pairs_file=cor_pairs,
                                            expression=expr_file,
                                            ref_gtf=ref_gtf,
                                            outfile=outfile,
                                            seed=int(PARAMS['seed']),
                                            submit=True)


@follows(getRandomPairs,
         getCorrelatePairs,
         getProximalPairs,
         getAntiCorrelations,
         plotCountsPerLnc)
@collate("correlations.dir/*.tsv",
         regex("correlations.dir/(.+)-(.+)-(.+)_pairs.tsv"),
         r"images.dir/\1-\2-correlation_distribution.png")
def plotCorrelations(infiles, outfile):
    '''
    Plot distributions of correlations of lncRNA:gene pairs
    '''
    cor_file = [x for x in infiles if re.search("correlation", x)][0]
    prox_file = [p for p in infiles if re.search("proximal", p)][0]
    rand_file = [r for r in infiles if re.search("random", r)][0]
    anti_file = [a for a in infiles if re.search("anti", a)][0]

    PipelineProject036.plotCorrelations(cor_file=cor_file,
                                        rand_file=rand_file,
                                        prox_file=prox_file,
                                        anti_file=anti_file,
                                        outfile=outfile,
                                        submit=True)

####################################################################
####################################################################
####################################################################


@follows(getCorrelatePairs,
         getRandomPairs,
         getProximalPairs,
         getAntiCorrelations,
         mkdir("merged_gtf.dir"))
@collate("*gtf.gz",
         regex("(.+).gtf.gz"),
         r"merged_gtf.dir/merged.gtf.gz")
def mergeGTFs(infiles, outfile):
    '''
    merge protein-coding gene and lncRNA GTFs together
    '''

    files = " ".join(infiles)
    statement = ''' zcat %(files)s | gzip > %(outfile)s'''

    P.run()


@follows(getCorrelatePairs,
         getProximalPairs,
         getRandomPairs,
         getAntiCorrelations,
         mergeGTFs,
         mkdir("mre.dir"))
@transform("correlations.dir/*_pairs.tsv",
           regex("correlations.dir/(.+)-(.+)_pairs.tsv"),
           add_inputs(r"merged_gtf.dir/merged.gtf.gz"),
           r"lncRNA_classification.dir/\1-\2_pairs.gtf")
def pairs2GTF(infiles, outfile):
    '''
    Return GTF file of positively correlated intergenic lncRNAs
    and protein-coding genes
    '''

    infile = infiles[0]
    merge_gtf = infiles[1]
    in_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    genes_list = set(in_frame.index)
    lncs_list = set(in_frame['lncRNA_id'].values)
    # I HATE IMPLICIT IN-PLACE OPERATIONS!!!!
    lncs_list.update(genes_list)

    PipelineProject036.list2GTF(list_of_ids=lncs_list,
                                gtf_file=merge_gtf,
                                out_gtf=outfile,
                                submit=True,
                                job_options="-l mem_free=24G")

###################################################################
###################################################################
###################################################################


@follows(mkdir("fasta.dir"),
         pairs2GTF)
@transform(pairs2GTF,
           regex("lncRNA_classification.dir/(.+)-(.+)_pairs.gtf"),
           r"fasta.dir/\1-\2-separate_exons.fa.gz")
def gtf2Fasta(infile, outfile):
    '''
    Take GTF file and pull out fasta format sequence
    from the appropriate genome directory for each exon,
    merge later
    '''
    
    genome = os.path.join(PARAMS['genome_dir'],
                          PARAMS['genome'])

    job_options = "-l mem_free=4G"
        
    statement = '''
    cat %(infile)s |
    python %(scriptsdir)s/gff2fasta.py
    --is-gtf
    --genome-file=%(genome)s
    --log=%(outfile)s.log
    | gzip  > %(outfile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(gtf2Fasta)
@transform(gtf2Fasta,
           suffix("-separate_exons.fa.gz"),
           "-spliced.fa")
def makeSplicedTranscripts(infile, outfile):
    '''
    Concatenate sequences to form contiguous spliced
    transcript sequences
    '''

    PipelineProject036.makeSplicedFasta(infile,
                                        outfile)

###################################################################
###################################################################
###################################################################


@follows(gtf2Fasta,
         makeSplicedTranscripts,
         mkdir("target_scan.dir"))
@files("miR_Family_Info.txt",
       "target_scan.dir/10090_miRNA.seed.tsv")
def miRNASeedFile(infile, outfile):
    '''
    Generate target scan miRNA seed sequence file - single species
    only - no alignment with other species.
    '''

    statement = '''cat %(infile)s | cut -f 1,2,3 > %(outfile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(gtf2Fasta, miRNASeedFile)
@transform("fasta.dir/*-spliced.fa",
           regex("fasta.dir/(.+)-(.+)-spliced.fa"),
           r"target_scan.dir/\1-\2-targetscan.tsv")
def targetFile(infile, outfile):
    '''
    Generate target sequence file from input fasta format.
    '''

    statement = '''cat %(infile)s | cut -d ' ' -f 1 | cut -d '>' -f 2
                   | awk '{if(NR%%2 != 0) {printf("%%s\\t%%s\\t", $0, %(ncbi)s)} else {printf("%%s\\n", $0)}}'
                   > %(outfile)s'''
    P.run()
###################################################################
###################################################################
###################################################################


@follows(miRNASeedFile,
         targetFile)
@transform(targetFile,
           regex("target_scan.dir/(.+)-(.+)-targetscan.tsv"),
           add_inputs("target_scan.dir/10090_miRNA.seed.tsv"),
           r"mre.dir/\1-\2-targetScan_results.tsv")
def runTargetScan(infile, outfile):
    '''
    Predict miRNA recognition elements using TargetScan. Do not
    use aligned sequences, match for seed sequences in miRNA from
    input species.  This allows identification of conserved MREs
    at a later time point and comparison to miRanda output.
    The TargetScan script `targetscan_60.pl` needs to be in the
    working directory.
    '''
    
    miRNA_file = infile[1]
    target_file = infile[0]

    PipelineProject036.targetScanWrapper(miRNA_file=miRNA_file,
                                         target_file=target_file,
                                         outfile=outfile)

###################################################################
###################################################################
###################################################################


@follows(runTargetScan)
@transform(runTargetScan,
           regex("mre.dir/(.+)-(.+)-targetScan_results.tsv"),
           add_inputs(r"lncRNA_classification.dir/\1-\2_pairs.gtf"),
           r"mre.dir/\1-\2-targetScan-mre.gtf.gz")
def targetScanParse(infile, outfile):
    '''
    Parse resuls from targetScan into a gff of MREs.
    Requires original target gtf.
    '''

    lnc_file = infile[1]
    infile = infile[0]

    job_options = "-l mem_free=6G"

    statement = ''' cat %(infile)s | 
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/miranda2GTF.py
    --format=targetscan
    --GTF=%(lnc_file)s
    --log=%(outfile)s.log | gzip > %(outfile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(runTargetScan, targetScanParse)
@transform(targetScanParse,
           regex("mre.dir/(.+)-(.+)-mre.gtf.gz"),
           add_inputs("merged_gtf.dir/merged.gtf.gz"),
           r"mre.dir/\1-\2-mre.annotated.gtf.gz")
def annotateMreGTF(infile, outfile):
    '''
    Annotate MRE GTF file with target gene/transcript attributes

    NB: needs to wrapped in a script to make it compatible with
    working on the stream.  Fill also allow compression of output
    files.
    '''
    
    lnc_gtf = infile[1]
    mre_gtf = infile[0]

    job_memory = "16G"
    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/mre2mre.py
    --log=%(outfile)s.log
    --task=annotate
    --annotation-gtf-file=%(lnc_gtf)s
    %(mre_gtf)s | gzip > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################

###########################################################
# add in functions to stratify miRNA expression and filter
# on highest expression strata - look at Ana's latest paper
###########################################################


@follows(annotateMreGTF)
@transform(annotateMreGTF,
           regex("mre.dir/(.+)-(.+)-mre.annotated.gtf.gz"),
           r"mre.dir/\1-\2-mre.filtered.gtf.gz")
def filterMREs(infile, outfile):
    '''
    Filter out MREs based on a `.tsv` file of
    miRNA IDs.
    '''

    job_options = "-l mem_free=10G"

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/proj036/src/mre2mre.py
    --log=%(outfile)s.log
    --task=filter
    --method=list
    --filter-filename=%(mre_filter_file)s
    %(infile)s | gzip > %(outfile)s
    '''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(filterMREs)
@transform(filterMREs,
           suffix("filtered.gtf.gz"),
           "unique.gtf.gz")
def cleanMREs(infile, outfile):
    '''
    Sort MREs by gene_id and remove duplicates.
    Crop out multiple redundant MREs that share the same genome
    co-ordinates - often due to multiple very closely related miRNAs
    '''
        
    job_options = "-l mem_free=5G"

    statement = '''
    zcat %(infile)s |
    python %(scriptsdir)s/gtf2gtf.py
    --method=sort
    --sort-order=gene
    -v 0
    | python %(scriptsdir)s/gtf2gtf.py
    --method=remove-duplicates
    --duplicate-feature=gene
    --log=%(outfile)s.log | gzip > %(outfile)s '''
    
    P.run()
###################################################################
###################################################################
###################################################################


@follows(cleanMREs)
@transform(cleanMREs,
           suffix("unique.gtf.gz"),
           "sorted.gtf.gz")
def sortMREs(infile, outfile):
    '''
    Sort MREs into contig/position order
    '''

    job_options = "-l mem_free=4G"

    statement = '''zcat %(infile)s |
                   python %(scriptsdir)s/gtf2gtf.py
                   --method=sort
                   --sort-order=position
                   --log=%(outfile)s.log | gzip > %(outfile)s'''

    P.run()

###################################################################
###################################################################
###################################################################


@follows(sortMREs)
@transform(sortMREs,
           regex("mre.dir/(.+)-(.+)-targetScan-mre.sorted.gtf.gz"),
           add_inputs(r"lncRNA_classification.dir/\1-\2_pairs.gtf"),
           r"mre.dir/\1-\2-targetScan-mre_counts.tsv")
def countMREsOverLncs(infile, outfile):
    '''
    Count the number of non-redundant MREs overlapping
    lncRNA gene models in input
    '''

    mre_gtf = infile[0]
    lnc_gtf = infile[1]

    PipelineProject036.countMREsOverLncs(mre_gtf=mre_gtf,
                                         lnc_gtf=lnc_gtf,
                                         outfile=outfile,
                                         submit=True,
                                         job_options="-l mem_free=16G")


###################################################################
###################################################################
###################################################################


@follows(sortMREs)
@transform(sortMREs,
           regex("mre.dir/(.+)-(.+)-targetScan-mre.sorted.gtf.gz"),
           add_inputs([r"merged_gtf.dir/merged.gtf.gz",
                       r"correlations.dir/\1-\2_pairs.tsv"]),
           r"mre.dir/\1-\2_pairs-mre_shared.tsv")
def sharedMREs(infiles, outfile):
    '''
    Get shared MREs between highly correlated lncRNA:gene pairs
    output column descriptions:
    gene_id - ensembl gene id
    gene_shared_proportion - proportion of gene MREs shared with this lncRNA
    lncRNA_id - catalog lncRNA id
    lncRNA_shared_proportion - proportion of lncRNA MREs shared with this gene
    shared_miRNAs - miRNA ids with shared MREs between this gene and lncRNA
    total_shared - total number of different miRNAs with MREs shared between
    this gene and lncRNA
    '''

    mre_file = infiles[0]
    pairs_gtf = infiles[1][0]
    correlations = infiles[1][1]

    PipelineProject036.shareMREs(mre_file=mre_file,
                                 pairs_gtf=pairs_gtf,
                                 correlations=correlations,
                                 outfile=outfile,
                                 submit=True,
                                 job_options="-l mem_free=36G")


@follows(countMREsOverLncs)
@transform(countMREsOverLncs,
           regex("mre.dir/(.+)-(.+)-targetScan-mre_counts.tsv"),
           r"images.dir/\1-\2-mre_counts.histogram.png")
def plotMreCounts(infile, outfile):
    '''
    Plot histograms of MRE counts over genes and lncRNAs
    '''

    PipelineProject036.plotMreCounts(infile,
                                     outfile)


@follows(countMREsOverLncs,
         plotMreCounts,
         mkdir("stats.dir"),
         sharedMREs)
@collate(countMREsOverLncs,
         regex("mre.dir/(.+)-(.+)-targetScan-mre_counts.tsv"),
         r"images.dir/\1-MRE_density.png")
def plotMREDensity(infiles, outfile):
    '''
    Plot density of MRE per nt across lncRNA sets and test
    for differences between gene sets
    '''

    cor_file = [x for x in infiles if re.search("correlation", x)][0]
    prox_file = [p for p in infiles if re.search("proximal", p)][0]
    rand_file = [r for r in infiles if re.search("random", r)][0]
    anti_file = [a for a in infiles if re.search("anticorrelation", a)][0]

    refs = PARAMS['refs'].split(",")
    ref_gtf = [rg for rg in refs if re.search("refcoding", rg)][0]
    lnc_gtf = [lg for lg in refs if re.search("lncRNA", lg)][0]

    PipelineProject036.plotMreDensity(cor_file=cor_file,
                                      prox_file=prox_file,
                                      rand_file=rand_file,
                                      anti_file=anti_file,
                                      ref_gtf=ref_gtf,
                                      lnc_gtf=lnc_gtf,
                                      outfile=outfile,
                                      submit=True,
                                      job_options="-l mem_free=24G")


@follows(plotMREDensity)
@collate("stats.dir/*-MRE_density-stats.tsv",
         regex("stats.dir/(.+)-MRE_density-stats.tsv"),
         r"stats.dir/MRE_density-stats_summary.tsv")
def combineDensityStats(infiles, outfile):
    list_of_files = [l for l in infiles]
    PipelineProject036.mergeStatTable(list_of_files, outfile)


@follows(combineDensityStats)
@transform(combineDensityStats,
           suffix(".tsv"),
           ".load")
def loadDensityStats(infile, outfile):
    P.load(infile, outfile)


@follows(countMREsOverLncs)
@transform(countMREsOverLncs,
           regex("mre.dir/(.+)-(.+)-targetScan-mre_counts.tsv"),
           r"images.dir/\1-\2-mre_counts.violin.png")
def plotViolinMreCounts(infile, outfile):
    '''
    Plot violin plots of MRE counts over genes and lncRNAs
    '''

    PipelineProject036.plotViolinCounts(infile,
                                        outfile)


@follows(sharedMREs,
         loadDensityStats)
@transform(sortMREs,
           regex("mre.dir/(.+)-(.+)-targetScan-mre.sorted.gtf.gz"),
           add_inputs([r"lncRNA_classification.dir/\1-\2_pairs.gtf",
                       r"mre.dir/\1-\2_pairs-mre_shared.tsv"]),
           r"mre.dir/\1-\2_pairs-shared_counts.tsv")
def countSharedMREs(infiles, outfile):
    '''
    count the number of MREs for each shared miRNA element between each
    gene and lncRNA pair
    '''

    mre_file = infiles[0]
    pairs_gtf = infiles[1][0]
    shared_file = infiles[1][1]

    PipelineProject036.countSharedMREs(mre_file,
                                       pairs_gtf,
                                       shared_file,
                                       outfile,
                                       submit=True,
                                       job_options="-l mem_free=16G")


@follows(sharedMREs,
         countSharedMREs,
         plotMREDensity)
@collate(sharedMREs,
         regex("mre.dir/(.+)-(.+)lncRNA_merged-(.+)_pairs-mre_shared.tsv"),
         r"images.dir/\1-\2shared_counts.png")
def plotSharedMreCounts(infiles, outfile):
    '''
    plot distributions of counts of MREs for
    shared miRNA elements and test for differences between gene sets
    '''

    cor_file = [x for x in infiles if re.search("correlation", x)][0]
    prox_file = [p for p in infiles if re.search("proximal", p)][0]
    random_file = [r for r in infiles if re.search("random", r)][0]
    anti_file = [a for a in infiles if re.search("anticorrelation", a)][0]

    PipelineProject036.plotSharedCounts(cor_file=cor_file,
                                        prox_file=prox_file,
                                        random_file=random_file,
                                        anti_file=anti_file,
                                        outfile=outfile,
                                        submit=True,
                                        job_options="-l mem_free=16G")


@follows(plotSharedMreCounts,
         sharedMREs)
@collate("stats.dir/*-MRE_shared-stats.tsv",
         regex("stats.dir/(.+)-MRE_shared-stats.tsv"),
         r"stats.dir/MRE_shared-stats_summary.tsv")
def combineSharedStats(infiles, outfile):
    list_of_files = [l for l in infiles]
    PipelineProject036.mergeStatTable(list_of_files, outfile)


@follows(combineSharedStats)
@transform(combineSharedStats,
           suffix(".tsv"),
           ".load")
def loadSharedStats(infile, outfile):
    P.load(infile, outfile)


###################################################################
###################################################################
###################################################################


@follows(countSharedMREs, mkdir("maf.dir"))
@transform("/ifs/mirror/ucsc/mm10/multiz60way/chr*.maf.gz",
           regex("(.+)/(.+).maf.gz"),
           r"maf.dir/\2.maf")
def gunzipMaf(infile, outfile):
    '''
    Temporarily gunzip maf files for downstream functions
    '''

    statement = ''' zcat %(infile)s > %(outfile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(gunzipMaf)
@transform(gunzipMaf,
           regex("maf.dir/(.+).maf"),
           r"maf.dir/\1.mafindex")
def BuildMafIndices(infile, outfile):
    '''
    Build bio-python maf indices
    '''

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/mreConserve.py
    --log=%(outfile)s.log
    --maf-file=%(infile)s
    --build-index
    --maf-directory=maf.dir
    --species-build=%(genome)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################


@follows(annotateMreGTF,
         cleanMREs, sortMREs,
         BuildMafIndices, gunzipMaf)
@transform(sortMREs,
           regex("mre.dir/(.+)-(.+)-(.+)-targetScan-mre.sorted.gtf.gz"),
           add_inputs(r"lncRNA_classification.dir/\1-\2-\3_pairs.gtf"),
           r"mre.dir/\1-\2-\3-conservation.tsv.gz")
def getConsScores(infile, outfile):
    '''
    Retrieve phastCons scores, calculate p-value based
    on empirical cumulative distribution and return
    table of MREs, conservation score and p-value.
    P-value based on local sequence conservation, preserving
    dinucleotides and controlling for synteny and exon splice
    junctions.

    Need to limit number of jobs used.  Hardcode at 1000, for now.
    '''

    mre_file = infile[0]
    target_file = infile[1]
    job_options = "-l mem_free=2G"
    # set cluster priority queue
    if PARAMS['mre_priority'] == 1:
        c_queue = "all.q"
    elif PARAMS['mre_priority'] == 2:
        c_queue = "pairsdb.q"
    elif PARAMS['mre_priority'] == 3:
        c_queue = "mpi.q"
    else:
        c_queue = "mpi.q"
        E.info("no queue priority set, defaulting to %s" % c_queue)

    statement = '''zcat %(mre_file)s |
    python %(pipeline_scriptsdir)s/farm.py
    --split-at-lines=1500
    --output-header
    --cluster-queue=%(c_queue)s
    --log=%(outfile)s.farm.log
    "python /ifs/projects/proj036/pipeline_project036_timeseries/src/mreConserve.py 
    --log=%(outfile)s.log
    --target-file=%(target_file)s
    --maf-directory=maf.dir
    --cons-file=%(mafopt_cons_file)s
    --species-build=%(genome)s
    --permutations=%(mre_perms)s" | gzip > %(outfile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(getConsScores)
@transform("maf.dir/*.maf",
           suffix(".maf"),
           ".maf.gz")
def gzipMAF(infile, outfile):

    statement = '''gzip %(infile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(getConsScores)
@transform(getConsScores,
           suffix(".tsv.gz"),
           ".load")
def loadConsScores(infile, outfile):

    tmpfile = P.getTempFilename(shared=True)

    statement = '''zcat %(infile)s
    | awk '{if(NR==1) {printf("mre_id\\t%%s\\n", $0)} else {printf("%%s\\n", $0)}}'
    > %(tmpfile)s'''
    
    P.run()

    P.load(tmpfile, outfile)

    statement = '''rm -f %(tmpfile)s'''

    P.run()


@follows(sharedMREs,
         getConsScores)
@collate(getConsScores,
         regex("mre.dir/(.+)-lncRNA_merged-(.+)-conservation.tsv.gz"),
         r"images.dir/\1-sig_conserved_MREs.png")
def plotSigMREDist(infiles, outfile):
    '''
    plot distributions of counts of MREs for
    shared miRNA elements and test for differences between gene sets
    '''
    cor_file = [x for x in infiles if re.search("correlation", x)][0]
    prox_file = [p for p in infiles if re.search("proximal", p)][0]
    random_file = [r for r in infiles if re.search("random", r)][0]

    PipelineProject036.plotSigMREs(cor=cor_file,
                                   prox=prox_file,
                                   random=random_file,
                                   outfile=outfile,
                                   submit=True)


@follows(plotSigMREDist,
         sharedMREs,
         mkdir("stats.dir"))
@collate("stats.dir/*-sig_conserved-stats.tsv",
         regex("stats.dir/(.+)-sig_conserved-stats.tsv"),
         r"stats.dir/Sig_conserved-stats_summary.tsv")
def combineSigConsStats(infiles, outfile):
    list_of_files = [l for l in infiles]
    PipelineProject036.mergeStatTable(list_of_files, outfile)


@follows(combineSigConsStats)
@transform(combineSigConsStats,
           suffix(".tsv"),
           ".load")
def loadSigConsStats(infile, outfile):
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


# @follows(loadSharedStats,
#          loadDensityStats,
#          loadConsScores,
#          loadSigConsStats,
#          mkdir("ceRNA.dir"))
@follows(loadSharedStats,
         loadDensityStats,
         mkdir("ceRNA.dir"))
@collate("mre.dir/*.tsv",
         regex("mre.dir/(.+)-lncRNA_merged-correlation(.+).tsv"),
         r"ceRNA.dir/\1-candidates.tsv")
def ceRNACandidates(infiles, outfile):
    '''
    From clustering analysis
    Take all lncRNA:gene pairs from MRE density and miRNA shared.
    Take intersection of each as candidates for functional analysis
    '''

    counts_file = [i for i in infiles if re.search("counts", i)][0]
    shared_file = [s for s in infiles if re.search("shared", s)][0]
    header = counts_file.split("/")[-1].split("-")[0]
    expr = "-".join([header, "average_expression.tsv.gz"])
    expression = "/".join([os.getcwd(), "expression.dir", expr])
    annotations = "-".join([header, "lncRNA_merged-correlation_pairs.tsv"])
    annotations = "/".join([os.getcwd(),
                            "correlations.dir",
                            annotations])
    refs = PARAMS['refs'].split(",")
    lncs_gtf = [l for l in refs if re.search("lncRNA", l)][0]
    gene_gtf = [g for g in refs if re.search("refcoding", g)][0]

    threshold = float(PARAMS['lncrna_correlation_threshold'])
    PipelineProject036.candidateCeRNAs(counts_file,
                                       shared_file,
                                       lncs_gtf,
                                       gene_gtf,
                                       annotations,
                                       threshold,
                                       expression,
                                       outfile,
                                       submit=True,
                                       job_options="-l mem_free=8G")


@follows(ceRNACandidates)
@collate("mre.dir/*.tsv",
         regex("mre.dir/(.+)-(.+)_lncRNA_merged-[correlation|proximal](.+).tsv"),
         r"ceRNA.dir/\1_\2-candidates.tsv")
def SubsetceRNACandidates(infiles, outfile):
    '''
    From differential expression analysis
    Take lncRNA:gene pairs from MRE density and miRNA shared.
    Take intersection of each as candidates for functional analysis
    '''

    counts_file = [i for i in infiles if re.search("counts", i)][0]
    shared_file = [s for s in infiles if re.search("shared", s)][0]
    header = counts_file.split("/")[-1].split("-")[0]
    expr = "-".join([header, "average_expression.tsv.gz"])
    expression = "/".join([os.getcwd(), "expression.dir", expr])
    annotations = "-".join([header, "lncRNA_merged-correlation_pairs.tsv"])
    annotations = "/".join([os.getcwd(),
                            "correlations.dir",
                            annotations])
    refs = PARAMS['refs'].split(",")
    lncs_gtf = [l for l in refs if re.search("lncRNA", l)][0]
    gene_gtf = [g for g in refs if re.search("refcoding", g)][0]

    threshold = float(PARAMS['lncrna_correlation_threshold'])
    PipelineProject036.candidateCeRNAs(counts_file,
                                       shared_file,
                                       lncs_gtf,
                                       gene_gtf,
                                       annotations,
                                       threshold,
                                       expression,
                                       outfile,
                                       submit=True,
                                       job_options="-l mem_free=16G")


@follows(ceRNACandidates,
         SubsetceRNACandidates)
@transform("ceRNA.dir/*-candidates.tsv",
           suffix(".tsv"),
           ".load")
def loadCeRNACandidates(infile, outfile):
    if IOTools.isEmpty(infile):
        P.touch(outfile)
    else:
        P.load(infile, outfile)


@follows(loadCeRNACandidates)
@transform(ceRNACandidates,
           regex("ceRNA.dir/(.+)-candidates.tsv"),
           add_inputs([r"lncRNA_merged.gtf.gz",
                       r"mre.dir/\1-lncRNA_merged-correlation_pairs-mre_shared.tsv"]),
           r"ceRNA.dir/\1-ceRNAs_annotated.tsv")
def annotateCeRNAs(infiles, outfile):
    '''
    Annotate candidate ceRNA and ceRNAt partners
    '''

    ceRNA_file = infiles[0]
    lnc_gtf = infiles[1][0]
    mre_file = infiles[1][1]

    PipelineProject036.annotateCeRNAs(ceRNA_file=ceRNA_file,
                                      lnc_gtf=lnc_gtf,
                                      mre_file=mre_file,
                                      outfile=outfile,
                                      submit=True,
                                      job_options="-l mem_free=4G")


@follows(loadCeRNACandidates)
@transform(SubsetceRNACandidates,
           regex("ceRNA.dir/(.+)_(.+)-candidates.tsv"),
           add_inputs([r"lncRNA_merged.gtf.gz",
                       r"mre.dir/\1-\2_lncRNA_merged-correlation_pairs-mre_shared.tsv"]),
           r"ceRNA.dir/\1_\2-ceRNAs_annotated.tsv")
def annotateSubsetCeRNAs(infiles, outfile):
    '''
    Annotate candidate ceRNA and ceRNAt partners
    '''

    ceRNA_file = infiles[0]
    lnc_gtf = infiles[1][0]
    mre_file = infiles[1][1]

    PipelineProject036.annotateCeRNAs(ceRNA_file=ceRNA_file,
                                      lnc_gtf=lnc_gtf,
                                      mre_file=mre_file,
                                      outfile=outfile,
                                      submit=True,
                                      job_options="-l mem_free")


@follows(annotateCeRNAs)
@transform(annotateCeRNAs,
           regex("ceRNA.dir/(.+)-ceRNAs_annotated.tsv"),
           add_inputs(r"expression.dir/\1-average_expression.tsv.gz"),
           r"ceRNA.dir/\1-ceRNAs_upregulated.tsv")
def selectUpRegulatedCeRNAs(infiles, outfile):
    '''
    Select up-regulated ceRNAs as KO targets
    '''

    cerna_file = infiles[0]
    expression_file = infiles[1]
    PipelineProject036.filterAnnotatedCeRNAs(cerna_file=cerna_file,
                                             expression_file=expression_file,
                                             direction="up",
                                             outfile=outfile,
                                             submit=True,
                                             job_options="-l mem_free=4G")


@follows(annotateSubsetCeRNAs)
@transform(annotateSubsetCeRNAs,
           regex("ceRNA.dir/(.+)_(.+)-ceRNAs_annotated.tsv"),
           add_inputs(r"expression.dir/\1-average_expression.tsv.gz"),
           r"ceRNA.dir/\1_\2-ceRNAs_upregulated.tsv")
def selectUpRegulatedSubsetCeRNAs(infiles, outfile):
    '''
    Select up-regulated ceRNAs as KO targets
    '''

    cerna_file = infiles[0]
    expression_file = infiles[1]
    PipelineProject036.filterAnnotatedCeRNAs(cerna_file=cerna_file,
                                             expression_file=expression_file,
                                             direction="up",
                                             outfile=outfile,
                                             submit=True,
                                             job_options="-l mem_free=4G")


@follows(selectUpRegulatedSubsetCeRNAs,
         selectUpRegulatedCeRNAs)
@transform(selectUpRegulatedCeRNAs,
           regex("ceRNA.dir/(.+)-ceRNAs_upregulated.tsv"),
           add_inputs(r"expression.dir/\1-average_expression.tsv.gz"),
           r"ceRNA.dir/\1-ceRNAs_upregulated.genexpr.tsv")
def addGeneExpression(infiles, outfile):
    '''
    Annotated ceRNA candidates with expression profiles of
    partner genes
    '''

    cerna_file = infiles[0]
    expression_file = infiles[1]

    PipelineProject036.annotateGeneExpr(cerna_file,
                                        expression_file,
                                        outfile,
                                        submit=True,
                                        job_options="-l mem_free=4G")


@follows(selectUpRegulatedSubsetCeRNAs,
         selectUpRegulatedCeRNAs)
@transform(selectUpRegulatedSubsetCeRNAs,
           regex("ceRNA.dir/(.+)_(.+)-ceRNAs_upregulated.tsv"),
           add_inputs(r"expression.dir/\1-average_expression.tsv.gz"),
           r"ceRNA.dir/\1_\2-ceRNAs_upregulated.genexpr.tsv")
def addSubsetGeneExpression(infiles, outfile):
    '''
    Annotated ceRNA candidates with expression profiles of
    partner genes
    '''

    cerna_file = infiles[0]
    expression_file = infiles[1]

    PipelineProject036.annotateGeneExpr(cerna_file,
                                        expression_file,
                                        outfile,
                                        submit=True,
                                        job_options="-l mem_free=5G")

@follows(addGeneExpression)
@transform(addGeneExpression,
           regex("ceRNA.dir/(.+)-ceRNAs_upregulated.genexpr.tsv"),
           add_inputs([r"mre.dir/\1-lncRNA_merged-correlation-targetScan-mre.sorted.gtf.gz",
                       r"merged_gtf.dir/merged.gtf.gz"]),
           r"stats.dir/\1-miRNA_ceRNA-enrichment.tsv")
def testMirEnrichment(infiles, outfile):
    '''
    Test each miRNA for enrichment in shared ceRNA vs. partner genes
    '''

    cerna_file = infiles[0]
    mre_file = infiles[1][0]
    pairs_gtf = infiles[1][1]
    mirna_file = PARAMS['mre_filter_file']

    job_options = "-l mem_free=8G"

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/mirs2enrichment.py
    --log=%(outfile)s.log
    --mre-gtf-file=%(mre_file)s
    --transcript-gtf-file=%(pairs_gtf)s
    --mirna-tsv-file=%(mirna_file)s
    %(cerna_file)s
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
## Pipeline tasks associated with analysing and comparing to     ##
## publicly available expression data sets - ENCODE and FANTOM5  ##
###################################################################
###################################################################
ENCODE_DIR = PARAMS['annotations_encode_dir']
FANTOM5_DIR = PARAMS['annotations_fantom5_dir']
FASTQ_SUFFIXES = ("*.fastq.1.gz",
                  "*.fastq.2.gz",
                  "*.fastq.gz")
FASTQ_DIR = PARAMS['fastq_dir']
FASTQ_FILES = tuple([os.path.join(FASTQ_DIR, suffix_name)
                     for suffix_name in FASTQ_SUFFIXES])
FASTQ_REGEX = regex(r"%s/(\S+).fastq.1.gz" % FASTQ_DIR)
SE_REGEX = regex(r"%s/(\S+).fastq.gz" % FASTQ_DIR)
FASTQ_PAIR = r"%s/\1.fastq.2.gz" % FASTQ_DIR
PAIR_EXIST = any([True for pe in os.listdir("%s" % FASTQ_DIR) if re.search("fastq.2.gz",
                                                                           pe)])
#PAIR_EXIST = None

# get tpm for lncRNAs in current data
@follows(mkdir("transcripts.dir"))
@transform("*.gtf.gz",
           regex("(.+).gtf.gz"),
           r"transcripts.dir/\1.fa")
def makeRepTranscripts(infile, outfile):
    '''
    make a single representative transcript for each
    lncRNA - put into a multi-fasta file
    '''

    genome_file = "/".join([PARAMS['genome_dir'], PARAMS['genome']])

    statement = '''
    zcat %(infile)s | python %(scriptsdir)s/gtf2gtf.py --method=merge-transcripts
    --use-gene-id
    --log=%(outfile)s.merge.log |
    python %(scriptsdir)s/gff2fasta.py
    --genome-file=%(genome_file)s
    --is-gtf
    --log=%(outfile)s.log
    > %(outfile)s
    '''

    P.run()


@follows(makeRepTranscripts)
@transform(makeRepTranscripts,
           regex("transcripts.dir/(.+).fa"),
           r"transcripts.dir/\1.spliced.fa")
def makeSplicedCatalog(infile, outfile):
    '''
    make fasta file of spliced transcript sequences
    '''

    statement = '''
    python %(scriptsdir)s/cgat_fasta2cDNA.py
    --log=%(outfile)s.log
    %(infile)s
    > %(outfile)s
    '''

    P.run()


if PARAMS['transcript_program'] == "kallisto":
    @follows(mkdir("kallisto_index.dir"),
             makeSplicedCatalog)
    @transform(makeSplicedCatalog,
               regex("transcripts.dir/(.+).spliced.fa"),
               r"kallisto_index.dir/\1.index")
    def makeKallistoIndex(infile, outfile):
        '''
        Make a kallisto index file from a multi-fasta of
        spliced transcript sequences
        '''

        statement = '''
        python %(scriptsdir)s/fastq2tpm.py
        --method=make_index
        --prorgam=kallisto
        --index-fasta=%(infile)s
        --output-file=%(outfile)s
        --log=%(outfile)s.log
        '''

        P.run()


        if os.path.exists(FASTQ_PAIR):
            @follows(mkdir("tpm.dir"))
            @transform(FASTQ_FILES,
                       FASTQ_REGEX,
                       add_inputs([makeKallistoIndex,
                                   FASTQ_PAIR]),
                       r"tpm.dir/\1.abundance")
            def quantifyWithKallisto(infiles, outfile):
                '''
                Quantify gene/transcript expression with Kallisto
                '''

                fastq1 = infiles[0]
                fastq2 = infiles[1][1]

                fastqs = ",".join([fastq1, fastq2])
                index_file = infiles[1][0]
                out_dir = ".".join(outfile.split(".")[:-1])
                job_threads = 4
                job_memory = "1G"
    
                count_file = "/".join([out_dir, "abundance.tsv"])
                # make the output directory
                os.system("mkdir %s" % out_dir)

                statement = '''
                python %(scriptsdir)s/fastq2tpm.py
                --log=%(outfile)s.log
                --program=kallisto
                --method=quant
                --index-file=%(index_file)s
                --output-directory=%(out_dir)s
                --use-bias
                --bootstraps=26
                --seed=%(seed)i
                --threads=%(job_threads)s
                --just-text
                %(fastqs)s;
                cp %(count_file)s %(outfile)s
                '''

                P.run()


            @follows(quantifyWithKallisto)
            @collate(quantifyWithKallisto,
                     regex("tpm.dir/(.+)-(.+)-(.+).abundance"),
                     r"tpm.dir/\1.tpm")
            def mergeKallistoRuns(infiles, outfile):
                '''
                Merge all raw tpm estimates from kallisto across each
                condition
                '''

                infiles = " ".join(infiles)
                job_memory = "2G"

                statement = '''
                python %(scriptsdir)s/combine_tables.py
                --columns=1
                --take=5
                --use-file-prefix
                --regex-filename='(.+).abundance'
                --log=%(outfile)s.log
                %(infiles)s
                > %(outfile)s'''

                P.run()


elif PARAMS['transcript_program'] == "sailfish":
    @follows(mkdir("sailfish_index.dir"),
             makeSplicedCatalog)
    @transform(makeSplicedCatalog,
               regex("transcripts.dir/(.+).spliced.fa"),
               r"sailfish_index.dir/sa.bin")
    def makeSailfishIndex(infile, outfile):
        '''
        Make a sailfish index file from a multi-fasta of
        spliced transcript sequences
        '''

        outdir = "/".join(outfile.split("/")[:-1])
        job_threads = 8
        job_memory = "12G"
        statement = '''
        python %(scriptsdir)s/fastq2tpm.py
        --method=make_index
        --program=sailfish
        --index-fasta=%(infile)s
        --kmer-size=31
        --threads=%(job_threads)s
        --output-directory=%(outdir)s
        --log=%(outfile)s.log
        '''

        P.run()


    if PAIR_EXIST:
        @follows(mkdir("tpm.dir"),
                 makeSailfishIndex)
        @transform(FASTQ_FILES,
                   FASTQ_REGEX,
                   add_inputs([FASTQ_PAIR,
                               makeSailfishIndex]),
                   r"tpm.dir/\1.tpm")
        def quantifyWithSailfish(infiles, outfile):
            '''
            Quantify gene/transcript expression with sailfish
            '''

            fastq1 = infiles[0]
            # need to check that fastq2 file exists
            # if not, run as single-end
            fastq2 = infiles[1][0]

            index_dir = "/".join(infiles[1][1].split("/")[:-1])
            out_dir = ".".join(outfile.split(".")[:-1])
            job_threads = 8
            job_memory = "6G"

            count_file = "/".join([out_dir, "quant.sf"])

            fastqs = ",".join([fastq1, fastq2])
            state1 = '''
            python %(scriptsdir)s/fastq2tpm.py
            --log=%(outfile)s.log
            --program=sailfish
            --method=quant
            --paired-end
            --index-file=%(index_dir)s
            --output-directory=%(out_dir)s
            --library-type=%(transcript_library)s
            --threads=%(job_threads)s
            %(fastqs)s;
            '''

            state2 = '''cat  %(count_file)s |
            awk 'BEGIN {printf("Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n")} 
            {if(NR > 11) {print $0}}' >  %(outfile)s'''
            statements = [state1, state2]

            P.run()


    else:
        @follows(mkdir("tpm.dir"),
                 makeSailfishIndex)
        @transform(FASTQ_FILES,
                   SE_REGEX,
                   add_inputs(makeSailfishIndex),
                   r"tpm.dir/\1.tpm")
        def quantifyWithSailfish(infiles, outfile):
            '''
            Quantify gene/transcript expression with sailfish
            '''

            fastqs = infiles[0]
            # need to check that fastq2 file exists
            # if not, run as single-end
            index_dir = "/".join(infiles[1].split("/")[:-1])
            out_dir = ".".join(outfile.split(".")[:-1])
            job_threads = 8
            job_memory = "6G"

            count_file = "/".join([out_dir, "quant.sf"])

            state1 = '''
            python %(scriptsdir)s/fastq2tpm.py
            --log=%(outfile)s.log
            --program=sailfish
            --method=quant
            --index-file=%(index_dir)s
            --output-directory=%(out_dir)s
            --library-type=%(transcript_library)s
            --threads=%(job_threads)s
            %(fastqs)s;
            '''

            state2 = '''cat  %(count_file)s |
            awk 'BEGIN {printf("Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n")} 
            {if(NR > 11) {print $0}}' >  %(outfile)s'''
            statements = [state1, state2]

            P.run()



    @follows(quantifyWithSailfish)
    @collate(quantifyWithSailfish,
             regex("tpm.dir/(.+)-(.+)-(.+).tpm"),
             r"tpm.dir/\1.tpm")
    def mergeSailfishRuns(infiles, outfile):
        '''
        Merge all raw tpm estimates from kallisto across each
        condition
        '''

        infiles = " ".join(infiles)
        job_memory = "2G"

        statement = '''
        python /ifs/devel/michaelm/cgat/CGAT/scripts/combine_tables.py
        --columns=1
        --take=4
        --use-file-prefix
        --regex-filename='(.+).tpm'
        --log=%(outfile)s.log
        %(infiles)s
        > %(outfile)s'''

        P.run()


# filter fantom5 expression data for regions that are lifted over
@follows(mkdir("fantom5.dir"))
@files("%s/%s.cage_peak_phase1and2combined_tpm.osc.txt.gz" % (FANTOM5_DIR,
                                                              PARAMS['annotations_fantom5_genome']),
       "fantom5.dir/%s.cage_peak_tpm.bed.gz" % PARAMS['annotations_fantom5_genome'])
def getFantomPeaksBed(infile, outfile):
    '''
    Pull peak co-ordinates out of expression file into BED6 format
    Needs to be BED6 for later compatibility
    '''

    job_memory = "2G"

    # need to grep for chromosomes too
    statement = '''
    zcat %(infile)s | grep -v '##' | grep 'chr' | cut -f 1 | awk '{if(NR > 2) {print $0}}' |
    sed 's/:/\\t/g' | sed 's/\\.\\./\\t/g' | sed 's/,/\\t/g' | 
    awk '{printf("%%s\\t%%i\\t%%i\\t%%s:%%i..%%i\\t%%s\\t0\\t%%s\\n", $1,$2,$3,$1,$2,$3,$4,$4)}' | 
    gzip > %(outfile)s
    '''
    P.run()


@follows(getFantomPeaksBed)
@files("fantom5.dir/%s.cage_peak_tpm.bed.gz" % PARAMS['annotations_fantom5_genome'],
       "fantom5.dir/%s.cage_peak_tpm.bed.gz" % PARAMS['genome'])
def liftOverBed(infile, outfile):
    '''
    Use UCSC liftOver tool to convert fantom5 genome co-ordinates to desired genome
    '''

    job_memory = "2G"
    tmp_file = P.getTempFilename(shared=True)

    statement = '''
    liftOver %(infile)s %(annotations_chain_file)s %(tmp_file)s fantom5.dir/unlifted.bed;
    cat %(tmp_file)s | gzip > %(outfile)s;'''

    P.run()


@follows(liftOverBed)
@files(["fantom5.dir/%s.cage_peak_tpm.bed.gz" % PARAMS['annotations_fantom5_genome'],
        "fantom5.dir/unlifted.bed"],
       "fantom5.dir/%s.cage_peak_tpm.filtered.bed.gz" % PARAMS['annotations_fantom5_genome'])
def filterLiftedBed(infiles, outfile):
    '''
    Filter out unlifted intervals from original bed file
    '''

    job_memory = "2G"
    bed_file = infiles[0]
    unlifted = infiles[1]

    statement = '''
    bedtools subtract -f 0.75 -s -N -a %(bed_file)s -b %(unlifted)s
    | gzip > %(outfile)s'''

    P.run()


@follows(filterLiftedBed)
@files(["fantom5.dir/%s.cage_peak_tpm.filtered.bed.gz" % PARAMS['annotations_fantom5_genome'],
        "fantom5.dir/%s.cage_peak_tpm.bed.gz" % PARAMS['genome']],
       "fantom5.dir/merged.bed.gz")
def mergeBedGenomeCoordinates(infiles, outfile):
    '''
    Generate a merged file of old and new genome co-ordinates
    '''

    job_memory = "2G"
    old_genome = infiles[0]
    new_genome = infiles[1]

    old_temp = P.getTempFilename(shared=True)
    new_temp = P.getTempFilename(shared=True)

    statement = '''
    zcat %(old_genome)s > %(old_temp)s;
    zcat %(new_genome)s > %(new_temp)s;
    paste %(old_temp)s %(new_temp)s | gzip > %(outfile)s;
    rm -rf %(old_temp)s %(new_temp)s
    '''

    P.run()

@follows(mergeBedGenomeCoordinates)
@files(["%s/%s.cage_peak_phase1and2combined_tpm.osc.txt.gz" % (FANTOM5_DIR,
                                                               PARAMS['annotations_fantom5_genome']),
       "fantom5.dir/merged.bed.gz"],
       "fantom5.dir/%s.cage_peak_tpm.filtered.txt.gz" % PARAMS['genome'])
def assignNewCoordinates(infiles, outfile):
    '''
    Match up and replace old genome expression table genome coordinates
    with new genome co-ordinates
    '''

    infile = infiles[0]
    map_file = infiles[1]
    job_memory = "7G"

    statement = '''
    python /ifs/devel/projects/proj036/bed2tpm.py
    --log=%(outfile)s.log
    --task=merge_coordinates
    --mapped-bed-file=%(map_file)s
    %(infile)s
    | gzip > %(outfile)s
    '''

    P.run()


@follows(assignNewCoordinates)
@files("lncRNA_merged.gtf.gz",
       "lncRNA_merged.bed.gz")
def gtf2Bed(infile, outfile):
    '''
    Merge into a single gene interval, then convert gtf file to bed
    format. Extend 5' region by 100bp
    to account for poor TSS annotation in lncRNA models
    '''

    genome_file = "/".join([PARAMS['genome_dir'], PARAMS['genome']])

    statement = '''
    zcat %(infile)s | python %(scriptsdir)s/gtf2gtf.py
    --method=merge-transcripts -v 0 |
    python %(scriptsdir)s/gff2bed.py
    --is-gtf --set-name=gene_id -v 0 | 
    sort -k1,1 -k2,2n |
    python %(scriptsdir)s/bed2bed.py
    --method=merge --merge-by-name -v 0 |
    python %(scriptsdir)s/bed2bed.py
    --log=%(outfile)s.extend.log --method=extend
    --offset=100 --genome-file=%(genome_file)s.fasta |
    gzip  > %(outfile)s'''

    P.run()


@follows(assignNewCoordinates,
         gtf2Bed)
@transform("fantom5.dir/%s.cage_peak_tpm.bed.gz" % PARAMS['genome'],
           regex("fantom5.dir/(.+).cage_peak_tpm.bed.gz"),
           add_inputs("lncRNA_merged.bed.gz"),
           r"fantom5.dir/\1.lncRNA_intervals.cage.bed.gz")
def intersectCageAndLncRNA(infiles, outfile):
    '''
    Intersect lncRNA and cage tag peak intervals - output are lncRNA intervals
    that intersect with cage peaks and overlap.
    '''

    lnc_file = infiles[1]
    bed_file = infiles[0]

    job_memory = "6G"

    statement = '''
    bedtools intersect -wo -a %(lnc_file)s -b %(bed_file)s
    | gzip > %(outfile)s'''

    P.run()


@follows(intersectCageAndLncRNA)
@transform(intersectCageAndLncRNA,
           regex("fantom5.dir/(.+).lncRNA_intervals.cage.bed.gz"),
           add_inputs(r"fantom5.dir/\1.cage_peak_tpm.filtered.txt.gz"),
           r"fantom5.dir/\1-lncRNA_merged-cage_tpm.tsv.gz")
def mergeTpmByLncRNA(infiles, outfile):
    '''
    Select the tpm data for intervals that intersect with lncRNAs
    '''

    job_memory = "4G"

    intervals_file = infiles[0]
    expr_file = infiles[1]

    statement = '''
    python /ifs/devel/projects/proj036/bed2tpm.py
    --task=intersect_transcripts
    --transcript-bed-file=%(intervals_file)s
    --log=%(outfile)s.log
    %(expr_file)s
    | gzip > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################


@follows(mkdir("Ig_csr.dir"))
@transform("*.bam",
           regex("(.+)-(.+)-(.+).bam"),
           r"Ig_csr.dir/\1-\2-\3-csr_counts.tsv")
def countClassSwitched(infile, outfile):
    '''
    Using paired reads, count the number of reads falling
    into properly class switch recombined Ig heavy chain
    transcripts
    '''

    statement = '''
    python /ifs/projects/proj036/pipeline_project036_timeseries/src/csr_counting.py
    --log=%(outfile)s.log
    --j-gtf=%(csr_j_gtf)s
    --c-gtf=%(csr_c_gtf)s
    --ig-coordinates=%(csr_locus)s
    %(infile)s
    > %(outfile)s'''

    P.run()


@follows(countClassSwitched)
@collate("Ig_csr.dir/*-csr_counts.tsv",
         regex("Ig_csr.dir/(.+)-(.+)-(.+)-csr_counts.tsv"),
         r"Ig_csr.dir/\1-csr_counts.tsv")
def aggregateClassSwitched(infiles, outfile):
    '''
    aggregate counts over class switched Ig transcripts
    for each condition
    '''

    infiles = " ".join(infiles)

    statement = '''
    python %(scriptsdir)s/combine_tables.py
    --columns=1
    --use-file-prefix
    --regex-filename='(.+)-csr_counts.tsv'
    --log=%(outfile)s.log
    %(infiles)s
    > %(outfile)s'''

    P.run()


@follows(clusterGeneOntology,
         geneOntologyCoreGenes,
         geneOntologyConditionGenes,
         clusterGeneOntology)
def go_enrichment():
    pass


@follows(aggregateClassSwitched)
def csr():
    pass
###################################################################
###################################################################
###################################################################


@follows(loadConsScores,
         countMREsOverLncs,
         ceRNACandidates)
def find_mre():
    pass


@follows(selectUpRegulatedCeRNAs,
         selectUpRegulatedSubsetCeRNAs)
def candidate_ceRNAs():
    pass
###################################################################
###################################################################
###################################################################

# get FANTOM5 expression levels that intersect lncRNAs
@follows(mergeTpmByLncRNA)
def fantom_expression():
    pass
###################################################################
###################################################################
###################################################################


@follows(csr,
         go_enrichment,
         find_mre)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish():
    '''publish report and data.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
