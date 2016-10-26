'''
lnc2correlate.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Classification of lncRNAs based on their correlation of expression
with protein-coding genes and cluster eigengenes

TODO:
* bring in line with CGAT scripts style guide
* remove redundancy in functions and tasks

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import pandas as pd
from rpy2.robjects import pandas2ri
import CGAT.Experiment as E
import PipelineProject036 as P36


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--infile", dest="infile", type="string",
                      help="input file path")

    parser.add_option("--task", dest="task", type="string",
                      help="analysis task to be executed")

    parser.add_option("--time", dest="timepoints", type="string",
                      help="a comma-separated list of time points measured")

    parser.add_option("--cor-threshold", dest="cor_threshold", type="string",
                      help="threshold at which to filter lncRNA expression"
                      "correlation with eigengene expression")

    parser.add_option("--output-gtf", dest="output_gtf", type="string",
                      help="output gtf file to write attributes to")

    parser.add_option("--summary-file", dest="summary_file", type="string",
                      help="summary file")

    parser.add_option("--lncRNA-class", dest="lncrna_class", type="choice",
                      choices=("intergenic",
                               "antisense",
                               "sense_downstream",
                               "sense_upstream",
                               "antisense_upstream",
                               "antisense_downstream"),
                      help="classification of lncRNA to use for"
                      " correlation analysis")

    parser.add_option("--cor-direction", dest="corr_direction", type="choice",
                      choices=("positive", "negative"), help="direction of "
                      "correlation to classify lncRNAs on")

    parser.add_option("--method", dest="method", type="choice",
                      choices=("temporal", "cross-correlation"),
                      help="correlation method to use, either cross- or"
                      " temporal")

    # add common options (-h/--help, ...) and parse command line

    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    if options.task == "correlate-lncs":
        files = infile.split(",")
        eigenf = files[0]
        lncFile = files[1]
        threshold = float(options.cor_threshold)
        df = P36.correlate_lncsEigen(lncFile=lncFile,
                                     eigenFile=eigenf,
                                     correlation="cross-correlation",
                                     lag=0)

    elif options.task == "filter-correlation":
        threshold = float(options.cor_threshold)
        df = P36.filter_correlation(infile=infile,
                                    threshold=threshold)

    elif options.task == "classify-lncrna":
        files = infile.split(",")
        lnc_file = files[0]
        lnc_gtf = files[1]
        threshold = float(options.cor_threshold)

        df = P36.classifyLncRNA(lnc_list=lnc_file,
                                lnc_gtf=lnc_gtf,
                                lnc_class=options.lncrna_class,
                                direction=options.corr_direction,
                                threshold=threshold)

    elif options.task == "correlate-direction":

        files = infile.split(",")
        cluster_file = files[0]
        gene_express = files[1]
        lncRNA_express = files[2]
        lnc_gtf = files[3]

        threshold = float(options.cor_threshold)
        lncrna_class = options.lncrna_class
        corr_direction = options.corr_direction

        df = P36.classCorrLncRNA(cluster_file,
                                 gene_express,
                                 lncRNA_express,
                                 lnc_gtf,
                                 threshold,
                                 lncrna_class,
                                 corr_direction)

    elif options.task == "correlate-cluster-genes":
        files = infile.split(",")
        cor_file = files[0]
        expr_file = files[1]
        gene_clusters = files[2]
        lnc_clusters = files[3]

        df = P36.correlateLncRNAsWithClusterGenes(cor_file=cor_file,
                                                  expr_file=expr_file,
                                                  set1_clusters=gene_clusters,
                                                  set2_clusters=lnc_clusters,
                                                  correlation="cross-correlation")

    elif options.task == "correlate-eigens":
        infiles = infile.split(",")
        ref_file = infiles[0]
        lnc_file = infiles[1]

        df = P36.correlateEigengenes(ref_eigens=ref_file,
                                     lnc_eigens=lnc_file,
                                     correlation="cross-correlation",
                                     lag=0)

    elif options.task == "filter-correlations":
        df = P36.filterCorrelations(infile)

    elif options.task == "filter-pairs":
        threshold = float(options.cor_threshold)
        df = P36.correlationPairs(infile,
                                  threshold)

    elif options.task == "max-correlation":
        df = P36.maxCorrelationPairs(infile)

    elif options.task == "correlate-all":
        if infile.split(".")[-1] == "gz":
            comp = "gzip"
        else:
            comp = None

        expr_df = pd.read_table(infile, sep="\t",
                                header=0, index_col=0,
                                compression=comp)

        df = P36.correlateLncRNAs(lnc_frame=expr_df,
                                  gene_frame=expr_df,
                                  correlation="cross-correlation",
                                  lag=0)

    df.to_csv(options.stdout,
              sep="\t",
              header=True,
              index_label="gene_id")

    # Write footer and output benchmark information.
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
