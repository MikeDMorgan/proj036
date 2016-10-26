'''
mirs2enrichment.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Test for enrichment of miRNAs amongst ceRNAs and their partner genes

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import PipelineProject036

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

    parser.add_option("--mre-gtf-file", dest="mre_file", type="string",
                      help="GTF of predicted MREs")

    parser.add_option("--transcript-gtf-file", dest="pairs_gtf", type="string",
                      help="GTF of input lncRNAs and genes")

    parser.add_option("--mirna-tsv-file", dest="mirna_file", type="string",
                      help="list of expressed miRNAs to test for enrichment")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    cerna_file = argv[-1]

    mre_file = options.mre_file
    pairs_gtf = options.pairs_gtf
    mirna_file = options.mirna_file

    res_out = PipelineProject036.mirEnrichment(cerna_file=cerna_file,
                                               mre_file=mre_file,
                                               pairs_gtf=pairs_gtf,
                                               mirna_file=mirna_file)                                     

    options.stdout.write("miRNA\toddsRatio\t95CI_low\t95_CI_hi\tpvalue\t"
                         "ncernas\tnpartners\ttcernas\ttpartners\tqvalue\n")
    for result in res_out:
        options.stdout.write("\t".join(map(str, result[0] + result[1])) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
