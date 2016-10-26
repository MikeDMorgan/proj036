'''
mre2mre.py - script for processing MRE prediction files
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----
Processing, filtering and annotating of or with in silico
predicted miRNA recognition elements.

Example::

   python mre2mre.py

Type::

   python mre2mre.py --help

for command line help.

Command line options
--------------------

'''

import sys
import pandas as pd
from rpy2.robjects import pandas2ri
import CGAT.Experiment as E
import PipelineProject036 as P36
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
pandas2ri.activate()


def filterMREsTSV(input_file, filter_set):
    '''
    Filter MREs in a GFF file based on a list of
    miRNA IDs.  Return a generator object.
    '''

    mre_file = IOTools.openFile(input_file, "rb")
    for x in GTF.transcript_iterator(GTF.iterator(mre_file)):
        for mre in x:
            if mre.asDict()['miRNA'] in filter_set:
                entry = GTF.Entry()
                entry.copy(mre)
                yield entry
            else:
                pass
    mre_file.close()


def filterMREsGTF(input_file, filter_file):
    '''
    Filter MREs in a GFF file with intersecting features/genes
    in a supplied GFF or GTF file.  Overlap must be complete,
    i.e. MRE must be contained within feature X.
    MRE.start >= X.start AND MRE.end <= X.end
    '''

    mre_file = IOTools.openFile(input_file, "rb")
    feature_file = IOTools.openFile(filter_file, "rb")

    # use index filter file, select MREs that completely
    # overlap with exons
    mre_it = GTF.transcript_iterator(GTF.iterator(mre_file))
    gff_idx = GTF.readAndIndex(GTF.iterator(feature_file))

    for mre in mre_it:
        gff_src = gff_idx.get(mre[0].contig,
                              mre[0].start,
                              mre[0].end)
        gff_entry = [x for x in gff_src]
        if (gff_entry):
            m_start = mre[0].start
            m_end = mre[0].end
            gff_start = gff_entry[0][0]
            gff_end = gff_entry[0][1]

            if (m_start >= gff_start) and (m_end <= gff_end):
                yield mre[0]
            else:
                pass
    mre_file.close()
    feature_file.close()


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--method", dest="method", type="string",
                      help="method to filter MREs")

    parser.add_option("--filter-filename", dest="filter_file", type="string",
                      help="input filter file path")

    parser.add_option("--task", dest="task", type="string",
                      help="analysis task to be executed")

    parser.add_option("--annotation-gtf-file", dest="annot_gtf", type="string",
                      help="GTF file containing transcripts for MRE "
                      "annotation")

    # add common options (-h/--help, ...) and parse command line

    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    # Write footer and output benchmark information.
    E.Stop()

    if options.task == "filter" and options.method == "list":
        filter_set = set()
        filter_file = IOTools.openFile(options.filter_file, "rb")

        for each in filter_file.readlines():
            filter_set.add(each.rstrip("\n"))
        filter_file.close()

        out_gen = filterMREsTSV(input_file=infile,
                                filter_set=filter_set)

        for x in out_gen:
            options.stdout.write("%s\n" % x)

    elif options.task == "filter" and options.method == "gtf":
        for each in filterMREsGTF(input_file=infile,
                                  filter_file=options.filter_file):
            options.stdout.write("%s\n" % each)

    elif options.task == "annotate":
        # annotate GTF of MREs with target transcript information
        transcript_gtf = options.annot_gtf
        for mre in P36.annotateMreGTF(lnc_gtf=transcript_gtf,
                                      mre_gtf=infile):
            options.stdout.write("%s\n" % mre)
    else:
        pass

if __name__ == "__main__":
    sys.exit(main(sys.argv))
