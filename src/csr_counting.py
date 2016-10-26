'''
csr_counting.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Count reads over class switched immunoglobulin heavy chain transcripts

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
import pysam
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import os


def gtfPositions(infile):
    '''
    collect all genomic positions across a geneset
    '''

    pos = set()
    with IOTools.openFile(infile, "rb") as open_file:
        gtf_it = GTF.transcript_iterator(GTF.iterator(open_file))
        for it in gtf_it:
            for trans in it:
                refs = range(trans.start - 1,
                             trans.end + 1)
                pos.update(refs)

    return pos


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

    parser.add_option("--j-gtf", dest="j_gtf", type="string",
                      help="gtf file of IgH J genes")

    parser.add_option("--c-gtf", dest="c_gtf", type="string",
                      help="gtf file of IgH constant genes")

    parser.add_option("--ig-coordinates", dest="locus", type="string",
                      help="reference coordinates for IgH locus in "
                      "aligned genome")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    # check for bam file index, otherwise make one
    try:
        assert os.path.exists(infile + ".bai")
    except AssertionError:
        E.info("No index file exists for %s" % infile)
        E.info("generating index for bam file %s" % infile)
        os.system("samtools index %s" % infile)

    j_pos = gtfPositions(options.j_gtf)
    c_pos = gtfPositions(options.c_gtf)

    coords = options.locus.split(":")
    contig = coords[0]
    start = int(coords[-1].split("-")[0])
    end = int(coords[-1].split("-")[1])

    samfile = pysam.AlignmentFile(infile, "rb")
    read_cache = {}
    for read in samfile.fetch(reference=contig,
                              start=start,
                              end=end):
        if read.is_proper_pair:
            try:
                read_cache[read.query_name].append(read)
            except KeyError:
                read_cache[read.query_name] = []
                read_cache[read.query_name].append(read)
        else:
            pass

    switched_reads = set()
    for pair in read_cache.values():
        # select unique pairs
        if len(pair) > 2:
            pass
        elif len(pair) == 2:
            read1 = pair[0]
            read2 = pair[1]
            r1_ref = set(read1.get_reference_positions())
            r2_ref = set(read2.get_reference_positions())
            if (r1_ref.intersection(j_pos)) and (r2_ref.intersection(c_pos)):
                switched_reads.add((read1, read2))
            elif (r1_ref.intersection(c_pos)) and (r2_ref.intersection(j_pos)):
                switched_reads.add((read1, read2))
        else:
            pass

    # add Ig constant genes to dictionary
    ig_dict = {}
    with IOTools.openFile(options.c_gtf) as gfile:
        for gene in GTF.transcript_iterator(GTF.iterator(gfile)):
            for trans in gene:
                pos = set(range(trans.start, trans.end))
                symbol = trans.asDict()['transcript_name'].split("-")[0]
                try:
                    ig_dict[symbol].update(pos)
                except KeyError:
                    ig_dict[symbol] = pos

    ig_count = {}
    for pair in switched_reads:
        all_refs = set()
        all_refs.update(pair[0].get_reference_positions())
        all_refs.update(pair[1].get_reference_positions())
        for gene in ig_dict.keys():
            inter = len(ig_dict[gene].intersection(all_refs))
            try:
                if inter:
                    ig_count[gene] += 1
                else:
                    ig_count[gene] += 0
            except KeyError:
                if inter:
                    ig_count[gene] = 1
                else:
                    ig_count[gene] = 0

    options.stdout.write("Ig_isotype\tcounts\n")
    for each in ig_count.keys():
        options.stdout.write("%s\t%i\n" % (each,
                                           ig_count[each]))
    # cache reads based on read name

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
