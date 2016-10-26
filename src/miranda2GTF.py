'''
miranda2GTF.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------
Parse output from miRNA binding site prediction softward
miRanda - a convoluted format and can produce
very large output files (~ GBs).

Improvements/to do:
:: Needs to either parse it in chunks or
set aside a specific portion of memory
to hold the entire file.

:: PEP8 compliant

Usage
-----

Example::

   python miranda2GTF.py

Type::

   python cgat_script_templatemiranda2GTF.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import itertools
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools


def attributeFinder(line, regex, delim):
    '''
    Find a term line by line from a fileSplitter generator
    '''
        
    for x in line:
        if re.search(regex, x):
            return getAttribute(x, regex, delim)


def getAttribute(attributeList, attribute, delim):
    '''
    return the particular attribute value split by the known delimiter
    e.g. gene_id="Gm12838" would be split on `=` and return the value
    "Gm12838"
    '''
    
    for j in attributeList.split(delim):
        if re.search(attribute,  j):
            if len(j.split(delim)[-1]) != 0:
                return j.split(" ")[-1]
        else:
            pass


def mirandaParse(infile, align, percent):
    '''
    A method to handle and parse output files from the miRNA binding site
    prediction software miRanda.
    
    The first 22 lines are header and comments.
    Lines 23-32 are the run settings for that execution of miRanda.
    
    
    A new miRNA record begins with:
    
    Read Sequence:ID=`miRNA_id`;Alias=`miRNA_alias`;Name=`miRNA_name`;
    Derives_from=`miRNA_derived` `contig`:`strand`:`start`-`end`(`length` nt)

    The alignments with the query sequence start thus:
    
    Read Sequence:gene_id `gene_id`; transcript_id `transcript_id`; `contig`:
    `strand`:`start`-`end`(`length` nt)
    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    Performing Scan: ID=`miRNA_id`;Alias=`miRNA_alias`;Name=`miRNA_name`;
    Derives_from=`miRNA_derived` vs gene_id
    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    '''

    line_counter = 0
    for line in infile:
        line_counter += 1
        if line_counter <= 32:
            pass

        elif line.split(" ")[0] == "Performing":
            miRNA_name = line.split("=")[3].rstrip(";Derives_from")
                    
        elif line.split(":")[0] == "Read Sequence":
            if line.split(":")[1].split(" ")[0] == "exon_number":
                exon_number = line.split(';')[0].split(":")[1].split(" ")[1]
                exon_number = exon_number.lstrip('"').rstrip('"')
                exon_status = line.split(';')[1].lstrip(" exon_status ")
                exon_status = exon_status.rstrip('"').lstrip('"')
                target_id = line.split(';')[3].lstrip(" gene_id ")
                target_id = target_id.lstrip('"').rstrip('"')
                contig = line.split(';')[-1].split(":")[0]
                contig = contig.lstrip(" ")
                strand = line.split(';')[-1].split(":")[1]
                target_start = line.split("; ")[-1].split(":")[2]
                target_start = int(target_start.split("-")[0])
            else:
                pass

        elif line.split("=")[0] == ">ID":
            align_length = int(line.split("\t")[4].split(" ")[1]) - int(line.split("\t")[4].split(" ")[0])
            align_pc = float(line.split("\t")[-1].rstrip("%\n"))
            start = int(line.split("\t")[5].split(" ")[0]) + target_start
            end = start + align_length
            
            entry = GTF.Entry()
                
            entry.contig = contig
            entry.start = start
            entry.end = end
            entry.strand = strand
            entry.feature = "MRE"
            entry.score = "."
            entry.frame = "."
            entry.gene_id = "%s_%s:%i-%i" % (miRNA_name,
                                             contig,
                                             start,
                                             end)
            entry.transcript_id = "%s_%s:%i-%i" % (target_id,
                                                   contig,
                                                   start,
                                                   end)
            entry.addAttribute("target", target_id)
            entry.addAttribute("miRNA", miRNA_name)
            entry.addAttribute("align_pc", align_pc)
            entry.addAttribute("exon_number", exon_number)
            entry.addAttribute("exon_status", exon_status)

            if align_length > align and align_pc > percent:
                yield entry
            else:
                pass


def targetScanParse(infile, lnc_gtf):
    '''
    Parse results from targetScan into GTF
    '''
    gtf_dict = {}
    lnc_file = IOTools.openFile(lnc_gtf)
    for each in GTF.transcript_iterator(GTF.iterator(lnc_file)):
        for trans in each:
            entry = GTF.Entry()
            entry = entry.copy(trans)
            gtf_dict[entry.transcript_id] = entry

    lnc_file.close()
    counter = 0
    for line in infile:
        counter += 1
        line = line.split("\t")
        if counter > 1:
            MRE = GTF.Entry()
            gene_id = line[0].lstrip('"').rstrip('"')
            target = gtf_dict[gene_id]
            align_start = int(line[3])
            align_end = int(line[4])
            size = align_end - align_start
            miRNA = "mmu-%s" % line[1]
            seed_class = line[8]

            MRE.contig = target.contig
            MRE.feature = "MRE"
            MRE.start = target.start + align_start
            MRE.end = MRE.start + size
            MRE.source = target.source
            MRE.strand = target.strand
            MRE.addAttribute('miRNA', miRNA)
            MRE.addAttribute('target', gene_id)
            try:
                MRE.addAttribute('exon_number',
                                 target.asDict()['exon_number'])
            except KeyError:
                E.info("No exon number data in GTF for %s" % gene_id)
                MRE.addAttribute('exon_number', '.')

            if target.source == "protein_coding":
                MRE.addAttribute('exon_status', "protein_coding")
            else:
                try:
                    MRE.addAttribute('exon_status',
                                     target.asDict()['exon_status'])
                except KeyError:
                    E.info("No exon status data in GTF for  %s" % gene_id)
                    MRE.addAttribute('exon_status', '.')

            MRE.transcript_id = "%s_%s:%i-%i" % (gene_id,
                                                 MRE.contig,
                                                 MRE.start,
                                                 MRE.end)
            MRE.gene_id = "%s_%s:%i-%i" % (miRNA,
                                           MRE.contig,
                                           MRE.start,
                                           MRE.end)
            MRE.addAttribute('seed_class', seed_class)

            yield MRE


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

    parser.add_option("-a", "--align-length", dest="align", type="string",
                      help="alignment length minimum")

    parser.add_option("-s", "--align-score", dest="score", type="string",
                      help="percent nts aligned threshold")

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("miranda", "targetscan"), help="program "
                      "used to predict MREs")

    parser.add_option("--GTF", dest="gtf", type="string",
                      help="GTF file used to generate targets for targetScan")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = options.stdin

    parser.set_defaults(align=6,
                        score=80.0)

    E.info("converting to GTF format")

    if options.format == "miranda":

        length_threshold = int(options.align)
        percent_threshold = float(options.score)
        
        for x in mirandaParse(infile=infile,
                              align=length_threshold,
                              percent=percent_threshold):
            options.stdout.write("%s\n" % x)

    elif options.format == "targetscan":
        lnc_gtf = options.gtf
        for x in targetScanParse(infile=infile,
                                 lnc_gtf=lnc_gtf):
            options.stdout.write("%s\n" % x)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
