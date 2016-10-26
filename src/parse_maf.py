
'''
parse_maf.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.Pipeline as P

class seqRecord(object):
    '''
    A class to hold individual maf alignment entries
    '''

    def __init__(self):
        self.species = ""
        self.source = ""
        self.start = 0
        self.size = 0
        self.source_size = 0
        self.sequence = ""
def maf_parse(maf_file,
              fasta_dir,
              species_list,
              ref_species):
              
    '''
    Output alignments to fasta format and conservation scores
    for the reference species to bed format
    '''

    scores = open(maf_file.rstrip(".maf") + (".bed"), "w")
    scores.write('track name="multiz60wayScores" description="Multiz60way conservation scores"\n')
    align = []
    with open(maf_file, "r") as infile:
        E.info("Converting maf alignment entries to fasta format")
        fasta_counter = 0
        for line in infile:
            line = line.split(" ")
            if line[0] == 'a':
                score = line[1].split("=")[1]
                if len(align) > 0:
                    fasta_file = open("%s/%s-%i.fa" %(fasta_dir,
                                                      maf_file.rstrip(".maf"),
                                                      fasta_counter), "w")
                    for entry in align:
                        fasta_file.write(">%s_%s_%s_%s\n" % (entry.species,
                                                             entry.source,
                                                             entry.start,
                                                             entry.size))
                        fasta_file.write("%s\n" % entry.sequence)
                    fasta_file.close()
                    fasta_counter += 1
                    align = []
                    if entry.species == ref_species:
                        scores.write("%s\t%s\t%i\t%s\n" % (entry.source,
                                                           entry.start,
                                                           int(entry.start)+int(entry.size)-1,
                                                           score.rstrip("\n")))
                else:
                    pass
            elif line[0] == 's':
                line[:] = [x for x in line if x != '']
                entry = seqRecord()
                entry.species = line[1].split(".")[0]
                if entry.species in species_list:
                    entry.source = line[1].split(".")[1]
                    entry.start = int(line[2]) # maf files are 0-based like bed files
                    entry.size = line[3]
                    entry.strand = line[4]
                    entry.source_size = line[5]
                    entry.sequence = ("".join(line[-1].split("-"))).rstrip("\n")
                    align.append(entry)
                else:
                    pass

    E.info("Conservation scores written to %s" % scores.name)
    scores.close()
    return "Fasta files created: %i" % fasta_counter
    

def gtf2_maf(infile,
             maf_dir):
    '''
    pull out all of the maf alignments for a given input
    gtf file and merge them into a single output file
    '''
    
    maf_pattern = r"(.+).maf$"
    chroms = []
    
    # regex for whole chromosome contigs, discounting
    # unnumbered chromosomes and unassigned contigs/scaffolds
    
    for x in os.listdir(maf_dir):
        match = re.search(maf_pattern, x)
        if match != None:
            [chroms.append(str(x)) for x in match.groups()]
        else:
            pass
    
    maf_files = ["%s/%s.maf" % (maf_dir, x) for x in chroms]
    suffix = "%s-%s" % (infile.split("-")[0],
                        infile.split("-")[1])

    statement_list = []
    for maf in maf_files:
        chrom = maf.split("/")[-1].rstrip(".maf")
        statement = "maf_parse --features %s %s " % (infile,
                                                     maf)       
        statement_list.append(statement)

    chain = "".join(["%s; " % x for x in statement_list])
       
    return chain
   

def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--maf-directory", dest="maf_dir", type="string",
                      help="directory containing .maf files")

    parser.add_option("--fasta-directory", dest="fasta_dir", type="string",
                      help="directory to output fasta files to")
    
    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--task", dest="task", type="string",
                      help="task to run")

    parser.add_option("--align-algorithm", dest="algorithm", type="string",
                      help="alignment method to use")

    parser.add_option("--species", dest="species",type="string",
                      help="a comma separated listed of taxonomies")
    
    parser.add_option("--ref-species", dest="ref_species", type="string",
                      help="target/reference alignment species")
    
    parser.add_option("--output-format", dest="outformat", type="string",
                      help="output format for MAFFT")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    
    infile = argv[-1]
    
    species_list = (options.species).split(",")

    if options.task == "parse":
        res = maf_parse(maf_file=infile,
                        fasta_dir=options.fasta_dir,
                        species_list=species_list,
                        ref_species=options.ref_species)

        options.stdout.write(res)
    else:
        pass
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
