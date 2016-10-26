'''
mreConserve.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Calculate the conservation score of MREs in a gtf file
.. and the probability the predicted score is derived from
.. a null distribution. Uses phastCons scores from UCSC in
.. bigWig format.

Options
-------

.. Options for the script, with detail of how they are combined
.. to provide intended functionality

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
import os
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
from bx.bbi.bigwig_file import BigWigFile
import numpy as np
import rpy2.robjects as ro
import itertools
import pandas as pd
import random
from Bio.AlignIO import MafIO

# require two gtf files - one of gene models and one of predicted mres
# the assumption is that the predicted mres have a `target` attritbute which
# corresponds to the lncRNA gene_id as well as the exon number to which
# it maps
# other input include a bigWig file of phastCons scores and a matching multiple
# alignment file


def keySearch(dictionary, contig, start, end):
    '''
    Search for overlap between an input search on contig, and co-ordinates
    and gtf entries in a dictionary with key: (contig, start, end).
    Return dictionary key.
    '''

    keys = (dictionary.keys())
    keys.sort()

    for i in keys:
        if (i[0] == contig):
            q_range = set(range(start, end+1))
            t_range = set(range(start, end+1))
            if q_range.intersection(t_range):
                return i
            else:
                pass
        else:
            pass

    return False


def fetchRandomMRE(MRE, bw_index, lnc_index, maf_index, species, perms):
    '''
    Return a set of pseudo-MREs from an indexed genome.
    pseudo-MREs are generated from the same exon
    to which the MRE maps
    '''

    target_cord = keySearch(lnc_index, MRE.contig, MRE.start, MRE.end)
    if target_cord:
        exon_dict = {}
        overlaps = lnc_index[target_cord]

        exon_number = MRE.asDict()['exon_number']

        start = overlaps.start
        end = overlaps.end

        mask = (MRE.start, MRE.end)
        exon_dict[exon_number] = genRandomMRE(MRE.contig,
                                              start,
                                              end,
                                              bw_index,
                                              mask,
                                              maf_index,
                                              species,
                                              perms)

        if len(exon_dict):
            return exon_dict[MRE.asDict()['exon_number']]
        else:
            pass
    else:
        raise IndexError("MRE does not match any position in this GTF file")


def genRandomMRE(contig,
                 start,
                 end,
                 bw_index,
                 mask,
                 maf_index,
                 species,
                 perms):
    '''
    Generate a set of pseudo-MREs with phastCons score
    '''

    mre_size = mask[1] - mask[0]
    mask_start = start - mask[0]
    mask_end = mask_start + mre_size

    scores = bw_index.get_as_array(contig, start, end)

    # some nucleotides have no conservation score - set these to zero
    scores[np.isnan(scores)] = 0
    size = len(scores)
    rand_scores = []
    check_maf = []
    for n in range(start, end):
        check_maf.append(checkMaf(maf_index, species, n, n+2))

    E.info("Generating null distribution of phastCons scores")
    if size <= 10:
        for i in xrange(perms):
            rand_scores.append(0.0)
    else:
        for i in xrange(perms):
            # randomly subsample dinucleotides across the length of the
            # pseudo-MRE, matched for length of predicted MRE.  Exclude 5nts
            # around exon splice sites based on Wilfried's conservation of
            # lncRNAs the MRE may take up the entire exon (minus the 5nts arond
            # the splice junctions) - set score to 0.
            nts = []
            mean_score = []
            nts = (random.randint(5, size - 5) for j in xrange(mre_size))

            # don't subsample across the same region as the predicted MRE
            for k in nts:
                if k >= mask_start and k <= mask_end:

                    # add a heuristic - if k is > half way into the MRE
                    # subsample from before, else subsample afterwards.  If k
                    # is exactly halfway subsample before unless start + 5nt
                    # is > start of MRE.

                    if (start + 5) > mask[0]:
                        k = random.randint(mask_end + 1, size - 5)
                    elif k < (mask_start + mask_end)/2.0:
                        try:
                            k = random.randint(5, mask_start - 1)
                        except(ValueError):
                            k = random.randint(mask_end + 1, size - 5)
                    score = sum(scores[k:k+1])/2.0
                    mean_score.append(score)
                else:
                    if check_maf[k]:
                        score = 0.0
                    else:
                        score = sum(scores[k:k+1])/2.0
                    mean_score.append(score)
            rand_scores.append(np.mean(mean_score))

    yield rand_scores


def fetchPhastCons(contig, start, end, bw_index, maf_index, species):
    '''
    Return the phastCons score for a given interval
    Control for gaps in multiple alignment
    '''

    E.info("Retrieving MRE phastCons score")

    score_array = bw_index.get_as_array(contig, start, end)
    score_array[np.isnan(score_array)] = 0
    mean_score = []
    for i in range(start, end):
        if checkMaf(maf_index, species, i, i+1):
            score = 0.0
        else:
            score = score_array[end - i - 1]
        mean_score.append(score)
    fin_score = np.mean(mean_score)

    return fin_score


def checkMaf(maf_index, species, start, end):
    '''
    Checks a multiple alignment file for gaps.
    If a gap exists at the queried position return 1,
    else 0.  This will cause the phastCons scores
    to be set to 0 in fetchPhastCons and genRandomMRE
    '''

    # this relies on the first species in the maf being
    # the one of interest.  If the target species
    # is not the same as `species` then use a different
    # maf file or change this to check for species id.
    
    mafs = [x for x in maf_index.search([start], [end])]
    try:
        rec = mafs[0][0]
        if rec.id.startswith(species):
            seq_start = start - rec.annotations['start']
            seq_end = seq_start + (end - start)
            for i in rec.seq[seq_start:seq_end]:
                if i == "-":
                    return 1
                else:
                    return 0
        else:
            E.warn("Check target species is compatible with input species")

    except(IndexError):
        E.warn("Region %i-%i not found in %s maf index" % (start,
                                                           end,
                                                           maf_index))
        return 0


def buildMafIndex(contig, maf_dir, species):
    '''
    Index a .maf file for a given contig and target species
    Uses BioPython MafIO module
    If mafindex exists returns the relevant index,
    otherwise it will build the index.
    '''

    E.info("Checking maf index for contig: %s" % contig)

    idx = maf_dir + "/" + contig + ".mafindex"
    species_space = species + "." + contig
    maf_path = maf_dir + "/" + contig + ".maf"

    maf_idx = MafIO.MafIndex(idx,
                             maf_path,
                             species_space)

    return maf_idx


def calcPval(pred, null):
    '''
    Based on the empirical cumulative distribution
    calculate the probability that this predicted
    MRE conservation score is derived from the
    null distribution
    '''

    E.info("Calculating p-value for conservation score")

    null_dist = itertools.chain.from_iterable(null)
    null_dist = [x for x in null_dist]
    null_dist = ro.FloatVector(null_dist)

    ecdf = ro.r('ecdf')
    null_cdf = ecdf(null_dist)
    pval = 1.0000 - null_cdf(np.float64(pred))[0]

    return pval


def MRE_results(mre_gtf, target_file, bw_file, maf_dir, species, perms):
    '''
    Iterate over GTF file, calculate phastCons and probability of
    high conservation and return as a dataframe of results
    '''

    mre_dict = {}
    trans_it = GTF.transcript_iterator(GTF.iterator(mre_gtf))
    mre_it = (x[0] for x in trans_it)

    maf_dict = {}
    maf_index = None

    target_index_dict = {}
    for entry in GTF.iterator(IOTools.openFile(target_file)):
        index_key = (entry.contig, entry.start, entry.end)
        target_index_dict[index_key] = entry

    bw_open = IOTools.openFile(bw_file)
    bw_index = BigWigFile(bw_open)

    for MRE in mre_it:
        res_dict = {}

        # if maf already exists it will create a MafIndex object,
        # else it will build the index then return a MafIndex object.
        # This is the default behaviour for MafIO.MafIndex()

        maf_idx = maf_dir + "/" + MRE.contig + ".mafindex"
        maf_index = buildMafIndex(MRE.contig,
                                  maf_dir,
                                  species)
        maf_dict[MRE.contig] = maf_index

        # check maf index exists
        try:
            assert os.path.exists(maf_idx)
        except(AssertionError):
            maf_dict[MRE.contig] = buildMafIndex(MRE.contig,
                                                 maf_dir,
                                                 species)

        phast = fetchPhastCons(MRE.contig,
                               MRE.start,
                               MRE.end,
                               bw_index,
                               maf_index,
                               species)

        null_dist = fetchRandomMRE(MRE,
                                   bw_index,
                                   target_index_dict,
                                   maf_index,
                                   species,
                                   perms)

        res_dict['target'] = MRE.asDict()['target']
        res_dict['exon'] = MRE.asDict()['exon_number']
        res_dict['exon_status'] = MRE.asDict()['exon_status']
        res_dict['source'] = MRE.source
        res_dict['contig'] = MRE.contig
        # need to add +1 correction for GTF coordinates
        res_dict['start'] = MRE.start + 1
        res_dict['end'] = MRE.end
        res_dict['size'] = MRE.end - MRE.start
        res_dict['score'] = phast
        res_dict['pCons'] = calcPval(phast, null_dist)
        mre_dict[MRE.gene_id] = res_dict

    mre_frame = pd.DataFrame(mre_dict).T

    col_order = ['contig',
                 'start',
                 'end',
                 'size',
                 'target',
                 'source',
                 'exon',
                 'exon_status',
                 'score',
                 'pCons']

    mre_frame = mre_frame[col_order]

    bw_open.close()

    E.info("Processing MRE")

    return mre_frame


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

    parser.add_option("--target-file", dest="targets", type="string",
                      help="gtf file of targets to which MREs map")

    parser.add_option("--maf-directory", dest="maf_dir", type="string",
                      help="directory containing multiple alignment files")

    parser.add_option("--cons-file", dest="cons_score", type="string",
                      help="bigWig file containig conservation scores")
    
    parser.add_option("--species-build", dest="species", type="string",
                      help="species genome build gtf, bigWig and maf"
                      "co-ordinates correspond to")

    parser.add_option("--permutations", dest="perm", type="string",
                      help="number of permutations used to calculated"
                      "null distribution of conservation scores for"
                      "each predicted MRE")

    parser.add_option("--build-index", dest="index", action="store_true",
                      help="build maf index for contig supplied in"
                      "input maf file")

    parser.add_option("--maf-file", dest="maf_file", type="string",
                      help="maf file to create index for")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    parser.set_defaults(index=False)

    if options.index:
        infile = options.maf_file
        contig = infile.split("/")[1].split(".")[0]
        maf_dir = options.maf_dir
        species = options.species
        buildMafIndex(contig,
                      maf_dir,
                      species)
    else:
        mre_file = options.stdin
        target_file = options.targets
        maf_dir = options.maf_dir
        bw_file = options.cons_score
        species = options.species
        perms = int(options.perm)

        out_frame = MRE_results(mre_file,
                                target_file,
                                bw_file,
                                maf_dir,
                                species,
                                perms)
        out_frame.index = [x.split("_")[0] for x in out_frame.index.tolist()]
                            
        out_frame.to_csv(options.stdout,
                         sep="\t",
                         header=True,
                         index_label="miRNA")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
