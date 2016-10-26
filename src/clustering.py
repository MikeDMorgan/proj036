'''
clustering.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Provide consensus clustering quality metrics

Options
-------

.. Options for the script, with detail of how they are combined
.. to provide intended functionality

Usage
-----

.. Example use case

Example::

   python clustering.py

Type::

   python clustering.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import pandas as pd
import itertools
import re
import os


def clusterConcordia(data1, data2, genes):
    '''
    Calculate clustering agreement counts:
    a = number of pairs in the same cluster in clustering 1 and clustering 2
    b = number of pairs in different clusters in clustering 1 and clustering 2
    c = number of pairs of clusters the same in clustering 1 and
    differenent in clustering 2
    d = number of pairs of clusters different in clustering 1 and the same
    in clustering 2
    Uses sets as these are much faster and more efficient than lists
    '''

    # set up all possible gene-pairs and set containers

    comp_iter = itertools.combinations_with_replacement(genes, 2)
    complete_pairs = set([x for x in comp_iter])
    
    a = set()
    b = set()
    c = set()
    d = set()
    total = set()
    
    # iterate over each pair-wise combination of cluster assignments
    for key1 in data1.keys():
        for key2 in data2.keys():
            set1 = data1[key1]
            set2 = data2[key2]
            c.update(set1.difference(set2))
            d.update(set2.difference(set1))
            a.update(set2.intersection(set1))
            total.update(set1.union(set2))
    b.update(complete_pairs.difference(total))
    
    concordance = {}
    concordance['a'] = len(a)
    concordance['b'] = len(b)
    concordance['c'] = len(c)
    concordance['d'] = len(d)

    return concordance


def concordanceMetric(concord_dict):
    '''
    Calculate cluster quality and concordance metrics:
    Rand Index, Adjusted Rand Inded, F-measure, Jaccard Index
    '''

    metric_dict = {}
    a = concord_dict['a']
    b = concord_dict['b']
    c = concord_dict['c']
    d = concord_dict['d']

    Rand = (a + b) / float(a + b + c + d)
    AdjRand = 2 * ((a*b) - (c+d)) / float(((a+d)*(d+b)) + ((a+c)*(c+b)))
    pos = a/float(a+c)
    recall = a/float(a+d)
    beta = 1
    Fstat = ((beta**2 + 1)*(pos*recall)) / float(((beta**2) * pos) + recall)
    Jaccard = a / float(a + c + d)
    Total = sum([a, b, c, d])

    metric_dict['Rand'] = Rand
    metric_dict['AdjRand'] = AdjRand
    metric_dict['Positive'] = pos
    metric_dict['Recall'] = recall
    metric_dict['F_measure'] = Fstat
    metric_dict['Jaccard'] = Jaccard
    metric_dict['Total'] = Total

    return metric_dict


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
    
    parser.add_option("--directory", dest="sample_dir", type="string",
                      help="directory containing clustering files")

    parser.add_option("--filename-regex", dest="file_reg", type="string",
                      help="regex for filenames from the same samples")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    files = os.listdir(options.sample_dir)
    assert files
    pattern = r"%s" % options.file_reg
    file_match = [x for x in files if re.search(pattern, x)]
    file_combos = itertools.combinations(file_match, 2)
    results_dict = {}

    for file_pair in file_combos:
        infile1 = "%s/%s" % (options.sample_dir, file_pair[0])
        infile2 = "%s/%s" % (options.sample_dir, file_pair[1])
        
        data1 = pd.read_table(infile1, sep=" ", index_col=0, header=0)
        data2 = pd.read_table(infile2, sep=" ", index_col=0, header=0)
        
        genes = data1.index.tolist()
        data1.columns = ['cluster']
        data2.columns = ['cluster']

        # make a set for each cluster that contains all gene-pairs in that
        # cluster test the intersection of each cluster-wise comparison for
        # concordance between two clusterings

        data1_groups = data1.groupby(by='cluster')
        data2_groups = data2.groupby(by='cluster')

        cl_dict1 = {}
        cl_dict2 = {}

        for name, group in data1_groups:
            cl_dict1[name] = group.index.tolist()

        for name, group in data2_groups:
            cl_dict2[name] = group.index.tolist()

        cluster_sets1 = {}
        for key in cl_dict1.keys():
            cl_g = itertools.combinations_with_replacement(cl_dict1[key], 2)
            cl_set = set([x for x in cl_g])
            cluster_sets1[key] = cl_set

        cluster_sets2 = {}
        for key in cl_dict2.keys():
            cl_g = itertools.combinations_with_replacement(cl_dict2[key], 2)
            cl_set = set([x for x in cl_g])
            cluster_sets2[key] = cl_set
        
        concordance = clusterConcordia(cluster_sets1,
                                       cluster_sets2,
                                       genes)
        metrics = concordanceMetric(concordance)
        results_dict[file_pair] = metrics

    res_frame = pd.DataFrame(results_dict).T
    
    res_frame.to_csv(options.stdout,
                     sep="\t",
                     index_label='cluster_pair')

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
