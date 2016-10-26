################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
'''
cgat_script_template.py
=============================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Generates a CORT distance matrix using DTW with a user definded
value of `k` for the adaptive tuning function.

Input is a single time-series expression data set with no replicates.
The input file must be tab-delimited, and time points must be in order.

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
import math
import itertools
import pandas as pd
import numpy as np
from rpy2.robjects import r as R
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

import CGAT.Experiment as E


def CORT(series1, series2):
    '''
    Calculate the temporal correlation according to Chouakira & Nagabhushan
    Assumes both time series are of the same length
    '''

    series1 = list(series1)
    series2 = list(series2)
        
    sum_prod = []
    sum_usq = []
    sum_vsq = []
    for i in range(len(series1)-1):
        u = float(series1[i+1]) - float(series1[i])
        v = float(series2[i+1]) - float(series2[i])
        prod = u * v
        sum_prod.append(prod)
        sq_u = u**2
        sq_v = v**2
        sum_usq.append(sq_u)
        sum_vsq.append(sq_v)
    
    nume = sum(sum_prod)
    denom = math.sqrt(sum(sum_usq)) * math.sqrt(sum(sum_vsq))
    
    if denom != 0:
        return(nume/float(denom))
    else:
        return 0


def adaptiveTune(value, k):
    '''
    Calculate the adaptive tuning function from Chouakira & Nagabhushan
    '''
    
    return (2/(1 + math.exp(k*abs(value))))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--k", dest="k", type="int",
                      help="value of k to adjust adaptive tuning function")
    
    parser.add_option("--out", dest="outfile", type="string",
                      help="output file name")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    DTW = importr("dtw")
    infile = args[0]

    data = pd.read_csv(infile, sep="\t",
                       index_col=0,
                       header=0)
    data.convert_objects(convert_numeric=True)

    # data should already be sorted in time-series order
    # the time and replicate columns needs to be dropped to ensure only the
    # gene data is passed into the DTW function
    try:
        data.drop(['Unnamed: 0'], inplace=True)
    except ValueError:
        data.drop(['gene_id'], inplace=True)
    data.drop(['times'], inplace=True)
    data.drop(['replicates'], inplace=True)
    genes = data.index
   
    # create a data frame of zeros of size number of genes x number of genes
    # fill it with the calculated distance metric for each pair wise comparison

    df_ = pd.DataFrame(index=genes,
                       columns=genes).fillna(0).astype(np.float64)

    # fill the array with dtw-distance values
    pandas2ri.activate()

    # iterate over the genes list in nested loops to get
    # all pair-wise combinations.

    def dtw_wrapper(data, i, j):
        series1 = data.loc[i].tolist()
        series2 = data.loc[j].tolist()
        DTW_value = (R.dtw(series1,
                           series2)).rx('distance')[0][0]
        cort_value = CORT(series1, series2)
        tuned_value = adaptiveTune(cort_value, k=int(options.k))
        time_dist = DTW_value * tuned_value
        df_[i][j] = float(time_dist)
        df_[j][i] = float(time_dist)

    for i in genes:
        for j in genes:
            dtw_wrapper(data, i, j)

    df_.to_csv(options.outfile, sep="\t")
    
    # write footer and output benchmark information
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

