
def DESeqNormalize(infile, outfile):

    # generates a lists for the design data frame
    # of the proper length
    # these need to be rpy2 objects to be parsed
    # properly in the string formatting

    time_rep_comb = [x for x in itertools.product(time_points, replicates)]
    time_cond = ro.StrVector([x[0] for x in time_rep_comb])
    rep_cond = ro.StrVector([x[1] for x in time_rep_comb])
    
    R('''countsTable = read.delim('%(infile)s',
    header=T,
    row.names=1)''' % locals())

    R('''design <- data.frame(row.names = colnames(countsTable), condition=%s, replicates=%s)''' % (time_cond.r_repr(),
                                                                                                   rep_cond.r_repr()))

    # create the count data set and normalize to library size
    # transform with variance stabilizing transformation
    # only select genes with an average of ten reads mapping

    R('''cds <- newCountDataSet(countsTable, design)''')
    R('''cds_size <- estimateSizeFactors(cds)''')
    R('''cds_disp <- estimateDispersions(cds_size, method="blind")''')
    R('''notZero <- (rowMeans(counts(cds))>10)''')
    R('''vst <- varianceStabilizingTransformation(cds_disp)''')
    R('''vst <- vst[notZero, ]''')
    # format data set to long format with condition and replicate labels
    # convert to a numpy array
    R('''replicates <- c(%s)''' % rep_cond.r_repr())
    R('''times <- c(%s)''' % time_cond.r_repr())
    R('''trans_vst = data.frame(t(exprs(vst)), times, replicates)''')
    
    # output the expression data as a tab-delimited file
    
    R('''write.table(trans_vst, file='%(outfile)s', sep="\t")''' % locals())

def diffExpression(infile, outfile):
    ref_gtf = str(infile).split("-")[1]
    data_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    design_dict = {}
    for x in data_frame.index.values:
        sample_dict = {}
        condition = str(x).split(".")[0]
        sample_dict[condition] = 1
        sample_dict['times'] = int(str(x).split(".")[1])
        sample_dict['replicates'] = str(x).split(".")[2]
        design_dict[x] = sample_dict
    design_frame = pd.DataFrame(design_dict)
    design_frame = design_frame.T
    cols = ['times', 'replicates', condition]
    design_frame = design_frame[cols]
    design_file = "deseq.dir/%s-%s-design.tsv" % (condition, ref_gtf)
    design_frame.to_csv(design_file, sep="\t")
    data_file = "deseq.dir/%s-%s-data.tsv" % (condition, ref_gtf)
    results_file = "deseq.dir/%s-%s-maSigPro.tsv" % (condition, ref_gtf)

    # data frame columns must be in the order time-replicate-condition
    # for maSigPro
    # define the numnber of higher-order terms included in the models

