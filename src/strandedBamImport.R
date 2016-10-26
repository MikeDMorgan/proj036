strandedBamImport = function (file, selection) {
  # courtesy of @sidderb - https://gist.github.com/sidderb/e485c634b386c115b2ef
  require(Rsamtools)
  
  if (!file.exists(paste(file, "bai", sep = "."))) {
    stop("Unable to find index")
  }
  
  # get pos strand pairs (F2R1):
  param_f2 = ScanBamParam(what = c("pos", "qwidth", "strand"), 
                          which = selection, 
                          flag = scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isFirstMateRead = FALSE, isMinusStrand = FALSE))
  x_f2 = scanBam(file, param = param_f2)[[1]]
  gr_f2 = GRanges(strand=x_f2[["strand"]], 
                  ranges=IRanges(x_f2[["pos"]], width = x_f2[["qwidth"]]), seqnames=seqnames(selection)[1])
  
  param_r1 = ScanBamParam(what = c("pos", "qwidth", "strand"), 
                          which = selection, 
                          flag = scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isFirstMateRead = TRUE, isMinusStrand = TRUE))
  x_r1 = scanBam(file, param = param_r1)[[1]]
  gr_r1 = GRanges(strand=x_r1[["strand"]], 
                  ranges=IRanges(x_r1[["pos"]], width = x_r1[["qwidth"]]), seqnames=seqnames(selection)[1])
  
  gr_f2r1 = c(gr_f2,gr_r1)
  
  # get rev strand reads (F1R2):
  param_f1 = ScanBamParam(what = c("pos", "qwidth", "strand"), 
                          which = selection, 
                          flag = scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isFirstMateRead = TRUE, isMinusStrand = FALSE))
  x_f1 = scanBam(file, param = param_f1)[[1]]
  gr_f1 = GRanges(strand=x_f1[["strand"]], 
                  ranges=IRanges(x_f1[["pos"]], width = x_f1[["qwidth"]]), seqnames=seqnames(selection)[1])
  
  param_r2 = ScanBamParam(what = c("pos", "qwidth", "strand"), 
                          which = selection, flag = scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isFirstMateRead = FALSE, isMinusStrand = TRUE))
  x_r2 = scanBam(file, param = param_r2)[[1]]
  gr_r2 = GRanges(strand=x_r2[["strand"]], 
                  ranges=IRanges(x_r2[["pos"]], width = x_r2[["qwidth"]]), seqnames=seqnames(selection)[1])
  
  gr_f1r2 = c(gr_f1,gr_r2)
  
  
  # calc coverage on both strands:
  cov_list = list("Forward" = coverage(ranges(gr_f2r1),
                                       width=end(selection)),
                  "Reverse" = coverage(ranges(gr_f1r2), width=end(selection)))
  pos = sort(unique(unlist(lapply(cov_list, function(y) c(start(y), end(y))))))
  
  # build final GR
  stranded_cov_gr = GRanges(seqnames = seqnames(selection)[1], ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)),
                            plus=as.numeric(cov_list[["Forward"]][head(pos, -1)]),
                            minus=-as.numeric(cov_list[["Reverse"]][head(pos, -1)]))
  
  return(stranded_cov_gr)
}