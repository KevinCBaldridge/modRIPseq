#' #write to have option for high (suppl fig 10c) or low (suppl fig 9b) plot
#' #'@title Make volcano plot for modRIPseq results
#' #'@description This function will create a volcano plot for the DESeq results for modRIPseq data
#' #'@details
#' #'@param
#' #'@export
#' #'@return
#' #'@examples
#' #'@seealso
#'
#' #not in paper? but part of data QC in pre-paper workup
#' #'@title Make independent filtering plot
#' #'@description
#' #'@details
#' #'@param
#' #'@export
#' #'@return
#' #'@examples
#' #'@seealso
#'
#'
#' #suppl figure 4a - can we get
#' #'@title Make Venn diagram plot of up/down regulated/oxidized transcripts
#' #'@description
#' #'@details
#'
#'
#' #suppl figure 8
#' #'@title Make log2FC 8OG enrichment plot with cutoff
#' #'
#' #'
#'
#' #suppl fig 9a
#' #'@title
#' #'
#' #'
#'
#'
#' #make option for low (suppl fig 3a) vs high (suppl fig 10b)
#' #'@title Make venn diagram comparing high vs low level overlap
#' #'
#' #'
#'
#' #fig 3c
#' #'@title Make log2FC heatmap with inset of top N in both high vs low
#' #'
#' #'
#'
#'
#' #not in paper? but part of data QC in pre-paper workup
#' #'@title Make MA plot
#' #'
#' #'

#EnrichR plots, have it pass gene list in so it can be same function for any set of genes in analysis
#one example would be fig 2e and another is fig 3b
#'@title Make enrichr functional enrichment bar plot

#plot PCA, not in paper
#'@title Make PCA plot to assess sample spread
#'@description This function wraps the [DESeq2::plotPCA()]
#'  function to examine data spread of modRIPseq data
#'  Note that for looping over ddsTfmList object, to display the plot you'll need to wrap this in a call to [base::plot()] as in the third example
#'@details More details...
#'@param ddstfm the DESeqTransform object to be visualized on PCA plot
#'@param intgroup the intgroup option to pass into [DESeq2::plotPCA()]. if not defined, the function will set the intgroup option from the "contrast" column of the passed ddstfm object. Default is NULL
#'@export
#'@return ggplot object containing the PCA plot of passed DESeqTransform object
#'@examples
#'makePCA(ddsTfm)
#'makePCA(myDdsTfm,intgroup="myDesignCol")
#'for (NAME in names(ddsTfmList)){plot(makePCA(ddsTfmList[[NAME]]))}
#'p <- makePCA(); p + theme(plot.title = element_text(hjust = 0.5))
#'@seealso [DESeq2::plotPCA()] which this function wraps
#'@seealso [modRIPseq::ddsTfm()] which should be used to generate a DESeqTransform object for PCA plotting
makePCA <- function(ddstfm=ddsTfm,intgroup=NULL){
  if (is.null(intgroup)){
    igroup <- unique(ddstfm$contrast)
  }
  else {
    igroup <- intgroup
    print("setting intgroup from ddsTfm$contrast")
  }
  p <- DESeq2::plotPCA(ddstfm,intgroup=igroup)
  p <- p + ggtitle(NAME,subtitle =paste0("comparison: ",igroup) )
  return(p)
}
