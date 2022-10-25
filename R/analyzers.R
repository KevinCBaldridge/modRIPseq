#'@title Run DESeq2 analysis in the modRIPseq pipeline
#'@description This function will run the DESeq analysis of the DESeqDataSet objects prepared in the modRIPseq pipeline to provide the resulting analyzed data.
#'@details This function provides a wrapper to the DESeq() function that analyzes DESeqDataSet objects.
#'  The test used in this wrapper is hard-coded as the likelihood-ratio test, comparing a reduced model and a full model
#'  The formula list parameter requires a specific structure, wherein for each name, the corresponding value is a vector with the first element as the reduced model and the second element as the full model in mosaic notation, i.e. formulaList[['IPc']]<- c('~Replicate','~Replicate+AbTreatment').
#'@export
#'@param ddsobjlist The named list of DESeq DataSet objects that have been created already by one of the DESeqDataSetFrom*** functions. Defaults to ddsObjList
#'@param formulalist The named list of vectors containing full and reduced formulas. The length of this vector should be . Defaults to formulaList
#'@param ... optionally, can pass named keyword parameters such as DESeq2 options passed to [DESeq2::DESeq()] or options like BPPARAM or parallel=TRUE for using biocparallel. for parallel=TRUE, you need to have registered the params already in your environment
#'@return named list of DESeq results objects that contain the DESeq2 analyzed DESeqDataSet objects
#'@examples
#'runDEseqSet(ddsObjList,formulaList)
#'runDEseqSet(ddsObjList,formulaList,parallel=TRUE)
#'runDEseqSet(myDDSobjList,myFormulaList)
#'runDEseqSet(ddsObjList,formulaList,test="Wald",fitType='local')
#'@seealso [modRIPseq::splitAndBuildDEobjs()] which is the step before this function in the 8OGseq pipeline
#'@seealso [DESeq2::DESeq()] which this function wraps
#'@seealso [BiocParallel] which can be used to employ parallel processing with DESeq2, and therefore with this package
#'@seealso [BiocParallel::register()] which must be used before passing parallel=TRUE parameter to this function
runDEseqSet <- function(ddsobjlist=ddsObjList, ...){
  ddsDEobjList <- c()
  ddrObjList <- c()
  for (i in 1:length(ddsobjlist)){
    NAME <- names(ddsobjlist[i])

    ddsDEobjList[[NAME]] <- DESeq2::DESeq(ddsobjlist[[NAME]],
                                  ...
                                  )
    ddrObjList[[NAME]] <- DESeq2::results(ddsDEobjList[[NAME]],...)
  }
  return(ddrObjList)
}

#'@title Sort the DESeqResults object by column of choice
#'@description Arrange the DESeqResults object by values in the chosen column(s)
#'@details This function will operate on a single DESeqResults object to sort the results by the chosen column(s).
#'  It is suggested to sort on padj, the adjusted p-value, with secondary sorting on the log2FoldChange from high to low.
#'  Note that you should loop over your ddrObjList and pass each list entry therein to this function.
#'@param res the DESeqResults object you wish to sort. Defaults to ddr
#'@param colSort the column you would like to sort on, will be passed to [dplyr::arrange] and so may be used as descending, i.e. desc(colSort). Defaults to padj
#'@param ... Additional arguments to pass into [dplyr::arrange] for sorting the df, i.e. desc(log2FoldChange). No default value
#'@examples
#'ddrSorted <- sortRes(ddr,padj)
#'sortRes(myddr,padj,desc(log2FoldChange))
#'for (NAME in names(ddrObjList)){ddrObjList[[NAME]]<-sortRes(ddrObjList[[NAME]])}
#'@export
#'@return DESeqResults object that has been sorted according to chosen columns
#'@seealso [dplyr::arrange()] that this method wraps
#'@seealso [[modRIPseq::runDEseqSet()]] that generates the list of DESeqResults objects as shown in the examples
#'@seealso [[DESeq2::DESeqResults]] that is the S4 object class this function will operate on
sortRes <- function(res=ddr,colSort=padj,...){
  #add error check with explanation if passed arguments are not columns of the DESeqResults object?
 tmpdf <- data.frame(res@listData)
 rownames(tmpdf) <- res@rownames
 tmpdf <- tmpdf %>% arrange({{colSort}},...)
 for (NAME in names(tmpdf)){
   res@listData[[NAME]]<-purrr::as_vector(tmpdf[[NAME]])
   }
 res@rownames <- rownames(tmpdf)
 return(res)
}



#'@title Filter top N most significant diff transcripts of DESeqResults object
#'@description subset the top N entries in the results of your modRIPseq analysis
#'@details this function cuts your results to the top N entries after sorting
#'  Run [modRIPseq::sortRes()] to sort on your desired criteria before running this function to subset
#'  Note that you should loop over your ddrObjList and pass each list entry therein to this function.
#'@param res the DESeqResults object from which to subset the top entries. Defaults to ddr
#'@param num the number of top entries to subset from the results object. Defaults to 10% of the transcripts that have padj<0.1
#'@export
#'@return the top N subsetted DESeqResults object
#'@examples
#'ddrTop10pctSig <- filterTopNpadj(ddr,0.1*nrow(ddr[!is.na(ddr$padj)&ddr$padj<0.1,]))
#'filterTopNpadj(myddr,0.05*nrow(myddr[!is.na(myddr$padj)&myddr$padj<0.1,]))
#'for (NAME in names(ddrObjList)){ddrObjList[[NAME]]<-filterTopNpadj(ddrObjList[[NAME]])}
#'@seealso [modRIPseq::sortRes()] which should be run before this function to arrange the results according to your desired criteria, i.e. padj or log2FoldChange
#'@seealso [modRIPseq::runDEseqSet()] that generates the list of DESeqResults objects as shown in the examples
#'@seealso [DESeq2::DESeqResults] that is the S4 object class this function will operate on
filterTopNpadj <- function(res=ddr, num=0.1*nrow(res[!is.na(res$padj)&res$padj<0.1,])){
  filtered <- head(res,num)
  return(filtered)
}
