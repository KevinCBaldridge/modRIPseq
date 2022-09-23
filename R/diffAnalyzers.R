
#thisfunction, you can optionally pass the BPPARAM option with your parameter set, or parallel=TRUE if you registered the params, if you desire to use parallel processing
#is it worth writing a dedicated function for looping over dds/ddr obj list named object entries? then that could be called,
#make sure to add error-check for length of named list parameters?

#'@title Run DESeq2 analysis in the modRIPseq pipeline
#'@description This function will run the DESeq analysis of the DESeqDataSet objects prepared in the modRIPseq pipeline to provide the resulting analyzed data.
#'@details This function provides a wrapper to the DESeq() function that analyzes DESeqDataSet objects.
#'  The test used in this wrapper is hard-coded as the likelihood-ratio test, comparing a reduced model and a full model
#'  The formula list parameter requires a specific structure, wherein for each name, the corresponding value is a vector with the first element as the reduced model and the second element as the full model in mosaic notation, i.e. formulaList[['IPc']]<- c('~Replicate','~Replicate+AbTreatment').
#'@export
#'@param ddsobjlist The named list of DESeq DataSet objects that have been created already by one of the DESeqDataSetFrom*** functions. Defaults to ddsObjList
#'@param formulalist The named list of vectors containing full and reduced formulas. The length of this vector should be . Defaults to formulaList
#'@param ... optionally, can pass named keyword parameters such as BPPARAM or parallel=TRUE for using biocparallel. for parallel=TRUE, you need to have registered the params already in your environment
#'@return named list of DESeq results objects that contain the DESeq2 analyzed DESeqDataSet objects
#'@examples
#'runDEseqSet(ddsObjList,formulaList)
#'runDEseqSet(ddsObjList,formulaList,parallel=TRUE)
#'runDEseqSet(myDDSobjList,myFormulaList)
#'@seealso [modRIPseq::splitAndBuildDEobjs()] which is the step before this function in the 8OGseq pipeline
#'@seealso [DESeq2::DESeq()] which this function wraps
#'@seealso [BiocParallel] which can be used to employ parallel processing with DESeq2, and therefore with this package
#'@seealso [BiocParallel::register()] which must be used before passing parallel=TRUE parameter to this function
runDEseqSet <- function(ddsobjlist=ddsObjList, ...){
  ddsDEobjList <- c()
  ddrObjList <- c()
  for (i in 1:length(ddsobjlist)){
    NAME <- names(ddsobjlist[i])
    if (! str_detect(NAME,"In")){
    ddsDEobjList[[NAME]] <- DESeq2::DESeq(ddsobjlist[[NAME]],
                                  test='LRT',
                                  full=~replicate+abTreatment,
                                  reduced=~replicate,
                                  ...
                                  )

    }
    else {
      ddsDEobjList[[NAME]] <- DESeq2::DESeq(ddsobjlist[[NAME]],
                                    ...
                                    )
    }
    ddrObjList[[NAME]] <- DESeq2::results(ddsDEobjList[[NAME]],...)
  }
  return(ddrObjList)
}

