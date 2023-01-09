if(getRversion() >= "2.15.1") utils::globalVariables("dds")

#'@title use biomaRt to translate annotation
#'@description this function wraps biomaRt to get attributes
#'  from the rownames of DESeqResults objects,
#'  and adds them to the listData attribute
#'@param obj DESeqResults object to transform. Defaults to object named ddr
#'@param attrvec vector of character strings corresponding to biomaRt attributes to fetch for annotating your results object. Defaults to c('hgnc_symbol','ensembl_transcript_id')
#'@param filtvec vector of character strings you will use for the key to select biomaRt entries for fetching annotations. Defaults to c('ensembl_transcript_id')
#'@export
addAttrCol <- function(obj=ddr,
                       attrvec = c('hgnc_symbol','ensembl_transcript_id'),
                       filtvec=c('ensembl_transcript_id')){
  vals <- rownames(obj) %>% str_remove(string=.,pattern="\\..*")
  mart <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL','hsapiens_gene_ensembl')
  annotDF <- getBM(attrvec,filters=filtvec,values=vals,mart)
  annotDF <- annotDF[match(str_remove(string=rownames(obj),pattern = "\\..*"),
                annotDF$ensembl_transcript_id),]
    obj <- cbind(obj,annotDF)
  return(obj)
}




#'@title transform expression data
#'@description wrapper function for transformations implemented by DESeq2 on count data for downstream analyses like PCA
#'@details wrapper function for transformations implemented by DESeq2,
#'  either [DESeq2::rlog()] or [DESeq2::varianceStabilizingTransformation()],
#'  with options defaulting to those of 8OGseq paper usage (namely, blind=FALSE)
#'  Note that you should loop over your ddsObjList as shown in the third example.
#'@param data DESeqDataSet object to transform, Default is dds
#'@param tfm the transformation to perform, must be one of 'rlog' (for [DESeq2::rlog()]) or 'varStbl' (for [DESeq2::varianceStabilizingTransformation()]). Default is 'rlog'
#'@param blind whether the transformation is blinded to experimental design, see documentation for [DESeq2::rlog()] or [DESeq2::varianceStabilizingTransformation()] for more information. Default in this package is "FALSE"
#'@return the transformed DESeqDataSet object as class DESeqTransform
#'@export
#'@param ... any option that can be passed to [DESeq2::rlog()] or [DESeq2::varianceStabilizingTransformation()]
#'@examples
#'dds.tfm <- ddsTfm(dds,'rlog',blind=FALSE)
#'dds.tfm <- ddsTfm(myDds,'varStbl')
#'ddsTfmList=c();for (NAME in names(ddsObjList)){ddsTfmList[[NAME]]<-ddsTfm(ddsObjList[[NAME]])}
#'@seealso [modRIPseq::splitAndBuildDEobjs()] which should be run before this to build the DESeqDataSet
#'@seealso [DESeq2::rlog()] which this function wraps
#'@seealso [DESeq2::varianceStabilizingTransformation()] which this function wraps
#'
ddsTfm <- function(data=dds,tfm='rlog',blind=FALSE,...){
  if (tfm == 'rlog'){dds.tfm <- DESeq2::rlog(data,blind=FALSE,...)}
  else if (tfm == 'varStbl'){dds.tfm <- DESeq2::varianceStabilizingTransformation(data,blind=FALSE,...)}
  else error('second parameter, tfm, must be set to either "rlog" or "varStbl"')
  return(dds.tfm)
}


#'@title add contrast information to ddsTfmList
#'@description This function operates on a list of DESeqTransform object to add contrast information for intgroup option of [DESeq2::plotPCA()]
#'@details This function will use the passed ddrObjList parameter to identify
#'  the contrast used for the corresponding entry
#'  in the passed ddsTfmList parameter. This contrast which will be set as an
#'  additional list entry in the DESeqTransform object, and this matches the
#'  column name of the variable used in the design of the DESeq analysis.
#'  The primary purpose here is to generate a linked value for
#'  the intgroup option of PCA plot
#'@param ddstfmlist the ddsTfmList on which you'd like to set the contrast values
#'@param ddrobjlist the ddrObjList to use for finding contrast to assign to the ddstfmlist. Defaults to ddrObjList
#'@export
#'@return the ddsTfmList object with each list entry containing the DESeqTransform object with an added column for contrast
#'@examples
#' ddsTfmList <- setTfmContrast(ddsTfmList,ddrObjList)
#' myDdsTfmlist <- setTfmContrast(myDdsTfmlist,myDdrObjList)
#'@seealso [modRIPseq::ddsTfm()] which should be used to build the list of DESeqTransform objects this function operates on
#'@seealso [modRIPseq::runDEseqSet()] which should be used to build the list of DESeqResults objects this function uses to get contrast data
setTfmContrast <- function(ddstfmlist=ddsTfmList,ddrobjlist=ddrObjList){
  for (NAME in names(ddstfmlist)){
    ddr <- ddrobjlist[[NAME]]
    contraster <-ddr@elementMetadata$description[2] %>%
      stringr::str_remove(.,".*: ") %>%
      str_remove(.," .*")
    ddstfmlist[[NAME]]$contrast <- contraster
  }
  return(ddstfmlist)
}
