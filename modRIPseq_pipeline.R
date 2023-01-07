

#load factor table
factorTbl <- loadFactors("./inst/extdata/factorsAll.tbl")

#load files
fileList <- loadFileList("./inst/extdata/inputfilelistAll.txt")

#generate list of DESeqDataSet objects
ddsObjList <- splitAndBuildDEobjs()

#generate list of DESeqResults objects
ddrObjList <- runDEseqSet()

#generate list of DESeqTransform objects for ML techniques
ddsTfmList=c();
for (NAME in names(ddsObjList)){
  ddsTfmList[[NAME]]<-ddsTfm(ddsObjList[[NAME]])}

#add metadata for the intgroup option of plotPCA
ddsTfmList <- setTfmContrast()

#plot PCAs to assess sample spread
for (NAME in names(ddsTfmList)){plot(makePCA(ddsTfmList[[NAME]]))}

#add annotation metadata from biomaRt, fetching only once for the biggest and going from there
#note this function must feed ddr object, not others - build the others from that annotation
ddrObjList$Input_Low <-addAttrCol(ddrObjList$Input_Low)
for (NAME in names(ddrObjList)){
  if (NAME=="Input_Low") {
    next
    }
  else {
    for (colname in c("hgnc_symbol","ensembl_transcript_id","ensembl_gene_id","transcript_biotype")){
      tmp <- ddrObjList[["Input_Low"]][match(stringr::str_remove(rownames(ddrObjList[[NAME]]),pattern="\\..*"),ddrObjList$Input_Low$ensembl_transcript_id),]
      ddrObjList[[NAME]][[colname]] <-   tmp[[colname]]
    }

      }
}
#drop all nas in hgnc_symbol to drop discontinued transcripts
#for (NAME in names(ddrObjList)){ddrObjList[[NAME]] <- ddrObjList[[NAME]][!is.na(ddrObjList[[NAME]]$hgnc_symbol),]}




#get all significantly differentially expressed/oxidized genesets
#(padj<0.1) and for oxidized also impose log2>2 filter
#in a new list of results objects
#make this into its own function?
ddrSigList <- c()
for (NAME in names(ddrObjList)){
  res <- ddrObjList[[NAME]]
  if (grepl(".*In.*",NAME,ignore.case = TRUE)){
    ddrSigList[[NAME]] <- res[!is.na(res$padj)&res$padj<0.1,]
  }
  else {
    ddrSigList[[NAME]] <- res[!is.na(res$padj)&res$padj<0.1&(res$log2FoldChange>2),]
  }
}


for (i in 1:ncol(combn(x=names(ddrSigList),m=3))){makeVenn(samples=combn(x=names(ddrSigList),m=3)[,i])}




for (NAME in names(ddrObjList)){
  print(paste0("Get results summary of DESeq analysis for ",NAME))
  summary(ddrObjList[[NAME]])
  print(
    paste0(
      "alpha value for ",NAME,": ",
      metadata(ddrObjList[[NAME]])$alpha,
      ". Filter threshold for ",NAME,": ",
      metadata(ddrObjList[[NAME]])$filterThreshold
      )
    )
}
