

#load factor table
factorTbl <- loadFactors("./inst/extdata/factorsAll.tbl")

#load files
fileList <- loadFileList("./inst/extdata/inputfilelistAll.txt")

#split and build DE object list
#this will also drop single rep rows
#this will also drop columns with univariate entries
#it calls
ddsObjList <- splitAndBuildDEobjs()
ddrObjList <- runDEseqSet()

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
