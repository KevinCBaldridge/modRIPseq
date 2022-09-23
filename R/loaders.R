# library(tidyverse)
# library(tximport)
# library(DESeq2)

#essentially global variables I will need:
#factor.tbl - path to file containing tab-separated table, with column names as first row, where the first column should be the filenames, second third and fourth columns are factors with appropriate first-row labels
#inputfilelist.txt - path to file containing newline-separated list of files that are the output files from RSEM at the transcript level. These filenames should match the filenames in the factor.tbl first column. (will add error checking to confirm this)
#done - note actually that the default filename output from rsem contains "isoform" or "gene", so we can have that be set from filenames automatically
#note that the functions in this file, none should use the biocparallel stuff...

#To do:
#DONE - add fucntionality for transcrtipt vs whole gene level
#add functionality for other count tables besides RSEM
#rework it all to be in line with my requirements for factors.tbl
#add processing of path so dont have to change wd partway through
#finish setddsnames in line with new requirements
#build initial formatting/error checking into all functions
#write the script for running all these functions, vignette?
#write all documentation Roxygen stuff




#load treatment factors (Ab and Exposure) file into table
#'@title Load factors table for count table files
#'@description This function will load the factors table file for annotating samples, which will be used in other package functions to set up analysis of count data from gene expression
#'@details The input file to this function should be a tab-separated table with:
#'  first row containing titles for each column
#'  first column should be the filename (i.e. filename or folder/filename), coresponding to the files listed on each line of the inputfile.txt file fed into the loadFileList function.
#'  second column should be the antibody/pulldown treatment factor (i.e. Yes vs No)
#'  third column should be the exposure condition (i.e. CAC/clean air/no vs OAM/exposure/yes)
#'  fourth column should be the exposure level (i.e. high vs low), if you have multiple levels of exposure in your data
#'@param filepath Path to the file containing your table of factors. Defaults to example file "./data/factors.tbl"
#'@export
#'@return Tibble object containing the factors table corresponding to filenames
#'@examples
#'loadFactors('./data/factors.tbl')
#'loadFactors('mydata/myfactorsTable.txt')
#'@seealso [modRIPseq::splitAndBuildDEobjs()] which follows this function in the modRIPseq pipeline
loadFactors <- function(filepath='./data/factors.tbl'){
  dir <- dirname(filepath)
  factorTbl <- read_table(filepath)
  for (i in seq_along(factorTbl)) factorTbl[[i]] <- as_factor(factorTbl[[i]])
  factorTbl$file <- file.path(dir,factorTbl$file)
  return(factorTbl)
}

#'@title Load a file containing list of count table files
#'@description This function will load into memory a list of file paths for loading the actual count table files
#'@details The file input to this function should be a simple list of file paths to the output files from RSEM,
#'  organized with each file path name on a newline, without any header row in the file.
#'  the file paths in each newline should correspond to the first column in the factors table file passed as a parameter to the loadFactors function
#'@param filepath Path to the file containing list of filenames for RSEM output files. Defaults to "data/inputfilelist.txt"
#'@export
#'@return Vector of filenames pointing to the files containing RSEM count data
#'@examples
#'loadFileList("./data/inputfilelist.txt")
#'loadFileList("./myData/myFileList.txt")
#'@seealso [[modRIPseq::loadFactors()]] which is the function that loads factor tables with first column corresponding to the fileList returned by this function
loadFileList <- function(filepath='./data/inputfilelist.txt'){
  dir <- dirname(filepath)
  fileList <- scan(file=filepath,sep="\n",what="")
  fileList <- basename(fileList)
  fileList <- file.path(dir,fileList)
  return(fileList)
}



#'@title Load RSEM output files
#'@description This function will load the output count files from RSEM
#'@details Before running this, you should first execute loadFileList() to read the list of filenames.
#'  In your inputfilelist.txt, ensure that your input files are all of the same output type, - or gene-level summaries.
#'  This function will check that the filenames are consistent to ensure you have consistent analyses
#'  Additionally, this function can call a secondary function for row name checking to ensure that all individual components
#'  (counts, abundances, and lengths) are consistently named correctly within each file and across files
#'@param filelist Vector of filenames with each element containing a string that indicates the path to an rsem output count table file. Defaults to fileList.
#'@param checkRowNaming (OPTIONAL) Boolean flag to check all row names within files and across files to ensure row name consistency in counts, lengths, and abundances. Defaults to FALSE.
#'@export
#'@return List of matrices corresponding to count, abundance, and lengths for all the count tables indicated in your input file list (same as return from tximport)
#'@examples
#'loadRSEMs(fileLlist)
#'loadRSEMs(myFileList,TRUE)
#'@seealso [tximport::tximport()] which this function wraps
#'@seealso [modRIPseq::loadFileList()] which is the precursor to this function in the modRIPseq pipeline
#note that the invisible doesnt work right here when tximport is loaded
#note also that as it stands currently, need to navigate to subdirectory
#  - add a test if there's a subfolder in the parameter fed to inputfilelist.txt, check if it matches the filenames or not, and then navigate if necessary to that directory before loading
loadRSEMs <- function(filelist=fileList,checkRowNaming=FALSE){
  iso <- unique(grepl('isoform',filelist))
  gene <- unique(grepl('gene',filelist))
  if(!length(iso)==1){stop('inputfilelist.txt appears to have a mixture of isoform- and gene-level summaries, please ensure consistency')}
  if(length(iso)==1 && iso){
    print("automatically identified as isoform-level summary from rsem-calculate-expression")
    txi.rsem <- invisible(tximport(files=filelist,type='rsem',txOut=TRUE))
  }
  if(length(gene)==1 && gene){
    print("automatically identified as gene-level summary from rsem-calculate-expression")
    txi.rsem <- invisible(tximport(files=filelist,type='rsem',txOut=FALSE))
  }
  if(checkRowNaming==TRUE){
    print("Checking row naming...")
    txi.rsem <- checkRowNaming(txi.rsem,filelist)
  }
  return(txi.rsem)
}


#'@title Row naming checks for tximport returned object
#'@description This function will check row naming for the txi.rsem object loaded by loadRSEMs() function
#'@details The function checks the first column of the matrices in the list returned by loadRSEMs (or another tximport return object) to ensure consistency in row names among files and the various matrices (counts, abundance, lengths)
#'  Additionally, if the row names don't check out properly, this will rename them accordingly, but you should manually check that things are correct/consistent among the various matrices (txi.rsem$counts, txi.rsem$lengths, txi.rsem$abundance) in the list
#'  This function generally should not be called directly, use the optional parameter checkRowNaming=TRUE into loadRSEMs() to call indirectly
#'@param txobj the txi.rsem object loaded by loadRSEMs(). No default
#'@param filelist the filelist passed into loadFileList() for checking rownames against. No default
#'@export
#'@return List of matrices corresponding to count, abundance, and lengths for all the count tables indicated in your input file list (same as return from tximport), with corrected row names if they did not pass the checks
#'@examples
#'checkRowNaming(txi.rsem,fileList)
#'@seealso [modRIPseq::loadRSEMs()] which is the function that calls this validation function
checkRowNaming <- function(txobj,filelist){
#  print(filelist)
#  print(head(txobj))
#should I have this loop over all files rather than just the first, checking for consistency among files as well as the various components of the txi object?
  conversiontable <- read.table(sep="\t",header=TRUE,file=filelist[1],stringsAsFactors=FALSE)
#  print(head(conversiontable[,1]))
#  print(head(txobj$counts))
  tmp <- unique(conversiontable[,1]==rownames(txobj$counts))
#  print(tmp)
  if(length(tmp)==1){
    if(tmp!=TRUE) {
      warning("Rownames for txi.rsem$counts are not the same as first column of your tximport file, auto-assigning now... but you should manually confirm after!")
      rownames(txobj$counts) <- conversiontable[,1]
    }
  }
  tmp <- unique(conversiontable[,1]==rownames(txobj$lengths))
  if(length(tmp)==1){
    if(tmp!=TRUE) {
      warning("Rownames for txi.rsem$lengths are not the same as first column of your tximport file, auto-assigning now... but you should manually confirm after!")
      rownames(txobj$lengths) <- conversiontable[,1]
    }
  }
  tmp <- unique(conversiontable[,1]==rownames(txobj$abundances))
  if(length(tmp)==1){
    if(tmp!=TRUE) {
      warning("Rownames for txi.rsem$abundances are not the same as first column of your tximport file, auto-assigning now... but you should manually confirm after!")
      rownames(txobj$abundances) <- conversiontable[,1]
    }
  }
  print("Row names check out")
  return(txobj)
}

#'@title Reorders factor tibble to match file list
#'@description This function will ensure the correct order of factors in the tibble according to the inputfilelist to properly match for DESeq analysis
#'@details This function is only called generally from within the [modRIPseq::splitAndBuildDEobjs()] function
#'@param factortibble The tibble of factors from loadFactors() that has annotation corresponding to RSEM output file names. Defaults to factorTbl
#'@param filecol The column header in the tibble of factors that corresponds to the file names. Should be acceptable quoted or unquoted as passed parameter. Defaults to file
#'@param filelist The filelist vector returned from loadFileList() for reordering the tibble of factors to ensure proper annotation. Defaults to fileList
#'@export
#'@return The tibble of factors from loadFactors() that has been rearranged to match order of the file list vector from loadFileList()
#'@examples
#'reorderFactorsByFile(factorTbl,file,fileList)
#'@seealso [modRIPseq::loadFileList()]
#'@seealso [modRIPseq::loadFactors()]
reorderFactorsByFile <- function(factortibble=factorTbl,filecol=file,filelist=fileList){
  factortibble <- factortibble %>% arrange(factor({{filecol}},levels=filelist))
  return(factortibble)
}

#
# #'@title Set names for DESeq object list
# #'@description This function will dynamically set names for each DESeqDataSet based on factor tibble
# #'@details This function will dynamically set names for the DESeqDataSet objects within reason, but it requires some loose structure around column naming in factor tibble
# #'  Specifically:
# #'    One column must be named something in the following list (case insensitive) corresponding to the replicate/trial number:
# #'      replicate
# #'      experiment (entries in column must be integers for this nomenclature to prevent confusion with exposure condition possibilities)
# #'      trial
# #'      sample
# #'    One column must be named something in the following list (case insensitive) corresponding to the antibody treatment/pulldown/IP condition:
# #'      AntibodyTreatment
# #'      pulldown
# #'      immunoprecipitation
# #'      IP
# #'      Ab
# #'      Antibody
# #'
# #'@param factortibble The tibble of factors that will be used to set DESeqDataSet object names in
# #'@return The list of names which will be used for denoting individual formulaList and ddsObjList entries
# #'@examples
# #'@seealso
# setDDSnames <- function(factortibble = factorTbl){
#   namelist <- c()
#   namelist[["Input"]] <- "Input"
#   namelist[["IP_control"]] <- "IP_control"
#   uniqentries <- factortibble %>%
#     select(exposureCondition) %>%
#     unique() %>%
#     as_vector() %>%
#     as.character()
#   uniqentries <- uniqentries[! uniqentries %in% c("control")]
#   for (i in 1:length(uniqentries)){
#     namelist[[paste0("IP_",uniqentries[i])]] <- paste0("IP_",uniqentries[i])
#     }
#   return(namelist)
# }


#'@title Set names for DESeq object list
#'@description This function will dynamically set names for each DESeqDataSet based on factor tibble
#'@details This function will dynamically set names for the DESeqDataSet objects according to the requirements of the input files.
#'@param factortibble The tibble of factors that will be used to set DESeqDataSet object names in
#'@return The list of names which will be used for denoting individual formulaList and ddsObjList entries
#'@examples
#'setDDSnamesDose(factorTbl)
#'setDDSnames
#'@seealso [[modRIPseq::loadFactors()]] which is a prerequisite for this function
#'@seealso [[modRIPseq::splitAndBuildDEobjs()]] which calls this function
setDDSnamesDose <- function(factortibble = factorTbl){
  namelist <- c()
  #namelist[["Input"]] <- "Input"
  #namelist[["IP_control"]] <- "IP_control"
  uniqentries <- factortibble %>%
    filter(abTreatment=="none") %>%
    distinct() %>%
    select(exposureLevel) %>%
    as_vector %>%
    as.character()
  for (i in 1:length(uniqentries)){
    namelist[[paste0("Input_",uniqentries[i])]] <- paste0("Input_",uniqentries[i])
  }
  uniqentries <- factortibble %>%
    select(abTreatment,exposureCondition,exposureLevel) %>%
    filter(abTreatment != "none") %>%
    distinct() %>%
    mutate(newcol=paste0(abTreatment,"_",exposureCondition,"_",exposureLevel)) %>%
    select(newcol) %>%
    as_vector %>%
    as.character()
  #uniqentries <- uniqentries[!grepl("control",uniqentries)]
  for (i in 1:length(uniqentries)){
    namelist[[paste0(uniqentries[i])]] <- paste0(uniqentries[i])
  }
  return(namelist)
}



#need to specify that factorTbl must have something indicating replicate number...
# dropSingleReps <- function(factortibble = factorTbl){
#   ncolsval <- 2^(ncol(factortibble)-2)
#   factortibble <- factortibble %>%
#     group_by(replicate) %>%
#     filter(!n()<ncolsval) %>%
#     ungroup()
#   return(factortibble)
# }
#
# dropSingleLevelFactors <- function(factortibble = factorTbl){
#   factortibble <- factortibble %>% select(!where(~length(unique(.))==1))
#   return(factortibble)
# }



#need to add a test and warning for if a factor only has one member of a level, and drop it, i.e. the R5 sample

#'@title Build DESeq2 objects for subsequent differential expression  & oxidation analysis
#'@description This function creates DESeq2 dataset objects for differential expression & oxidation analysis from the various file loading functions in the modRIPseq pipeline.
#'@details This function will use the other objects (factorTbl,fileList) created in the modRIPseq pipeline to create a set of DESeq dataset objects for subsequent analysis of exposure-induced differential expression and oxidation of transcripts.
#'@param factortibble The reordered tibble of factors. Defaults to factorTbl
#'@param filelist The vector of filenames for rsem output files, same as the vector returned from loadFileList()
#'@export
#'@return List of DESeqDataSet objects, named with "input","IPc","IPx"
#'@examples
#'splitAndBuildDEobjs(factorTbl,fileList)
#'splitAndBuildDEobjs(myFactorTbl,myFileList)
#'@seealso [DESeq2::DESeqDataSetFromTximport()] which this function wraps
#'@seealso [modRIPseq::loadFileList()] and [modRIPseq::loadFactors()] which are the functions earlier in the modRIPseq pipeline that prepare inputs for this function
splitAndBuildDEobjs <- function(factortibble=factorTbl,filelist=fileList){
  #this function will assume the factor.tbl column naming, see requirements
  #need to move these singlelevel factors and single reps to later stage of processing...
  factortibbleSet <- c()
#  if("exposureLevel" %in% colnames(factortibble)){
    #This block works now - need to move the singlerepremoval to a later stage though
    nameList <- setDDSnamesDose(factortibble)
    for (n in names(nameList)){
      parsedname <- str_split(n,"_",simplify=TRUE)
        if (parsedname[1]=="Input") {
          factortibbleSet[[n]] <- factortibble %>%
            filter(abTreatment == "none",exposureLevel==parsedname[2])
          #factortibbleSet[[n]] <- dropSingleLevelFactors(factortibbleSet[[n]])
          #factortibbleSet[[n]] <- dropSingleReps(factortibbleSet[[n]])
          factortibbleSet[[n]] <- factortibbleSet[[n]] %>%
            group_by(replicate) %>%
            filter(n()!=1) %>%
            ungroup()
        } else {
          factortibbleSet[[n]] <- factortibble %>%
            filter(exposureCondition==parsedname[2],
                   exposureLevel==parsedname[3])
          #factortibbleSet[[n]] <- dropSingleLevelFactors(factortibbleSet[[n]])
          #factortibbleSet[[n]] <- dropSingleReps(factortibbleSet[[n]])
          factortibbleSet[[n]] <- factortibbleSet[[n]] %>%
            group_by(replicate) %>%
            filter(n()!=1) %>%
            ungroup()
              }
          }

    # } else { #need to double check this block works right...
    # nameList <- setDDSnames(factortibble)
    # for (n in names(nameList)){
    #     parsedname <- str_split(n,"_",simplify=TRUE)
    #     if("IP"==parsedname[1]){
    #       filtervec <-parsedname[2]
    #       factortibbleSet[[n]] <- factortibble %>%
    #         filter(exposureCondition %in% filtervec)
    #       #factortibbleSet[[n]] <- dropSingleLevelFactors(factortibbleSet[[n]])
    #       #factortibbleSet[[n]] <- dropSingleReps(factortibbleSet[[n]])
    #       #print(factortibbleSet[[n]])
    #       factortibbleSet[[n]] <- factortibbleSet[[n]] %>%
    #         group_by(replicate) %>%
    #         filter(n()!=1) %>%
    #         ungroup()
    #     } else if ("Input" == parsedname[1]){
    #       factortibbleSet[[n]] <- factortibble %>%
    #         filter(abTreatment == "none")
    #       #factortibbleSet[[n]] <- dropSingleLevelFactors(factortibbleSet[[n]])
    #       #factortibbleSet[[n]] <- dropSingleReps(factortibbleSet[[n]])
    #       factortibbleSet[[n]] <- factortibbleSet[[n]] %>%
    #         group_by(replicate) %>%
    #         filter(n()!=1) %>%
    #         ungroup()
    #     } else {simpleError("Not matching correctly for naming stuff")}
    #   }
    # }
#can I set these dynamically based on column names? might be too much...
#well below doesnt work right yet, need to step through it...
designs <- c()
filelistlist <- c()
files<-c()
ddslist <- c()
txi.rsem <- c()
for (n in names(factortibbleSet)){
  if(str_detect(n,"In")){
    designs[[n]]<-~exposureCondition
  } else {
    designs[[n]]<-~replicate+abTreatment
    }
  files[[n]] <-filelist[filelist %in% factortibbleSet[[n]]$file]
  #print(files[[n]])
  #filelistlist[[n]] <-eval(sym(paste0("files",n)))
  #print(filelistlist[[n]])
  factortibbleSet[[n]] <- reorderFactorsByFile(factortibbleSet[[n]],file,filelistlist[[n]])
  txi.rsem[[n]] <- loadRSEMs(files[[n]])

  ddslist[[n]] <-   DESeqDataSetFromTximport(txi.rsem[[n]],
                                  colData = factortibbleSet[[n]],
                                  design=designs[[n]])
  }
    return(ddslist)
}






#note to self, I think I'd like to make this figure things out on its own the ordering, rather than feeding numbered rows into the DEseqdataset constructors, so to that end I'll have the factors.tbl first column be the filename, then the remainder be the factors themselves.
#will probably need rigid column naming for the factors though...


#note to self, can't find factors.tbl that is read in my real dissertation work. i have output logs, result files, etc, but can't find the script and inputs for some reason...
