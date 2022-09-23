# modRIPseq
R package for analyzing immunoprecipitation-RNAseq for identifying transcript-specific RNA modifications (i.e. 8OG) at the whole transcript level

This package will work to provide an easily usable to replicate the 8OGseq analysis performed in: Gonzalez-Rivera, Juan C., et al. "Post-transcriptional air pollution oxidation to the cholesterol biosynthesis pathway promotes pulmonary stress phenotypes." Communications biology 3.1 (2020): 1-16., based upon the original code as published in: Baldridge, Kevin Charles. Insights into the functions of RNA post-transcriptional modifications gained through studies in cellular stress. Diss. 2017.

This is a work in progress! Look forward to upcoming changes. 


(planned sections)

Limitations of this package
  This package is specifically designed to wrap the DESeq2 R package for analyzing RNAseq data using RNA-modification-specific antibody enrichment of unfragmented, rRNA-depleted RNA. Notably, the original application was for enrichment of oxidized transcripts in cell models of air pollution exposure, i.e. 8OGseq, but it SHOULD be fine for use with similar studies using different antibodies specific for other RNA modifications, provided that your input RNA was unfragmented+rRNA depleted. There are other available software packages for analyzing other RNA preparations, i.e. fragmented transcriptomes.

Usage
