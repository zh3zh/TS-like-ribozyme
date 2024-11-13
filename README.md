# TS-like ribozyme

# Introduction
This is modified version of CODA for RtcB ligation based deep mutational scanning. It is used to generate the base-pairing map from deep mutational sequencing data of RNA ribozyme, with the method introduced in the paper. We also provided the processed .ra files (L1-fl.ra, L1-mini.ra, OR-fl.ra) from different DMS experiments. 

# Requirement

Before running this pipline, make sure these programs were correctly installed.

  PEAR (https://cme.h-its.org/exelixis/web/software/pear/)
  
  python 2.7.15
  
  samtools 0.0.18
  
  java 1.7.0_85


# Usage

  `bash run.sh $DNAFILE1 $DNAFILE2 $RNAFILE1 $RNAFILE2 $FASTAFILE $OUTPUTPATH`

# Arguments:
  
  ```
  $DNAFILE1: first read DNA sequencing file, in gziped fastq format (fq.gz)
  
  $DNAFILE2: second read DNA sequencing file, in gziped fastq format (fq.gz)
  
  $RNAFILE1: first read RNA sequencing file, in gziped fastq format (fq.gz)
  
  $RNAFILE2: second read RNA sequencing file, in gziped fastq format (fq.gz)
  
  $FASTAFILE: base sequence file in fasta format (ATGC sequence)
  
  $OUTPUTPATH: all output files will be written here, empty file folder is recommended
```

# Outputs:

  var.count: cleaved read number of each variant, total number of each variant in the input DNA-seq library
  
  var.ra: organized relative activity of each variant
  
  var.pos.ra: organized relative activity of all mutants of each position pair
  
  pred.mtx: ps score matrix
  
  pred.ss: 100 predicted secondary structure in the bracket format with a consensus prediction
  
# Reference:

  Zhang et al. eLife 2023;13:RP90254. DOI: https://doi.org/10.7554/eLife.90254
