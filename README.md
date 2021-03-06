# Origami-HiChIP

## Overview

This is a beta version last modified on October, 2017. 

Origami-HiChIP is a pipeline for processing and calling high-confidence chromatin loops generated by HiChIP. The pipeline implements a semi-Bayesian two-compenent mixture model to learn which contacts within the data set appear to represent apparent structured chromatin loops in a cell population while correcting for important biases in HiChIP data (notably, linear genomic distance between the ends of a potential chromatin loop).

## Installation

To install origami, first run the configure script to test for the presence of mdost of the dependencies origami needs. After configure completes successfully, run make to compile the C++ code. Then  make install will copy all the necessary files to the installation directory.

(Note: at the present time, the configure script does not check if Bamtools is installed correctly, please make sure it is! The Makfile may need to be adjusted to adapt to where Bamtools API files are installed on your system.)

### Depedencies
* cutadapt (tested with version 1.8)
* samtools (tested with version 1.2)
* bamtools (tested with version 2.2.3)
* bowtie1
* MACS1 and MACS2
* R (also requires the R libraries GenomicsRanges, matrixStats)
* g++
* bash
* make
* autoconf


## Running the pipeline

The pipeline has three steps that need to be run:

* `origami-alignment` -- aligns the reads
* `origami-analysis` -- estimates the confidence in each interaction
* `origami-conversion` -- converts the origami output to popular genomic file formats for visualization


  
  ***Note: scripts were modified from work with Daniel Day during his postdoc in Rick Young's lab at Whitehead Institute. 
