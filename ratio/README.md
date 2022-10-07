# Fragmentomics
This is a tool used to calculate genome-wide fragmentation profiles of plasma samples using shallow whole genome sequencing. 

The scripts and concepts are based off of: 
Cristiano, S., Leal, A., Phallen, J. et al. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature 570, 385–389 (2019). https://doi.org/10.1038/s41586-019-1272-6

## Description
This tool calculates the ratio of short (90-150bp) to long (151-220bp) DNA fragments across the genome in 100kb bins. Areas excluded from the fragment count include centromeres, telomeres, genome gaps, varaible number tandem repeats, and repetitive regions described in the ENCODE blacklist.

## Usage
```
Rscript runFrag.R\
 --id <sample id>\
 --bam <path/to/bamfile>\
 --filters <path/to/filters>\
 --gaps <path/to/gaps>\
 --VNTRs <path/to/VNTRs>\
 --tiles <path/to/tiled_genome>\
 --healthy <path/to/healthy_median>\
 --outdir <path/to/output_directory>\
 --libdir <path/to/fragmentomics>
 ```
This tool is memory intensive for large bam files.
 
Filters, gaps, VNTRs, and tiles files for hg38 can be found in the extdata folder.
Scripts to generate your own filters, gaps, VNTRs can be found in the filtered_regions folder.
 
Healthy median was generated using healthy control plasmas from the pughlab.
Scripts to generate your own health median can be found in the healthy_median folder.
