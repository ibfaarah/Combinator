# Combinator

A tool creates and concatenates biological summary stats & metadata (V.0.0.1)

## Background

Due to high-throughput bioloigcal/medical data being generated, there are a couple issues. One problem is that, datasets of the same type are generated globally and is of great interest to for researchers to curate and combine them to carry out more analysis. Another issue, is that output for results are sometimes split into several files due to their size and would require concatenation after generation. Concatenation of files are essential as they allow for (a) optimised query look-ups, (e.g. for very variants/genes of interest) - or (b) for further insilico downstream analysis (e.g GWAS summary statistics, meta-analysis). 

Combinator is a command line interface (CLI), written in python, breifly, it allow for QC'ing, data harmonising and concatenation of datasets - so that further analysis can be undertaken. Combinator, runs over both standard (.csv) file formats as well gzipped file formats (.csv.gz) for large datasets, but requires the appropriate flag to used ('--gzip'). The output is a single gzipped file (.csv.gz).

## Dataset types

Currently, the program only takes 10 specific biological dataset types, which the user is expected state using ('--type', see details below). For each data type, there is a predefined format for the header (see detail below). First, the program checks to see  whether the dataset header matches the expected datatype '--type' specified, if there is a mismatch the dataset will be reomved, the user will then be notified and will be required to rectify in line with the correct header (specified below) before trying again. Furthermore, the program also qc's column entries to ensure incorrect data type (e.g. integer, string etc) and outliers are removed. If these are detected, the dataframe is removed and error message, specifying the file in question and what column contains the issue.

Currently (v0.0.1), Combinator can take in 10 specified datasets, they are the following:

1. **Study-related meta data**
2. **GWAS summary statistics**
3. eQTL summary statistics
4. pQTL summary statistics
5. mQTL summary statistics
6. EWAS summary statistics
7. Allele frequency reference cohorts
8. Allele MAP
9. Ethinicity codes
10. Variant annotation


## Usage