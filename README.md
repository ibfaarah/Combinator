# Combinator

A tool creates and concatenates biological summary stats & metadata (V.0.0.1)

## Background

Due to high-throughput bioloigcal/medical data being generated, there are a couple issues. One problem is that, datasets of the same type are generated globally and is of great interest to for researchers to curate and combine them to carry out more analysis. Another issue, is that output for results are sometimes split into several files due to their size and would require concatenation after generation. Concatenation of files are essential as they allwow for (a) optimised query look-ups, (e.g. for very variants/genes of interest) - or (b) for further insilico downstream analysis (e.g GWAS summary statistics, meta-analysis). 

<br/><br/>Combinator is a command line interface (CLI), written in python, breifly, it allow for QC'ing, data harmonising and concatenation of datasets - so that further analysis can be undertaken. Combinator, runs over both standard (.csv) file formats as well gzipped file formats (.csv.gz) for large datasets, but requires the appropriate flag to used (--gzip). The output is a single gzipped file (.csv.gz).

## Dataset types

Currently there are ten

## Usage