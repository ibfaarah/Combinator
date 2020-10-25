# Combinator

A tool creates and concatenates biological summary stats & metadata `V.0.0.5`

## Background

Due to high-throughput bioloigcal/medical data being generated, there are a couple issues. One issue is that, datasets of the same type are generated globally and is of great interest to for researchers to curate, qc and then combine - in order to carry out further  interrogation. Another issue, is that output for results are sometimes split into several files due to their size and would require concatenation after generation. Concatenation of files are essential as they allow for (a) optimised query look-ups, (e.g. for very variants/genes of interest) - or (b) for further insilico downstream analysis (e.g GWAS summary statistics, meta-analysis). 

Combinator is a command line interface (CLI), written in python, breifly, it allow for QC'ing, data harmonising and concatenation of datasets - so that further analysis can be undertaken. Combinator, runs over both standard (.csv) file formats as well gzipped file formats (.csv.gz) for large datasets, but requires the appropriate flag to used (`--gzip`). The output is a single gzipped file (.csv.gz).

## Modes
There are two modes which the programs runs on, *strict* (`--strict`) and *lenient* (default). The modes are based around the file column headers that a user is required to provide using (`--set_header`) command. This command takes two types of files in order to specify the final header. 1) Comma separated value file (.txt format) containing the list of column (ordered) 2) One of the files to be concatenated that contains the users preferred column header (.csv format). The column names provided using this parameter will serve as the column names in the final output file.

1. Strict mode (`--strict`)

This mode will, only concatenate files that contain the exact number of columns and names (case and space sensitive). Any files that do not meet this stringent requirements will be deleted, name of the file will be notified so that the user can correct this file. The argument `--strict` is necessary for this mode to run, no other value is required.

2. Lenient (default)

This is the default mode of the program. The user will still be required to provide the column names using the mandatory `--set_header`. In lenient mode, the program first checks that the list of column names provided using the `--set_header` parameter matches, if this is not the case then the program checks moves to the next step. Here the program checks to see whether the file to be concatenated by the user has the same number of columns as the columns provided in `--set_header` argument. If the number of columns match, the file will be concatenated despite the name of the columns not matching - however, the user will be presented with a warning to let them know that the concatenation will be based on 'number of columns' and not the 'column name' and may present some potential issue in later analysis. 


## Usage

Combinator has two main purposes which will be explained in this section. 1) Creating a summary file 2) Adding new summary file to an existing file. Each functions are similar, but are used in the different situations depending on the users need. 


**Creating a summary file**

This functions creates the initial summary file from all the existing files. All files to be concatenated will need to be stored in a specific directory, as the program will use all *files* in the stated directory to be processed, it hasn't got the capability to differentiate between files. In addition, all files are required to be either all compressed as (.csv.gz) or (.csv). For compressed files it necessary to pass the `--gzip` / `--gz` parameter. It is also required to provide an output name for the resulting file using `--out` / `--o`. 

*Example 1*

Processing several compressed GWAS summary statistics file stored in a directory.

```
python Combinator.py \
 --create_summary '/Users/Name/Directory/t2d_gwas_summary_folder' \
 --type gwas \
 --gzip \
 --strict \
 --out t2d_summary \
 ```

 The expected output is a compressed gzipped file irrespective of whether or not the initial directory contains compressed or standard format files.


 **Adding a summary file to existing summary file**

 This function will be used more frequently as more data sets are collected, each would be required to be concatenated to the master summary file. For this, two flags are required,  `--existing_summary` - this flag takes a single compressed (.csv.gz) file as input, the second flag is `--summary_to_add`, this also takes in a single file but can be either compressed (requires `--gzip` / `--gz`) or standard format. As with the previous example, this function also requires, `--out` / `--o`.

 *Example 2*

 Takes a single compressed GWAS summary statistics from *Example 1* and add a new T2D GWAS summary statistics.

 ```
 python Combinator.py \
 --existing_summary t2d_summary \
 --summary_to_add new_g \
 --type gwas \
 --out updated_t2d_summary \
 ```

## Dataset types

Currently for `v0.0.1`, the program takes 10 specific biological dataset types, which the user is expected state using (`--type`, see details below). For each data type, there is a predefined format for the header (see detail below). First, the program checks to see  whether the dataset header matches the expected datatype `--type` specified, if there is a mismatch the dataset will be reomved, the user will then be notified and will be required to rectify in line with the correct header (specified below) before trying again. Furthermore, the program also qc's column entries to ensure incorrect data type (e.g. integer, string etc) and outliers are removed. If these are detected, the dataframe is removed and error message, specifying the file in question and what column contains the issue.

In order for Combinator to carry out quality quality checks for columns, it is strongly recommended that you follow the column index rule for each dataset provided below. This allows for qc to be carried irrespective of the column name. 

Currently , Combinator can take in 10 specified datasets, they are the following:



1. **Study-related meta data**

When specifying this datatype on the CLI, use: `metadata`. Example column names in header: 

*'data_set_identifier', 'study', 'cohort', 'study_type', 'study_year', 'study_genotyping_platform', 'number_of_participants_study', 'number_of_cases_study', 'number_of_controls_study', 'ancestry', 'original_build', 'trait_type', 'efo_category', 'efo_term', 'icd10', 'icd9', 'opcs4', 'opcs3', 'gene', 'transcript', 'protein', 'metabolite', 'source', 'trait', 'trait_abbreviation', 'file_raw', 'time1_raw', 'time2_raw', 'file_processed', 'time1_processed'*

2. **GWAS summary statistics**

When specifying this datatype on the CLI, use: `gwas`. Example column names in header: 

*'data_set_identifier','efo_term','chromosome','position','reference_allele','alternative_allele','snp','strand','effect_allele_frequency','minor_allele_frequency','effect_estimate','standard_error','p','z','genotype_imputation_score','direction','number_of_participants','number_of_cases','number_of_controls','hetisq','hetdf','hetpval','hweq','original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score'*

#### GWAS column indexes

| Column Name              | Column Index |
|--------------------------|:------------:|
| Chromosome               | 2            |
| Position                 | 3            |
| Effect Allele Frequency  | 4            |
| Minor Allele Frequency   | 5            |
| Standard Error           | 11           |
| P-value                  | 12           |
| Strand                   | 7            |
| Direction                | 15           |
| Original Strand          | 25           |
| Original Direction       | 26           |


3. **eQTL summary statistics**

When specifying this datatype on the CLI, use:  `eQTL`. Example column names in header: 


*'data_set_identifier', 'gene', 'transcript', 'source', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'effect_allele_frequency', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'z', 'p', 'genotype_imputation_score', 'direction', 'number_of_participants', 'number_of_cases', 'number_of_controls', 'hetisq', 'hetdf', 'hetpval', 'hweq', 'original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score'*

#### eQTL column indexes

| Column Name              | Column Index |
|--------------------------|:------------:|
| Chromosome               | 4            |
| Position                 | 5            |
| Effect Allele Frequency  | 10           |
| Minor Allele Frequency   | 11           |
| Standard Error           | 13           |
| P-value                  | 15           |
| Strand                   | 9            |
| Direction                | 17           |
| Original Strand          | 27           |
| Original Direction       | 28           |

4. **pQTL summary statistics**

When specifying this datatype on the CLI, use: `pQTL`. Example column names in header: 

*'data_set_identifier', 'gene', 'protein', 'source', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'effect_allele_frequency', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'p', 'z', 'genotype_imputation_score', 'direction', 'number_of_participants', 'number_of_cases', 'number_of_controls', 'hetisq', 'hetdf', 'heptval', 'hweq', 'original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score', 'name'*

#### pQTL column indexes

| Column Name              | Column Index |
|--------------------------|:------------:|
| Chromosome               | 4            |
| Position                 | 5            |
| Effect Allele Frequency  | 10           |
| Minor Allele Frequency   | 11           |
| Standard Error           | 13           |
| P-value                  | 14           |
| Strand                   | 9            |
| Direction                | 17           |
| Original Strand          | 27           |
| Original Direction       | 28           |

5. **mQTL summary statistics**

When specifying this datatype on the CLI, use: `mQTL`. Example column names in header: 

*'data_set_identifier', 'metabolite', 'source', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'effect_allele_frequency', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'z', 'p', 'genotype_imputation_score', 'direction', 'number_of_participants', 'number_of_cases', 'number_of_controls', 'hetisq', 'hetdf', 'hetpval', 'hweq', 'original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score', 'gene', 'fdr', 'name'*

#### mQTL column indexes

| Column Name              | Column Index |
|--------------------------|:------------:|
| Chromosome               | 3            |
| Position                 | 4            |
| Effect Allele Frequency  | 9            |
| Minor Allele Frequency   | 10           |
| Standard Error           | 12           |
| P-value                  | 14           |
| Strand                   | 8            |
| Direction                | 16           |
| Original Strand          | 26           |
| Original Direction       | 27           |

6. **EWAS summary statistics**

When specifying this datatype on the CLI, use: `ewas`. Example column names in header: *'data_set_identifier', 'probe_name', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'probe_chromosome', 'probe_position', 'type', 'effect_allele_frequency', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'z', 'p', 'genotype_imputation_score', 'direction', 'number'*

7. **Allele frequency reference cohorts**

When specifying this datatype on the CLI, use: `allele_freq`. Example column names in header: 

*'chromosome', 'allele', 'cohort', 'ethnicity', 'genotyping_method', 'alternative_allele_frequency', 'minor_allele'*

8. **Allele MAP**

When specifying this datatype on the CLI, use: `allele_map`. Example column names in header: 

*'chromosome', 'position', 'allele', 'type'*

9. **Ethinicity codes**

When specifying this datatype on the CLI, use: `ethnicity`. Example column names in header: 

*'ETHNICITY', 'DESCRIPTION'*

10. **Variant annotation**

When specifying this datatype on the CLI, use: `variant`. Example column names in header: 

*'chromosome', 'position', 'allele', 'vep_sequence_chance', 'gene', 'feature', 'feature_type', 'consequence', 'cdna_position', 'cds_position', 'protein_position', 'amino_acids', 'codons', 'existing_variation', 'distance', 'strand', 'sift', 'polyphen', 'motif_name', 'motif_position', 'high_inf_position', 'motif_score_changes'*




## License
This project is licensed by Novo Nordisk Research Centre Oxford (NNRCO)

## Authors
Ibrahim Bashe Farah (zifz@novonordisk.com, NNRCO Computational Biology)
