# Creates and concatenates biological summary stats & metadata (V.0.0.1)

# August 2020 - NNRCO ZIFZ


# Import libaries
import time
import timeit
import os
import glob
import csv
import gzip
import numpy as np 
import pandas as pd
import argparse

# Import all functions
from combinator.clean_concatenate import read_compressed_files, read_compressed_filenames 
from combinator.clean_concatenate import check_headers, check_columns, concatenate_all_tables, add_data_to_existing


# Column names

'''
Contains all the names for dataset header
'''

metadata_column_names = ['data_set_identifier', 'study', 'cohort', 'study_type', 'study_year', 'study_genotyping_platform', 'number_of_participants_study', 'number_of_cases_study', 'number_of_controls_study', 'ancestry', 'original_build', 'trait_type', 'efo_category', 'efo_term', 'icd10', 'icd9', 'opcs4', 'opcs3', 'gene', 'transcript', 'protein', 'metabolite', 'source', 'trait', 'trait_abbreviation', 'file_raw', 'time1_raw', 'time2_raw', 'file_processed', 'time1_processed']
eQTL_column_names = ['data_set_identifier', 'gene', 'transcript', 'source', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'effect_allele_frequency', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'z', 'p', 'genotype_imputation_score', 'direction', 'number_of_participants', 'number_of_cases', 'number_of_controls', 'hetisq', 'hetdf', 'hetpval', 'hweq', 'original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score']
pQTL_column_names = ['data_set_identifier', 'gene', 'protein', 'source', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'effect_allele_frequency', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'p', 'z', 'genotype_imputation_score', 'direction', 'number_of_participants', 'number_of_cases', 'number_of_controls', 'hetisq', 'hetdf', 'heptval', 'hweq', 'original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score', 'name']
mQTL_column_names = ['data_set_identifier','metabolite', 'source', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'effect_allele_freqeuncy', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'z', 'p', 'genotype_imputation_score', 'direction', 'number_of_participants', 'number_of_cases', 'number_of_controls', 'hetisq', 'hetdf', 'hetpval', 'hweq', 'original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score', 'gene', 'fdr', 'name']
gwas_column_names = ['data_set_identifier','efo_term','chromosome','position','reference_allele','alternative_allele','snp','strand','effect_allele_frequency','minor_allele_frequency','effect_estimate','standard_error','p','z','genotype_imputation_score','direction','number_of_participants','number_of_cases','number_of_controls','hetisq','hetdf','hetpval','hweq','original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score']
ewas_column_names = ['data_set_identifier', 'probe_name', 'source', 'chromosome', 'position', 'reference_allele', 'alternative_allele', 'snp', 'strand', 'probe_chromosome', 'probe_position', 'type', 'effect_allele_frequency', 'minor_allele_frequency', 'effect_estimate', 'standard_error', 'z', 'p', 'genotype_imputation_score', 'direction', 'number_of_participants', 'number_of_cases', 'number_of_controls', 'hetisq', 'hetdf', 'hetpval', 'hweq', 'original_effect_allele', 'original_other_allele', 'original_strand', 'original_direction', 'original_effect_allele_frequency', 'statistics_imputation_score', 'gene', 'fdr', 'name']
allele_freq_ref_column_names = ['chromosome', 'allele', 'cohort', 'ethnicity', 'genotyping_method', 'alternative_allele_frequency', 'minor_allele']
allele_MAP_column_names = ['chromosome', 'position', 'allele', 'type']
ethnicity_codes = ['ETHNICITY', 'DESCRIPTION']
variant_annotation = ['chromosome', 'position', 'allele', 'vep_sequence_chance', 'gene', 'feature', 'feature_type', 'consequence', 'cdna_position', 'cds_position', 'protein_position', 'amino_acids', 'codons', 'existing_variation', 'distance', 'strand', 'sift', 'polyphen', 'motif_name', 'motif_position', 'high_inf_position', 'motif_score_changes']


# List of acceptable data types for '--type' argument
type_list = ['metadata','gwas', 'ewas','eQTL','mQTL', 'pQTL', 'allele_freq', 'allele_map', 'ethnicity', 'variant']

# Arguments 
parser = argparse.ArgumentParser(description='Combinator: This tool concatenates multiple summary datasets into a single comma-separated file')
parser.add_argument("--create_summary", "--cs", type=str, help= "Directory path that contains only files that will be concatenated")
parser.add_argument("--existing_summary", "--es", type=str, help= "Single existing summary .csv that will concatenate new data")
parser.add_argument("--summary_to_add", "--sta", type=str, help= "Single new summary file to be concatenated (.csv) ")
parser.add_argument("--type", "--t", type=str, help= "Select predetermined biological dataset.........."
                    'For Metadata summary statistics: metadata.............'
                    'For GWAS summary statistics: gwas.............' 
                    'For EWAS summary statistics: ewas.............'
                    'For Expression Quantitative Trait Loci summary statistics: eQTL.............'
                    'For Metabolite Quantitative Trait Loci: mQTL.............'
                    'For Protein Quatitative Trait Loci: pQTL.............'
                    'For Allele Frequency Reference Cohort: allele_freq.............'
                    'For Allele MAP: allele_map.............'
                    'For Ethnicity codes: ethnicity.............'
                    'For Variant Annotation: variant.............')
parser.add_argument("--gzip", "--gz", help= "Used when files are compressed (gzip)", default=None, nargs='?', const='gzip', required = False)
parser.add_argument("--out", "--o", type=str, help= "Output filename ")

MASTHEAD = "...\n"
MASTHEAD = "-------------------------------------------------------------------------\n"
MASTHEAD = "\n-------------------------------------------------------------------------\n"
MASTHEAD += "* Combinator\n"
MASTHEAD += "* Summary: creates and concatenates biological summary stats & metadata \n"
MASTHEAD += "* Version: v0.0.1\n"
MASTHEAD += "* Novo Nordisk Reseach Centre Oxford\n"
MASTHEAD += "-------------------------------------------------------------------------\n"
MASTHEAD += "-------------------------------------------------------------------------\n"



header = MASTHEAD
print(header)
time.sleep(5)

start = timeit.default_timer()
stop = timeit.default_timer() 


# Creates new summary '.csv' from multiple files in a directory    
        
if __name__ == '__main__':
    args = parser.parse_args()
    
  
    
    if args.type not in type_list:
        raise ValueError('\nDatatype not valid, please see documentation for more detail')
    else:
        if args.type is None:
            raise ValueError('\nThe --out flag is required')
    
    if args.create_summary:
        print('\nJOB: Creating new summary....')
        
        if args.gzip:
            list_of_files = read_compressed_files(args.create_summary, compression= True)
            filenames = read_compressed_filenames(args.create_summary, compression= True)
        else:
            list_of_files = read_compressed_files(args.create_summary, compression= None)
            filenames = read_compressed_filenames(args.create_summary, compression= None)
        
        files_key = dict(zip(filenames, list_of_files))
        list_of_files = check_headers(file_list = list_of_files, filenames=filenames, dataset = args.type, files_key = files_key)
        list_of_files = check_columns(file_list = list_of_files, filenames=filenames,dataset = args.type, files_key = files_key)
        
        if not list_of_files:
            raise ValueError('\n Invalid files.')
        concat_file = concatenate_all_tables(file_list = list_of_files, out=args.out)

        print('\n....Successful concatenation of' + " " + args.out)


    if args.existing_summary and args.summary_to_add:
        print('\nJOB: Adding new summary to existing file....')
        
        if args.gzip:
            if args.summary_to_add.endswith('.csv'):
                raise ValueError('\nFor files that are (.csv.gz) please use --gzip flag. For (.csv) files no flag required.')
            else:
                add_data_to_existing(existing_summary = args.existing_summary, summary_to_add = args.summary_to_add, dataset = args.type, compression = True, out = args.out)
            
        if not args.gzip:
            if args.summary_to_add.endswith('.csv.gz'):
                raise ValueError('\nFor files that are (.csv.gz) please use --gzip flag. For (.csv) files no flag required.')
            else:
                add_data_to_existing(existing_summary = args.existing_summary, summary_to_add = args.summary_to_add, dataset = args.type, compression = None, out = args.out)


execution_time = stop - start
print("\nProgram Executed in "+str(round(execution_time, 2)) + " " +"seconds\n") 


   







































