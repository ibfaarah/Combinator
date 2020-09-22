# Import libaries
import os
import glob
import re
import csv
import numpy as np 
import pandas as pd


# path = '/Users/ibrahim/Documents/Novo Nordisk/Combinator/data'
# basename = os.path.basename(path)

def read_compressed_files(path):
    '''
    Takes 'path' to directory containing all the datasets.
    Whil

    Tip: 
    (1) Add 'error_bad_lines=False' to pandas method.
    (2) When implementing this CycleCloud change '/*.csv' to '/*.gz'
    '''
    file_list = []

    directory = glob.glob(path + "/*.csv")

    for file in directory:
        df = pd.read_csv(file, engine='python', compression = 'infer', header = 0, sep=',')
        file_list.append(df)
    
    return file_list


# list_of_files = read_compressed_files(path)



def read_compressed_filenames(path):
    '''
    Takes 'path' to directory containing all the datasets.
    Whil

    Tip: 
    (1) Add 'error_bad_lines=False' to pandas method.
    (2) When implementing this CycleCloud change '/*.csv' to '/*.gz'
    '''

    filenames = glob.glob(path + "/*.csv")
    
    return filenames


# filenames = read_compressed_filenames(path)
# print(filenames)

# files_key = dict(zip(filenames, list_of_files))
# print(files_key)
# print(files_key['data/chr5_baso.csv'])
# test_key = list(files_key.keys())[0]
# print(test_key)


def check_columns(file_list, filenames, dataset, files_key):
    '''
    Checks all that data type for all columns are as expected:
    Use: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.dtypes.html
    '''
    
    logging_df = []
    error = []

    chr_list = ['1', '2', '3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y', 'MT']
    directions = ['+','-']
    # stats_imp_score = ['NA']
    
    
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        dataframe_to_remove= []
        
        # This is for 'GWAS summary statistics data
        if dataset == 'gwas' or 'eQTL' or 'pQTL' or 'mQTL' or 'ewas':
            # Checks that all entries in 'chromosome' column are correct, 1-22 & X, Y, MT
            if ~dataframes.chromosome.astype(str).isin(chr_list).all() == True:
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Chromosome column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass
            
            # Checks that all entries in 'position' column have correct data type, only integers (int64)
            if dataframes.position.dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Position column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass
               
            # Checks that all entries in 'effect allele frequency' column has values between (0 & 1)
            if dataframes.effect_allele_frequency.between(0,1).all() != True:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Effect allele frequency column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass

            # Checks that all entries in 'minor allele frequency' column has values between (0 & 0.5)    
            if dataframes.minor_allele_frequency.between(0,0.5).all() != True:
                dataframe_to_remove.append(index)
                
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Minor Allele Frequency column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                
                print(dataframe_to_remove)
            else:
                pass

            # Checks that all entries in 'Standard error' column has values greater than (0)    
            if dataframes.standard_error.any() <= 0:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Standard error column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                print('Standard error column, has incorrect entry')
    
            else:
                pass
            
            # Checks that all entries in 'p-value' column has values between (0 & 1)    
            if dataframes.p.between(0,1).any() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "P-value column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass
            
            # Checks that all entries in 'strand' column are either ('+' or '-' )
            if dataframes.strand.astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Strand column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass

            # Checks that all entries in 'direction' column are either ('+' or '-' )
            if dataframes.direction.astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Direction column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass
            
            # Checks that all entries in 'original strand' column are either ('+' or '-' )
            if dataframes.original_strand.astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Original strand column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass
            
            # Checks that all entries in 'original direction' column are either ('+' or '-' )
            if dataframes.original_direction.astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Original direction column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass
            
           
        elif dataset == 'allele_map' or 'variant':
            
            if ~dataframes.chromosome.astype(str).isin(chr_list).all() == True:
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Chromosome column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            
            else:
                pass

            # Checks that all entries in 'position' column have correct data type, only integers (int64)
            if dataframes.position.dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "Position column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass

        elif  dataset == 'variant':
            if dataframes.cdna_position.dtype != np.int64:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "cDNA position column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass

            if dataframes.cds_position.dtype != np.int64:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('ERROR found in "cds position column\"' + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
            else:
                pass
                   
            pass       
    return file_list

# check_columns = check_columns(file_list = list_of_files, filenames = filenames, dataset='gwas', files_key=files_key)
# print(check_columns)



def concatenate_all_tables(file_list, out):
    '''
    Takes list of dataframe that are correctly formatted and 
    concatentates them into a single dataframe
    '''
    final_frame = pd.concat(file_list)
    final_frame = final_frame.sort_values(by=['chromosome'])
    final_frame['statistics_imputation_score'] = 'NA'
    out = out + '.csv.gz'
    return final_frame.to_csv(out, compression = 'gzip')

# final_frame = concatenate_all_tables(file_list = list_of_files, out='heeeeey')



def add_data_to_existing(existing_summary, summary_to_add, type, out):
    '''
    Takes a new dataframe, correctly formatted and/or cleaned 
    thereafter add to an existing dataframe.
    '''
    list_of_files = []
    existing_summary = pd.read_csv(existing_summary, engine='python', compression = 'infer', header = 0, sep=',')
    list_of_files.append(existing_summary)
    existing_summary_filename = existing_summary

    summary_to_add = pd.read_csv(summary_to_add, engine='python', compression = 'infer', header = 0, sep=',')
    list_of_files.append(summary_to_add)

    summary_to_add_filename = existing_summary  
    files_key = dict(zip(summary_to_add_filename, list_of_files))

    new_summary_checked_columns = check_columns(file_list = list_of_files, filenames = summary_to_add_filename, dataset = type, files_key = files_key)
    concatenate_all_tables(file_list = list_of_files, out = out)
    return 

# final_frame = add_data_to_existing(existing_summary = '/Users/ibrahim/Documents/Novo Nordisk/Combinator/ppp.csv', summary_to_add = '/Users/ibrahim/Documents/Novo Nordisk/Combinator/data/to_append/chr6_baso.csv', type = 'gwas', out = 'Hi')


