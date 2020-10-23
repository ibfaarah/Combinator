# Import libaries
import os
import glob
import re
import csv
import numpy as np 
import pandas as pd

# path = '/Users/ibrahim/Documents/Novo Nordisk/Combinator/test_compressed'
# basename = os.path.basename(path)

def extract_column_names(path, type):
    '''
    Takes a text file as input. This will contain a list comma separated values.
    This will be the default header for all concatenated files.
    Input file can be used a summary statistic that will be used in the concatenation. Or, can be a text file containing just a list.
    Output will be a list of column names for a given dataset 
    '''
  
    if type == 'csv':
        column_file = pd.read_csv(path, compression = 'infer')
        column_list = list(column_file)
    
    elif type == 'txt':
        file = open(path, 'r')
        column_list = file.read().split(',')
        
    # df.columns = list_of_names // using this change column names in a dataframe using list

    return column_list
# print(extract_column_names(path = '/Users/ibrahim/Documents/Novo_Nordisk/Combinator/data/chr1_baso.csv', type = 'csv'))

def read_compressed_files(path, compression):
    '''
    Takes 'path' to directory containing all the datasets.
    

    Tip: 
    (1) Add 'error_bad_lines=False' to pandas method.
    (2) When implementing this CycleCloud change '/*.csv' to '/*.gz'
    '''
    file_list = []
    
    if compression is None:
        directory = glob.glob(path + "/*.csv")

        for file in directory:
            df = pd.read_csv(file, engine='python', header = 0, sep=',')
            file_list.append(df)

    if compression == True:
        directory = glob.glob(path + "/*.csv.gz")

        for file in directory:
            df = pd.read_csv(file, engine='python', compression = 'gzip', header = 0, sep=',')
            file_list.append(df)

            if not file_list:
                raise ValueError('For files that are (.csv.gz) please use --gzip flag')


    return file_list


# list_of_files = read_compressed_files(path=path, compression='gzip')
# print(list_of_files)



def read_compressed_filenames(path, compression=None):
    '''
    Takes 'path' to directory containing all the datasets.
    Whil

    Tip: 
    (1) Add 'error_bad_lines=False' to pandas method.
    (2) When implementing this CycleCloud change '/*.csv' to '/*.gz'
    '''
    if compression is None:
        filenames = glob.glob(path + "/*.csv")
    elif compression == True:
        filenames = glob.glob(path + "/*.csv.gz")
    
    return filenames


# filenames = read_compressed_filenames(path=path, compression='gzip')
# print(filenames)

# files_key = dict(zip(filenames, list_of_files))
# print(files_key)



def check_headers(file_list, filenames, files_key, dataset, column_names, strict):

    dataframe_to_remove= []
    dataframe_to_warn = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        if strict == True:
            # if dataframes.columns.values.tolist() != column_names:
            if dataframes.columns.str.strip().tolist() != column_names:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in header of' + " " + filenames[x])    
                    print('Please check requirements follow required column names')
                    print('This dataset has now been removed..')
                    del file_list[x]
        
        elif strict == False:
            
            if dataframes.columns.str.strip().tolist() != column_names:
                template_column = len(column_names)
                number_of_columns = len(dataframes.columns.str.strip().tolist())
                
                if number_of_columns == template_column:
                    dataframe_to_warn.append(index)
                    for x in sorted(dataframe_to_warn, reverse=True):
                        print('\nWARNING: Header mismatch found in' + " " + filenames[x])
                        print('File will be concatenated, however, beware of potential errors')
                
                else:
                    dataframe_to_remove.append(index)
                    for x in sorted(dataframe_to_remove, reverse=True):
                        print('\nERROR found in header of' + " " + filenames[x])
                        print('Please check requirements follow required column names')
                        print('This dataset has now been removed..')
                        del file_list[x]

    return file_list




def check_columns(file_list, filenames, dataset, files_key):
    '''
    Checks all that data type for all columns are as expected:
    Use: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.dtypes.html
    '''

    logging_df = []
    error = ['NA', np.nan]

    chr_list = ['1', '2', '3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y', 'MT']
    directions = ['+','-']
    # stats_imp_score = ['NA']


    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):

        dataframe_to_remove= []
        # print(file_list)
        # print(filenames)
        # This is for 'GWAS summary statistics data
        
        # gwas

        if dataset == 'gwas':
            
            '''
            chromosome - ~dataframes.iloc[:,[2]]
            position - dataframes.iloc[:,[3]]
            effect_allele_frequency - dataframes.iloc[:,[8]]
            minor_allele_frequency - dataframes.iloc[:,[9]]
            standard_error - dataframes.iloc[:,[11]]
            p - dataframes.iloc[:,[12]]
            strand - dataframes.iloc[:,[7]]
            direction - dataframes.iloc[:,[15]]
            original_strand - dataframes.iloc[:,[25]]
            original_direction - dataframes.iloc[:,[26]]
            '''

            # Checks that all entries in 'chromosome' column are correct, 1-22 & X, Y, MT
            if dataframes[dataframes.columns[2]].astype(str).isin(chr_list).all() == False:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                    print('\n This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass


            # Checks that all entries in 'position' column have correct data type, only integers (int64)
            if dataframes[dataframes.columns[3]].dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

        
            # Checks that all entries in 'effect allele frequency' column has values between (0 & 1)
            if dataframes[dataframes.columns[8]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[8]].between(0,1).all() != True:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            # Checks that all entries in 'minor allele frequency' column has values between (0 & 0.5)    

            if dataframes[dataframes.columns[9]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[9]].between(0,0.5).all() != True:
                dataframe_to_remove.append(index)
        
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
        
            else:
                pass


            # Checks that all entries in 'Standard error' column has values greater than (0)    

            if dataframes[dataframes.columns[11]].all() <= 0:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
                print('Standard error column, has incorrect entry')

            else:
                pass

            # Checks that all entries in 'p-value' column has values between (0 & 1)    
            if dataframes[dataframes.columns[12]].between(0,1).any() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass


            # Checks that all entries in 'strand' column are either ('+' or '-' )
            if dataframes[dataframes.columns[7]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            # Checks that all entries in 'direction' column are either ('+' or '-' )
            if dataframes[dataframes.columns[15]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
    
            # Checks that all entries in 'original strand' column are either ('+' or '-' )
            if dataframes[dataframes.columns[25]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
    
            # Checks that all entries in 'original direction' column are either ('+' or '-' )
            if dataframes[dataframes.columns[26]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
        

        # eqtl

        if dataset == 'eQTL':
        
         
            '''
            chromosome - ~dataframes.iloc[:,[4]]
            position - dataframes.iloc[:,[5]]
            effect_allele_frequency - dataframes.iloc[:,[10]]
            minor_allele_frequency - dataframes.iloc[:,[11]]
            standard_error - dataframes.iloc[:,[13]]
            p - dataframes.iloc[:,[15]]
            strand - dataframes.iloc[:,[9]]
            direction - dataframes.iloc[:,[17]]
            original_strand - dataframes.iloc[:,[27]]
            original_direction - dataframes.iloc[:,[28]]
            '''
        
            # Checks that all entries in 'chromosome' column are correct, 1-22 & X, Y, MT
            if dataframes[dataframes.columns[4]].astype(str).isin(chr_list).all() == False:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                    print('\n This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass


            if dataframes[dataframes.columns[5]].dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[10]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[10]].between(0,1).all() != True:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[11]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[11]].between(0,0.5).all() != True:
                dataframe_to_remove.append(index)
        
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
        
            else:
                pass   

            if dataframes[dataframes.columns[11]].all() <= 0:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
                print('Standard error column, has incorrect entry')

            else:
                pass

            if dataframes[dataframes.columns[15]].between(0,1).any() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass


            if dataframes[dataframes.columns[9]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[17]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[27]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
    

            if dataframes[dataframes.columns[28]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

        # mqtl

        if dataset == 'mQTL':

            '''
            chromosome - ~dataframes.iloc[:,[3]]
            position - dataframes.iloc[:,[4]]
            effect_allele_frequency - dataframes.iloc[:,[9]]
            minor_allele_frequency - dataframes.iloc[:,[10]]
            standard_error - dataframes.iloc[:,[12]]
            p - dataframes.iloc[:,[14]]
            strand - dataframes.iloc[:,[8]]
            direction - dataframes.iloc[:,[16]]
            original_strand - dataframes.iloc[:,[26]]
            original_direction - dataframes.iloc[:,[27]]
            '''
            # Checks that all entries in 'chromosome' column are correct, 1-22 & X, Y, MT
            if dataframes[dataframes.columns[3]].astype(str).isin(chr_list).all() == False:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                    print('\n This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass


            if dataframes[dataframes.columns[4]].dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[9]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[9]].between(0,1).all() != True:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[10]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[10]].between(0,0.5).all() != True:
                dataframe_to_remove.append(index)
        
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
        
            else:
                pass   

            if dataframes[dataframes.columns[12]].all() <= 0:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
                print('Standard error column, has incorrect entry')

            else:
                pass

            if dataframes[dataframes.columns[14]].between(0,1).any() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[8]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[16]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[26]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
    

            if dataframes[dataframes.columns[27]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
        
        # pqtl
        if dataset == 'pQTL':

            '''
            chromosome - ~dataframes.iloc[:,[4]]
            position - dataframes.iloc[:,[5]]
            effect_allele_frequency - dataframes.iloc[:,[10]]
            minor_allele_frequency - dataframes.iloc[:,[11]]
            standard_error - dataframes.iloc[:,[13]]
            p - dataframes.iloc[:,[14]]
            strand - p - dataframes.iloc[:,[9]]
            direction - dataframes.iloc[:,[17]]
            original_strand - dataframes.iloc[:,[27]] 
            original_direction - dataframes.iloc[:,[28]] 
            '''

            if dataframes[dataframes.columns[4]].astype(str).isin(chr_list).all() == False:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                    print('\n This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass


            if dataframes[dataframes.columns[5]].dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[10]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[10]].between(0,1).all() != True:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[11]].astype(str).isin(error).all() != True:
                pass
            elif dataframes[dataframes.columns[11]].between(0,0.5).all() != True:
                dataframe_to_remove.append(index)
        
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
        
            else:
                pass   

            if dataframes[dataframes.columns[13]].all() <= 0:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
                print('Standard error column, has incorrect entry')

            else:
                pass

            if dataframes[dataframes.columns[14]].between(0,1).any() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[9]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[17]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[27]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
    

            if dataframes[dataframes.columns[28]].astype(str).isin(directions).all() != True:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass


        # variant

        if  dataset == 'variant':

            '''
            chromosome - dataframes.iloc[:,[0]]
            position - dataframes.iloc[:,[1]]
            cds_position - dataframes.iloc[:,[9]]
            cdna_position - dataframes.iloc[:,[8]]
            '''

            if dataframes[dataframes.columns[0]].astype(str).isin(chr_list).all() == False:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                    print('\n This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
    
            if dataframes[dataframes.columns[1]].dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[9]].dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "cds position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass

            if dataframes[dataframes.columns[8]].dtype != np.int64:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "cDNA position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
        

        # allele_MAP
        
        if  dataset == 'allele_map':

            '''
            chromosome - dataframes.iloc[:,[0]]
            position - dataframes.iloc[:,[1]]
            '''           

            if dataframes[dataframes.columns[0]].astype(str).isin(chr_list).all() == False:
                dataframe_to_remove.append(index)
                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                    print('\n This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass
    
            if dataframes[dataframes.columns[1]].dtype != np.int64:
                dataframe_to_remove.append(index)

                for x in sorted(dataframe_to_remove, reverse=True):
                    print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                    print('This dataset has now been removed..')
                    del file_list[x]
                    del filenames[x]
            else:
                pass  
            
            
    return file_list

# check_columns = check_columns(file_list = list_of_files, filenames = filenames, dataset='gwas', files_key=files_key)
# print(check_columns)



def concatenate_all_tables(file_list, column_list, out):
    '''
    Takes list of dataframe that are correctly formatted and 
    concatentates them into a single dataframe
    '''
    final_frame = pd.concat(file_list)
    final_frame = final_frame.sort_values(by=['chromosome'])
    final_frame['statistics_imputation_score'] = 'NA'
    final_frame.columns = column_list
    out = out + '.csv.gz'
    return final_frame.to_csv(out, index=False, compression = 'gzip')

# final_frame = concatenate_all_tables(file_list = list_of_files, out='heeeeey')



def add_data_to_existing(existing_summary, summary_to_add, compression, dataset, set_header, strict, out):
    '''
    Takes a new dataframe, correctly formatted and/or cleaned 
    thereafter add to an existing dataframe.
    '''
    list_of_files = []
    existing_summary_csv = existing_summary
    
    if set_header.endswith('.csv'):
        column_list = extract_column_names(path = set_header, type = 'csv')
    else:
        column_list = extract_column_names(path = set_header, type = 'txt')
    
    
    if existing_summary_csv.endswith('.csv.gz'):
        existing_summary = pd.read_csv(existing_summary, engine='python', compression = 'gzip', header = 0, sep=',')
    else:
        existing_summary = pd.read_csv(existing_summary, engine='python', compression = 'infer', header = 0, sep=',')



    if compression is None:
        summary_to_add_dataframe = pd.read_csv(summary_to_add, engine='python', compression = 'infer', header = 0, sep=',')
        list_of_files.append(summary_to_add_dataframe)

    else:
        summary_to_add_dataframe = pd.read_csv(summary_to_add, engine='python', compression = 'gzip', header = 0, sep=',')
        list_of_files.append(summary_to_add_dataframe)
    
    if not list_of_files:
        raise ValueError('\nFor files that are (.csv.gz) please use --gzip flag. For (.csv) files no flag required.')
    
    summary_file_name = []
    basename = os.path.basename(summary_to_add)
    summary_file_name.append(basename)
    files_key = dict(zip(summary_file_name, list_of_files))

    list_of_files = check_headers(file_list = list_of_files, filenames = summary_file_name, dataset = dataset, files_key = files_key, column_names=column_list, strict= strict)
    list_of_files = check_columns(file_list = list_of_files, filenames = summary_file_name, dataset = dataset, files_key = files_key)

    if not list_of_files:
        print('\nCould not add' + " " + basename + " " + "..." + " " + " Please fix error and try again...")
    else:
        list_of_files.append(existing_summary)
        final_frame = pd.concat(list_of_files)
        final_frame = final_frame.sort_values(by=['chromosome'])
        final_frame.columns = column_list
        out = out + '.csv.gz'
        final_frame.to_csv(out, index=False, compression = 'gzip')
        print('....Successful addition of' + " " + os.path.basename(summary_to_add) + " " + "to" + " " + os.path.basename(existing_summary_csv) + ".")
    return 

# final_frame = add_data_to_existing(existing_summary = '/Users/ibrahim/Documents/Novo Nordisk/Combinator/ppp.csv', summary_to_add = '/Users/ibrahim/Documents/Novo Nordisk/Combinator/data/to_append/chr6_baso.csv', type = 'gwas', out = 'Hi')





