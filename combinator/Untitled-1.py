
import pandas as pd
import numpy as np
import glob
import re

#Conditions
error = ['NA', np.nan]
directions = ['-','+']
chr_list = ['1', '2', '3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y', 'MT']

#Individual files
# file = '/Users/ibrahim/Documents/Novo_Nordisk/Combinator/example/data/chr2_baso.csv'
# dataframes = pd.read_csv(file, keep_default_na=False, engine='python', header = 0, sep=',')

# if dataframes[dataframes.columns[7]].astype(str).isin(directions).all() != True:
#     print('ERROR found in "Strand column')    
# else:
#     print('no error')
#     pass

#List of files 

path = '/Users/ibrahim/Documents/Novo_Nordisk/Combinator/example/data'
directory =  glob.glob(path + "/*.csv")
filenames = glob.glob(path + "/*.csv")

file_list = []



for file in directory:
            df = pd.read_csv(file, keep_default_na=False, engine='python', header = 0, sep=',')
            file_list.append(df)








# For all datasets

def empty_columns(file_list, filenames):
    '''
    Function takes list of files and checks whether
    the column '7' contains only '+' or '-' 
    - removing it when it doesnt
    '''

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        dataframes[dataframes.columns[-1]].replace({None:'Empty_values'}, inplace=True)
        if dataframes[dataframes.columns[-1]].isin(['Empty_values']).any().any() == True:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR, rows are unequal in length in \"'+ " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

# file_list, filenames = empty_columns(file_list=file_list, filenames=filenames)

# Individual GWAS functions
def gwas_chromosome(file_list, filenames):
    '''
    Function takes list of files and checks whether
    the column '15' contains only '+' or '-' 
    - removing it when it doesnt
    '''

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[2]].astype(str).isin(chr_list).all() == False:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

# file_list, filenames = gwas_chromosome(file_list=file_list, filenames=filenames)

def gwas_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[3]].dtype != np.int64:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

# file_list, filenames = gwas_position(file_list=file_list, filenames=filenames)

def gwas_strand(file_list, filenames):
    '''
    Function takes list of files and checks whether
    the column '7' contains only '+' or '-' 
    - removing it when it doesnt
    '''

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[7]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

# file_list, filenames = gwas_strand(file_list=file_list, filenames=filenames)


def gwas_direction(file_list, filenames):
    '''
    Function takes list of files and checks whether
    the column '15' contains only '+' or '-' 
    - removing it when it doesnt
    '''

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[15]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Direction column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

# file_list, filenames =gwas_direction(file_list=file_list, filenames=filenames)


def gwas_effect_allele_frequency(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[8]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[8]].between(0,1).all() != True:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

# file_list, filenames = gwas_effect_allele_frequency(file_list=file_list, filenames=filenames)

def gwas_minor_allele_frequency(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[9]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[9]].between(0,0.5).all() != True:
            dataframe_to_remove.append(index)
    
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

# file_list, filenames = gwas_minor_allele_frequency(file_list=file_list, filenames=filenames)


# Main function for gwas

def check_columns(file_list, filenames, dataset):


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

        file_list, filenames = empty_columns(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_chromosome(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_direction(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_effect_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_minor_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_standard_error(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_p_value(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_original_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = gwas_original_direction(file_list=file_list, filenames=filenames)
    
    if dataset == 'eqtl':

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
        
        file_list, filenames = empty_columns(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_chromosome(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_position(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_effect_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_minor_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_standard_error(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_p_value(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_direction(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_original_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = eqtl_original_direction(file_list=file_list, filenames=filenames)

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

        file_list, filenames = empty_columns(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_chromosome(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_position(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_effect_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_minor_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_standard_error(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_p_value(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_direction(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_original_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = mqtl_original_direction(file_list=file_list, filenames=filenames)

        

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

        file_list, filenames = empty_columns(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_chromosome(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_position(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_effect_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_minor_allele_frequency(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_standard_error(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_p_value(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_direction(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_original_strand(file_list=file_list, filenames=filenames)
        file_list, filenames = pqtl_original_direction(file_list=file_list, filenames=filenames)    
    
    if dataset == 'variant':

        '''
        chromosome - dataframes.iloc[:,[0]]
        position - dataframes.iloc[:,[1]]
        cds_position - dataframes.iloc[:,[9]]
        cdna_position - dataframes.iloc[:,[8]]
        '''   

        file_list, filenames = empty_columns(file_list=file_list, filenames=filenames)     
        file_list, filenames = variant_chromosome(file_list=file_list, filenames=filenames)
        file_list, filenames = variant_position(file_list=file_list, filenames=filenames)
        file_list, filenames = variant_cds_position(file_list=file_list, filenames=filenames)
        file_list, filenames = variant_cdna_position(file_list=file_list, filenames=filenames)
        
    if dataset == 'allele_map':
        
        '''
        chromosome - dataframes.iloc[:,[0]]
        position - dataframes.iloc[:,[1]]
        '''
        
        file_list, filenames = empty_columns(file_list=file_list, filenames=filenames)
        file_list, filenames = allele_map_chromosome(file_list=file_list, filenames=filenames)
        file_list, filenames = allele_map_position(file_list=file_list, filenames=filenames)


    else:
        pass

    print(filenames)

x = check_columns(file_list=file_list, filenames=filenames, dataset='gwas')
print(x)  


def gwas_standard_error(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[11]].all() <= 0:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = gwas_standard_error(file_list=file_list, filenames=filenames)


def gwas_p_value(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[12]].between(0,1).any() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames =gwas_p_value(file_list=file_list, filenames=filenames)


def gwas_original_strand(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[25]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = gwas_original_strand(file_list=file_list, filenames=filenames)


def gwas_original_direction(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[26]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = gwas_original_direction(file_list=file_list, filenames=filenames)


# For eQTL datasets

def eqtl_chromosome(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[4]].astype(str).isin(chr_list).all() == False:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                print('\n This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = eqtl_chromosome(file_list=file_list, filenames=filenames)

def eqtl_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[3]].dtype != np.int64:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = eqtl_position(file_list=file_list, filenames=filenames)


def eqtl_effect_allele_frequency(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
    
        if dataframes[dataframes.columns[10]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[10]].between(0,1).all() != True:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = eqtl_effect_allele_frequency(file_list=file_list, filenames=filenames)

def eqtl_minor_allele_frequency(file_list, filenames):

    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):


        if dataframes[dataframes.columns[11]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[11]].between(0,0.5).all() != True:
            dataframe_to_remove.append(index)
    
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = eqtl_minor_allele_frequency(file_list=file_list, filenames=filenames)

def eqtl_standard_error(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):

        if dataframes[dataframes.columns[11]].all() <= 0:
            dataframe_to_remove.append(index)
            
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]
            print('Standard error column, has incorrect entry')

    return file_list, filenames

file_list, filenames = eqtl_standard_error(file_list=file_list, filenames=filenames)

def eqtl_p_value(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):

        if dataframes[dataframes.columns[15]].between(0,1).any() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]      

    return file_list, filenames

file_list, filenames = eqtl_p_value(file_list=file_list, filenames=filenames)

def eqtl_strand(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[9]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]  

    return file_list, filenames

file_list, filenames = eqtl_strand(file_list=file_list, filenames=filenames)

def eqtl_direction(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[27]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x] 

    return file_list, filenames

file_list, filenames = eqtl_direction(file_list=file_list, filenames=filenames)

def eqtl_original_strand(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[27]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x] 

    return file_list, filenames

file_list, filenames = eqtl_original_strand(file_list=file_list, filenames=filenames)

def eqtl_original_direction(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[28]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = eqtl_original_direction(file_list=file_list, filenames=filenames)


# For mQTL datasets

def mqtl_chromosome(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[3]].astype(str).isin(chr_list).all() == False:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                print('\n This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = mqtl_chromosome(file_list=file_list, filenames=filenames)

def mqtl_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[4]].dtype != np.int64:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = mqtl_position(file_list=file_list, filenames=filenames)

def mqtl_effect_allele_frequency(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
    
        if dataframes[dataframes.columns[9]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[9]].between(0,1).all() != True:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = mqtl_effect_allele_frequency(file_list=file_list, filenames=filenames)

def mqtl_minor_allele_frequency(file_list, filenames):

    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):


        if dataframes[dataframes.columns[10]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[10]].between(0,0.5).all() != True:
            dataframe_to_remove.append(index)
    
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = mqtl_minor_allele_frequency(file_list=file_list, filenames=filenames)

def mqtl_standard_error(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):

        if dataframes[dataframes.columns[12]].all() <= 0:
            dataframe_to_remove.append(index)
            
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]
            print('Standard error column, has incorrect entry')

    return file_list, filenames

file_list, filenames = mqtl_standard_error(file_list=file_list, filenames=filenames)

def mqtl_p_value(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):

        if dataframes[dataframes.columns[14]].between(0,1).any() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]      

    return file_list, filenames

file_list, filenames = mqtl_p_value(file_list=file_list, filenames=filenames)

def mqtl_strand(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[8]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]  

    return file_list, filenames

file_list, filenames = mqtl_strand(file_list=file_list, filenames=filenames)

def mqtl_direction(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[16]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Direction column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = mqtl_direction(file_list=file_list, filenames=filenames)

def mqtl_original_strand(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[26]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x] 

    return file_list, filenames

file_list, filenames = mqtl_original_strand(file_list=file_list, filenames=filenames)

def mqtl_original_direction(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[27]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = mqtl_original_direction(file_list=file_list, filenames=filenames)

# For pQTL datasets

def pqtl_chromosome(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[4]].astype(str).isin(chr_list).all() == False:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                print('\n This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = pqtl_chromosome(file_list=file_list, filenames=filenames)

def pqtl_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[5]].dtype != np.int64:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = pqtl_position(file_list=file_list, filenames=filenames)

def pqtl_effect_allele_frequency(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
    
        if dataframes[dataframes.columns[10]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[10]].between(0,1).all() != True:
            dataframe_to_remove.append(index)

            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Effect Allele Frequency column\"' + " "  + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = pqtl_effect_allele_frequency(file_list=file_list, filenames=filenames)

def pqtl_minor_allele_frequency(file_list, filenames):

    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):


        if dataframes[dataframes.columns[11]].astype(str).isin(error).all() != True:
            pass
        if dataframes[dataframes.columns[11]].between(0,0.5).all() != True:
            dataframe_to_remove.append(index)
    
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Minor Allele Frequency column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = pqtl_minor_allele_frequency(file_list=file_list, filenames=filenames)

def pqtl_standard_error(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):

        if dataframes[dataframes.columns[13]].all() <= 0:
            dataframe_to_remove.append(index)
            
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Standard error column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]
            print('Standard error column, has incorrect entry')

    return file_list, filenames

file_list, filenames = pqtl_standard_error(file_list=file_list, filenames=filenames)

def pqtl_p_value(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):

        if dataframes[dataframes.columns[14]].between(0,1).any() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "P-value column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]      

    return file_list, filenames

file_list, filenames = pqtl_p_value(file_list=file_list, filenames=filenames)

def pqtl_strand(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[9]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]  

    return file_list, filenames

file_list, filenames = pqtl_strand(file_list=file_list, filenames=filenames)

def pqtl_direction(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[17]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Direction column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x] 

    return file_list, filenames

file_list, filenames = pqtl_direction(file_list=file_list, filenames=filenames)

def pqtl_original_strand(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[26]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original strand column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x] 

    return file_list, filenames

file_list, filenames = pqtl_original_strand(file_list=file_list, filenames=filenames)

def pqtl_original_direction(file_list, filenames):
    
    dataframe_to_remove = []
    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[28]].astype(str).isin(directions).all() != True:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Original direction column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]

    return file_list, filenames

file_list, filenames = pqtl_original_direction(file_list=file_list, filenames=filenames)

# For variant datasets

def variant_chromosome(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[0]].astype(str).isin(chr_list).all() == False:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                print('\n This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = variant_chromosome(file_list=file_list, filenames=filenames)


def variant_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[1]].dtype != np.int64:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = variant_position(file_list=file_list, filenames=filenames)


def variant_cds_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[9]].dtype != np.int64:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "cds position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = variant_cds_position(file_list=file_list, filenames=filenames)


def variant_cdna_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[8]].dtype != np.int64:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "cDNA position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = variant_cdna_position(file_list=file_list, filenames=filenames)


# For allele MAP datasets

def allele_map_chromosome(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[0]].astype(str).isin(chr_list).all() == False:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\nERROR found in "Chromosome column\"' + " " + "in" + " " + filenames[x])    
                print('\n This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = allele_map_chromosome(file_list=file_list, filenames=filenames)


def allele_map_position(file_list, filenames):

    dataframe_to_remove = []

    for index, (filename, dataframes) in enumerate(zip(filenames, file_list)):
        
        if dataframes[dataframes.columns[1]].dtype != np.int64:
            dataframe_to_remove.append(index)
            for x in sorted(dataframe_to_remove, reverse=True):
                print('\n ERROR found in "Position column\"' + " " + "in" + " " + filenames[x])    
                print('This dataset has now been removed..')
                del file_list[x]
                del filenames[x]


    return file_list, filenames

file_list, filenames = allele_map_position(file_list=file_list, filenames=filenames)