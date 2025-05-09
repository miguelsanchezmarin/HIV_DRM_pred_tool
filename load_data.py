##Functions for loading input data and prediction models as dictionaries

import pandas as pd
import os
import pickle

##Drug used for this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']

def load_sample_data(pathtodir:str = 'input_samples'):
    '''Loads the sample data coverage.tsv and mutation_freq.tsv from the given directory in a single dictionary.'''
    
    sample_data = dict()
    for samples in os.listdir(pathtodir):
        mut_freq_data = pd.read_csv(f'{pathtodir}/{samples}/mutation_freq.tsv', sep='\t')
        coverage_data = pd.read_csv(f'{pathtodir}/{samples}/coverage_annotated.tsv', sep='\t')

        #we filter out mutations with Ref == Mut and with Mut == *
        mut_freq_data = mut_freq_data[mut_freq_data['Ref'] != mut_freq_data['Mut']]
        mut_freq_data = mut_freq_data[mut_freq_data['Mut'] != '*']

        sample_data[samples] = dict()
        for protein in ['PR', 'RT', 'IN']:
            mut_freq_data_protein = mut_freq_data[mut_freq_data['Prot'] == protein]
            coverage_data_protein = coverage_data[coverage_data['Protein'] == protein]
            sample_data[samples][protein] = dict()
            for position in mut_freq_data_protein['Position'].values:
                sample_data[samples][protein][position] = dict()
                sample_data[samples][protein][position]['Ref'] = mut_freq_data_protein[mut_freq_data_protein['Position'] == position]['Ref'].values[0]
                sample_data[samples][protein][position]['Mut'] = mut_freq_data_protein[mut_freq_data_protein['Position'] == position]['Mut'].values[0]
                sample_data[samples][protein][position]['Freq'] = mut_freq_data_protein[mut_freq_data_protein['Position'] == position]['Freq'].values[0]
                sample_data[samples][protein][position]['Coverage'] = coverage_data_protein[coverage_data_protein['CodonPosition'] == position]['Coverage'].values[0]
            
    return sample_data


def load_hivdb_data(pathtodir:str = 'HIVDB_rules'):
    '''Loads the HIVDB data from the given directory in a single dictionary.'''
    
    hivdb_data = dict()
    for dataset in ['INI', 'NNRTI', 'NRTI', 'PI']:
        hivdb_data[dataset] = dict()
        #we load the single muts score
        hivdb_data[dataset]['single_mut'] = dict()
        hivdb_muts = pd.read_csv(f'{pathtodir}/{dataset}_muts_score_Stanford_HIVDB', sep=',')
        for rule in hivdb_muts['Rule'].values:
            hivdb_data[dataset]['single_mut'][rule] = dict()
            for drug in hivdb_muts.columns[1:]:
                hivdb_data[dataset]['single_mut'][rule][drug.replace("/r", "")] = hivdb_muts[hivdb_muts['Rule'] == rule][drug].values[0]

        #we load the combinations score
        hivdb_data[dataset]['combination'] = dict()
        hivdb_comb = pd.read_csv(f'{pathtodir}/{dataset}_combinations_score_Stanford_HIVDB', sep=',')
        for rule in hivdb_comb['Combination Rule'].values:
            hivdb_data[dataset]['combination'][rule.strip()] = dict()
            for drug in hivdb_comb.columns[1:]:
                hivdb_data[dataset]['combination'][rule.strip()][drug.replace("/r", "")] = hivdb_comb[hivdb_comb['Combination Rule'] == rule][drug].values[0]
        
        #we load the comments
        hivdb_data[dataset]['comments'] = dict()
        hivdb_comments = pd.read_csv(f'{pathtodir}/{dataset}_comments_Stanford_HIVDB', sep=',')
        hivdb_comments.columns = hivdb_comments.columns.str.replace('\n','')
        for condition in hivdb_comments['Condition'].values:
            hivdb_data[dataset]['comments'][condition] = dict()
            for drug in hivdb_comments.columns[1:]:
                hivdb_data[dataset]['comments'][condition]["mut_type"] = hivdb_comments[hivdb_comments['Condition'] == condition]["Comment/Mutation Type"].values[0]
                hivdb_data[dataset]['comments'][condition]["comment"] = hivdb_comments[hivdb_comments['Condition'] == condition]["Comment"].values[0]

    return hivdb_data

def load_lsr_coef(pathtodir:str = 'linear_regression_coefficients'):
    '''Loads the LSR coefficients from the given directory in a single dictionary.'''
    
    lsr_coef = dict()
    for drug in INIs+NNRTIs+NRTIs+PIs:
        lsr_coef[drug] = dict()
        lsr_coef_data = pd.read_csv(f'{pathtodir}/OLS_{drug}_combinations_tsm_all_folds.txt')
        for n in range(lsr_coef_data.shape[0]):
            if n == 0:
                lsr_coef[drug]['intercept'] = lsr_coef_data.iloc[n,0].split(' ')[1]
            else:
                coef_line = lsr_coef_data.iloc[n,0].split(' ')
                mutation = coef_line[0][2:]
                coef = coef_line[1]
                lsr_coef[drug][mutation] = float(coef)
    
    return lsr_coef

def load_random_forest(pathtodir:str = 'random_forest_models'):
    '''Loads the Random Forest models from the given directory in a single dictionary.'''
    
    rf_models = dict()
    for drug in INIs+NNRTIs+NRTIs+PIs:
        rf_models[drug] = pickle.load(open(f'{pathtodir}/random_forest_python_{drug}_RF_model_allfolds.pkl', 'rb'))
    
    return rf_models

# print(load_hivdb_data())
# print(load_lsr_coef())
# print(load_random_forest())
# print(load_sample_data())