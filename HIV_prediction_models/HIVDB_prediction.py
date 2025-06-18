###Drug Resistance Mutation prediction with HIVDB sierrapy prediction tool

###sierrapy should be installed from Github sierra-client (https://github.com/hivdb/sierra-client.git).

import pandas as pd
import os
import logging
import shlex
import subprocess
import json

def HIVDB_prediction(path_to_mutations_df):
    """Runs ''sierrapy patterns'' from a dataframe with seqID, gene, CompMutList columns.
    Then reads the output json and returns a dataframe with the results.

    path_to_mutations_df: absolute path to the mutations dataframe

    """
    mutations_df = path_to_mutations_df

    seq_count = 0
    
    with open('pattern.txt', 'w') as h:   
        for n in range(mutations_df.shape[0]):
            seqID, muts, gene = str(mutations_df.iloc[n]['SeqID']), str(mutations_df.iloc[n]['CompMutList']), mutations_df.iloc[n]['gene']
            if muts == 'nan':   #we skip the sequences with no mutations
                print(f'Skipping {seqID}, no mutations')
                continue
            mut_list = muts.split(', ')
            mut_list = [ gene + ':'+mut for mut in mut_list]
            muts = ' + '.join(mut_list)
            h.write('>' + seqID + '\n' + muts + '\n')
            seq_count += 1
    
    #we first check how many .json files will be generated (1 every 100 rows)    
    if seq_count % 100 == 0:
        n_files = seq_count//100
    else:
        n_files = seq_count//100 + 1    

    print(f'{n_files} json files will be generated')

    cml = shlex.split('sierrapy patterns pattern.txt -o o.json')    #we run sierrapy patterns
    logging.debug(cml)
    subprocess.call(cml)

    #we read over the json files
    resistance_results = [] #we will store the results in a list of dictionaries

    for i in range(n_files):
        
        with open(f'o.{i}.json') as h:  #we read the json output
            patterns = json.load(h)
    
        for seq in patterns:
            seq_results = {}#we create a dictionary with the results for each sequence
            seqID = seq['name'] #we fetch the seqID
            seq_results['seqID'] = seqID
            for dr in seq['drugResistance']:
                gene_name = dr['gene']['name']
                for drugscore in dr['drugScores']:
                    drug = drugscore['drug']['name']
                    score = drugscore['score']
                    label = drugscore['text'] #we get the classification label
                    seq_results[drug] = label
            resistance_results.append(seq_results)
        
        os.remove(f'o.{i}.json')#we remove the json file
    
    resistance_df = pd.DataFrame(resistance_results)    #we create a dataframe with the results
    os.remove('pattern.txt')

    return resistance_df


###We load the datasets
ini_dataset = pd.read_csv("../datasets/INI_dataset.tsv", sep='\t')
nnrti_dataset = pd.read_csv("../datasets/NNRTI_dataset.tsv", sep='\t')
nrti_dataset = pd.read_csv("../datasets/NRTI_dataset.tsv", sep='\t')
pi_dataset = pd.read_csv("../datasets/PI_dataset.tsv", sep='\t')

ini_dataset['gene'] = 'IN'
nnrti_dataset['gene'] = 'RT'
nrti_dataset['gene'] = 'RT'
pi_dataset['gene'] = 'PR'


###We predict the resistance of all four databases and save the results
ini_preds = HIVDB_prediction(ini_dataset)
ini_preds.to_csv('../method_predictions/HIVDB/HIVDB_INI_predictions.tsv', sep='\t')

nnrti_preds = HIVDB_prediction(nnrti_dataset)
nnrti_preds.to_csv('../method_predictions/HIVDB/HIVDB_NNRTI_predictions.tsv', sep='\t')

nrti_preds = HIVDB_prediction(nrti_dataset)
nrti_preds.to_csv('../method_predictions/HIVDB/HIVDB_NRTI_predictions.tsv', sep='\t')

pi_preds = HIVDB_prediction(pi_dataset)
pi_preds.to_csv('../method_predictions/HIVDB/HIVDB_PI_predictions.tsv', sep='\t')
