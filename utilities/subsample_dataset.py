###Randomly subsample increasing missing percentages of sequences of the original dataset
import pandas as pd
import numpy as np
import re


test_drugs = ['RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
INIs=['RAL', 'EVG', 'DTG', 'BIC']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
drug_classes = ['INI', 'PI', 'NNRTI', 'NRTI']

pi_major_positions = [30, 32, 47, 48, 50, 54, 76, 82, 84, 88]
nrti_major_positions = [41, 65, 70, 74, 75, 151, 184, 210, 215]
nnrti_major_positions = {100: ['I'], 101: ['P'], 103: ['N'], 106: ['A','M'], 181: ['C','I','A'], 188: ['C','L','H'], 190: ['A','E','S','Q'], 230: ['L']}
ini_major_positions = {66: ['A','I','K'], 92: ['Q'], 118: ['R'], 143: [r'[A-Z]'], 148: ['H','R','K'], 155: ['H'], 263: ['K']}
nnrti_major_positions = [100, 101, 103, 106, 181, 188, 190, 230]
ini_major_positions = [66, 92, 118, 143, 148, 155, 263]

percentage_sequence = [0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.0]


for cl in drug_classes:

    if cl == "PI":
        major_positions = pi_major_positions
        seq_length = 99
    elif cl == "NRTI":
        major_positions = nrti_major_positions
        seq_length = 240
    elif cl == "NNRTI":
        major_positions = nnrti_major_positions
        seq_length = 240
    elif cl == "INI":
        major_positions = ini_major_positions
        seq_length = 288

    full_data_df = pd.DataFrame()
    for prc in percentage_sequence:
        for it in range(1): #iterations just in case we want to add at some point
            # print(prc, '%', it, 'it')
            data_set = pd.read_csv(f"../datasets/{cl}_dataset.tsv", sep='\t')
            ##we read the CompMutList column
            comp_mut_list = data_set['CompMutList'].tolist()

            comp_mut_list_keep, major_miss, minor_miss, keep_mutations_list = [], [], [], []
            for mut_list in comp_mut_list:

                #we randomly select a prc% of the positions of the sequence
                keep_positions = np.random.choice(range(1, seq_length+1), int(prc*seq_length), replace=False)
                keep_positions_str = ', '.join(np.sort(keep_positions).astype(str))

                if pd.isnull(mut_list):
                    comp_mut_list_keep.append('')
                    major_miss.append('')
                    minor_miss.append('')
                    keep_mutations_list.append(keep_positions_str)
                    continue


                mut_list = mut_list.split(', ')
                mut_list_keep, mut_list_miss_major, mut_list_miss_minor = [], [], []

                for mut in mut_list:
                    position = int(re.search(r'\d+', mut).group())
                    if position in keep_positions:
                        mut_list_keep.append(mut.strip())
                    elif position in major_positions:
                        mut_list_miss_major.append(mut.strip())
                    else:
                        mut_list_miss_minor.append(mut.strip())

                comp_mut_list_keep.append(', '.join(mut_list_keep))
                major_miss.append(', '.join(mut_list_miss_major))
                minor_miss.append(', '.join(mut_list_miss_minor))
                keep_mutations_list.append(keep_positions_str)

            data_set[f'CompMutList_keep'] = comp_mut_list_keep
            data_set[f'Miss_major'] = major_miss
            data_set[f'Miss_minor'] = minor_miss            
            #we add a colum with the iteration number and other with the seq_percent
            data_set[f'keep_seq_positions'] = keep_mutations_list
            data_set['iteration'] = [it for _ in range(len(comp_mut_list))]
            data_set['seq_percent'] = [prc for _ in range(len(comp_mut_list))]
            #we concatenate
            full_data_df = pd.concat([full_data_df, data_set], axis=0)
        
            print("Done")
    
    full_data_df.to_csv(f'../datasets/subsampled_datasets/{cl}_subsampled_dataset.tsv', sep='\t', index=False)
