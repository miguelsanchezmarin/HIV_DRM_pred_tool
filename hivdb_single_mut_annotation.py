import pandas as pd
import re

##Drugs per drug class used in this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']

def HIVDB_single_mut(mut_dic : dict, hivdb_data_dic: dict, cutoff: float):
    '''Generates a dictionary with the mutations and their HIVDB annotations for each drug class.
    INPUT:
    mut_dic: dictionary with the mutations and their frequencies. The keys are the mutations and the values are dictionaries with the frequencies.
    hivdb_data_dic: dictionary with the HIVDB data. The keys are the drug classes and the values are dictionaries with the mutations and their annotations.
    cutoff: minimum frequency to be considered as a mutation.
    
    OUTPUT:
    output: dictionary with the mutations and their HIVDB annotations and comments for each drug class.'''

    output = dict()
    output["cutoff"] = cutoff

    for prot in ["RT", "PR", "IN"]:
        output[prot] = dict()

    for dataset in ['INI', 'NNRTI', 'NRTI', 'PI']:

        if dataset == "INI":
            protein = 'IN'
        elif dataset == "NNRTI":
            protein = 'RT'
        elif dataset == "NRTI":
            protein = 'RT'
        elif dataset == "PI":
            protein = 'PR'

        #we filter the mutations by frequency
        mut_dic[protein] = freq_filter(mut_dic[protein], cutoff)

        for mutation in mut_dic[protein].keys():
            frequency = mut_dic[protein][mutation]["Freq"]

            #we add the single mutation annotations
            output[protein][mutation] = dict()
            output[protein][mutation]['Freq'] = frequency #save the frequency
            position, aa = mutation[:-1], mutation[-1]
            comment_in = False
            for comm in hivdb_data_dic[dataset]['comments'].keys():
                if re.match(r'^\d+[A-Z]+$', comm) and position==re.sub(r'\D+', '', comm):
                    comm_ambig = re.sub(r'\d+', '', comm)
                    if aa in comm_ambig:
                        output[protein][mutation]['Annotation'] = hivdb_data_dic[dataset]['comments'][comm]['mut_type'] #save HIVDB annotations
                        output[protein][mutation]['Comment'] = hivdb_data_dic[dataset]['comments'][comm]['comment']
                        comment_in = True
                        break
            if not comment_in:
                output[protein][mutation]['Annotation'] = "Unknown"
                output[protein][mutation]['Comment'] = "Unknown"

    return output


def freq_filter(mut_dict: dict, cutoff: float):
    '''Filters the mutations in the input dictionary by frequency.
    
    INPUT:
    mut_dict: dictionary with the mutations as keys and the frequencies as values.
    cutoff: float with the frequency cutoff value.
    
    OUTPUT:
    filtered_mut_dict: dictionary with the filtered mutations.
    '''
    filtered_mut_dict = mut_dict.copy()
    del_muts = []
    for mut in filtered_mut_dict.keys():
        if filtered_mut_dict[mut]['Freq'] < cutoff:
            del_muts.append(mut)
    #we remove it
    for mut in del_muts:
        del filtered_mut_dict[mut]

    return filtered_mut_dict

def get_single_mut_annotations_df(single_annotations_dic_list:list, sample_ids:list = None):
    '''Creates a pandas DataFrame from the single mutation annotations dictionary list of samples.
    INPUT:
    single_annotations_dic_list: list of dictionaries with the single mutation annotations for each sample.
    sample_ids: list of sample ids.
    OUTPUT:
    single_mut_df: pandas DataFrame with the single mutation annotations for each sample.
    '''
    df = []
    if sample_ids is None:
        cutoff = single_annotations_dic_list['cutoff']
        for prot in single_annotations_dic_list.keys():
            if prot != 'cutoff':
                #we sort the mutations by position
                sorted_muts = sorted(single_annotations_dic_list[prot].keys(), key=lambda x: int(re.sub(r'\D+', '', x)))
                for mut in sorted_muts:
                    mut_data = single_annotations_dic_list[prot][mut]
                    df.append([prot, mut, mut_data['Freq'], mut_data['Annotation'], mut_data['Comment'], cutoff])
        single_mut_df = pd.DataFrame(df, columns=['Protein', 'Mutation', 'Freq', 'Annotation', 'Comment', 'Cutoff'])
    else:
        for sample_ID, sample_annot in zip(sample_ids, single_annotations_dic_list):
            cutoff = sample_annot['cutoff']
            for prot in sample_annot.keys():
                if prot != 'cutoff':
                    #we sort the mutations by position
                    sorted_muts = sorted(sample_annot[prot].keys(), key=lambda x: int(re.sub(r'\D+', '', x)))
                    for mut in sorted_muts:
                        mut_data = sample_annot[prot][mut]
                        df.append([sample_ID, prot, mut, mut_data['Freq'], mut_data['Annotation'], mut_data['Comment'], cutoff])
        
        single_mut_df = pd.DataFrame(df, columns=['Sample_ID', 'Protein', 'Mutation', 'Freq', 'Annotation', 'Comment', 'Cutoff'])
    single_mut_df = single_mut_df.drop_duplicates()

    return single_mut_df