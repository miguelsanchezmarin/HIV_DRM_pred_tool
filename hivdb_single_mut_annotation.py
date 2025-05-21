import pandas as pd
import re

##Drugs per drug class used in this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']

# def HIVDB_singlemut_annot(mut_list, dataset: str):
#     '''Retrieves the HIVDB annotation for each mutation in a list of mutations.
#     INPUT:
#     mut_list: list of mutations i.e [M41L,L90M].
#     dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.

#     OUTPUT:
#     mut_annot_df: dataframe with the mutations and their HIVDB annotations.
#     '''
#     #we read the comments files for the dataset
#     comments_file = pd.read_csv(f'HIVDB_rules/{dataset}_comments_Stanford_HIVDB', sep=',')
#     comments_file.columns = comments_file.columns.str.replace('\n','')

#     mut_list = [mut[1:] for mut in mut_list] #we remove the first letter of the mutations
#     mut_annot = []
#     for mut in mut_list:
#         position = re.sub(r'\D+', '', mut)
#         mut_aa = re.sub(r'\d+', '', mut)
#         comment_in = False
#         for i, comm in enumerate(comments_file['Condition'].values):
#             if re.match(r'^\d+[A-Z]+$', comm) and position==re.sub(r'\D+', '', comm): #we check if the mutation is in the comment
#                 comm_ambig = re.sub(r'\d+', '', comm)
#                 if mut_aa in comm_ambig:
#                     mut_annot.append([mut, 
#                                      comments_file['Comment/Mutation Type'].values[i], 
#                                      comments_file['Comment'].values[i]])
#                     comment_in = True
#                     break
        
#         if not comment_in:
#             mut_annot.append([mut, "Unknown", "Unknown"])
                       
#     mut_annot_df = pd.DataFrame(mut_annot, columns=['Mutation', 'Annotation', 'Comment'])
#     return mut_annot_df

# def HIVDB_table(mut_tsv, lower_cutoff: float = 0.015):
#     '''Generates a table with the observed mutations and their HIVDB classifications for all drugs and drug class.
#     INPUT:
#     mut_tsv: tsv file output from annotate_vcf.py with the observed mutations. Columns are Position, Ref, Mut, Freq and Prot. 
#     lower_cutoff: minimum frequency to be considered as a mutation.

#     OUTPUT:
#     mutations_df: dataframe with the mutations for all drugs and drug classes and their HIVDB annotations.
#     '''
#     mutations_df = pd.read_csv(mut_tsv, sep='\t')#we read the tsv file
#     mutations_filter_freq = mutations_df[mutations_df['Freq'] > lower_cutoff]#we apply the frequency cutoff
#     mutations_filter_freq = mutations_filter_freq[mutations_filter_freq['Ref'] != mutations_filter_freq['Mut']]##We filter out the rows where Ref == Mut
#     PR_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'PR')
#     RT_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'RT')
#     IN_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'IN')
#     INI_annot_df = HIVDB_singlemut_annot(list(IN_mut_freq.keys()), 'INI')
#     NNRTI_annot_df = HIVDB_singlemut_annot(list(RT_mut_freq.keys()), 'NNRTI')
#     NRTI_annot_df = HIVDB_singlemut_annot(list(RT_mut_freq.keys()), 'NRTI')
#     PI_annot_df = HIVDB_singlemut_annot(list(PR_mut_freq.keys()), 'PI')
#     INI_annot_df['Freq'] = IN_mut_freq.values()
#     INI_annot_df['Prot'] = 'IN'
#     NNRTI_annot_df['Freq'] = RT_mut_freq.values()
#     NNRTI_annot_df['Prot'] = 'RT'
#     NRTI_annot_df['Freq'] = RT_mut_freq.values()
#     NRTI_annot_df['Prot'] = 'RT'
#     PI_annot_df['Freq'] = PR_mut_freq.values()
#     PI_annot_df['Prot'] = 'PR'
    
#     merged_NNRTI_NRTI = pd.merge(NNRTI_annot_df, NRTI_annot_df, on=['Mutation', 'Freq', 'Prot'], how='outer', suffixes=('_NNRTI', '_NRTI'))##We unify NNRTI and NRTI dataframes
#     merged_NNRTI_NRTI['Annotation'] = merged_NNRTI_NRTI.apply(lambda row: row['Annotation_NNRTI'] if row['Annotation_NNRTI'] != 'Unknown' else row['Annotation_NRTI'], axis=1) ##The NRTI and NNRTI commented positions do not overlap
#     merged_NNRTI_NRTI['Comment'] = merged_NNRTI_NRTI.apply(lambda row: row['Comment_NNRTI'] if row['Comment_NNRTI'] != 'Unknown' else row['Comment_NRTI'], axis=1)
#     merged_NNRTI_NRTI = merged_NNRTI_NRTI[['Mutation', 'Freq', 'Prot', 'Annotation', 'Comment']]
#     merged_NNRTI_NRTI["Pos"] = merged_NNRTI_NRTI['Mutation'].apply(lambda x: int(re.sub(r'\D+', '', x)))
#     merged_NNRTI_NRTI = merged_NNRTI_NRTI.sort_values(by=['Pos'], ascending = True).drop(columns=['Pos'], axis=1)


#     mutations_df = pd.concat([INI_annot_df, merged_NNRTI_NRTI, PI_annot_df], axis=0, ignore_index=True)
#     mutations_df = mutations_df[['Mutation', 'Freq', 'Prot', 'Annotation', 'Comment']]

#     return mutations_df

# def extract_prot_mut_freq(mutations_df, prot = ''):
#     '''Outputs a list of dictionaries, one for each alternative list of mutations.
#        Each dictionary presents the mutations for a given protein and their observed frequency, coming from a mutation_freq.tsv file.
#     INPUT:
#     mutations_df: tsv file output from annotate_vcf.py with the observed mutations. Columns are Position, Ref, Mut, Freq and Prot.
#     prot: protein to filter the mutations. PR, RT or IN. If not one of those, it will return all mutations.
    
#     OUTPUT:
#     prot_mut_freq: dictionary with the mutations for the given proteinas keys and the frequencies as values.
#     '''
#     if prot in ['PR', 'RT', 'IN']:
#         mutations_df = mutations_df[mutations_df['Prot'] == prot]
    
#     prot_mut_freq = {} #we get a dictionary with the mutations as keys and the frequencies as values
#     for i, row in mutations_df.iterrows():
        
#         full_mut = mutations_df.loc[i, 'Ref'] + mutations_df.loc[i, 'Position'].astype(str) + mutations_df.loc[i, 'Mut']

#         if full_mut not in prot_mut_freq:
#             prot_mut_freq[full_mut] = row['Freq']
#         else:
#             prot_mut_freq[full_mut] += row['Freq']
    
#     return prot_mut_freq

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