###Labels a list of mutations as Resistant or Susceptible for different anti-HIV drugs
import pandas as pd
import re
import pickle
import typer
from typing_extensions import Annotated, Optional
from itertools import product
import matplotlib.pyplot as plt

##Drugs per drug class used in this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']

def HIVDB_table(mut_tsv, lower_cutoff: float = 0.015):
    '''Generates a table with the observed mutations and their HIVDB classifications for all drugs and drug class.
    INPUT:
    mut_tsv: tsv file output from annotate_vcf.py with the observed mutations. Columns are Position, Ref, Mut, Freq and Prot. 
    lower_cutoff: minimum frequency to be considered as a mutation.

    OUTPUT:
    mutations_df: dataframe with the mutations for all drugs and drug classes and their HIVDB annotations.
    '''
    mutations_df = pd.read_csv(mut_tsv, sep='\t')#we read the tsv file
    mutations_filter_freq = mutations_df[mutations_df['Freq'] > lower_cutoff]#we apply the frequency cutoff
    mutations_filter_freq = mutations_filter_freq[mutations_filter_freq['Ref'] != mutations_filter_freq['Mut']]##We filter out the rows where Ref == Mut
    PR_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'PR')
    RT_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'RT')
    IN_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'IN')
    INI_annot_df = HIVDB_singlemut_annot(list(IN_mut_freq.keys()), 'INI')
    NNRTI_annot_df = HIVDB_singlemut_annot(list(RT_mut_freq.keys()), 'NNRTI')
    NRTI_annot_df = HIVDB_singlemut_annot(list(RT_mut_freq.keys()), 'NRTI')
    PI_annot_df = HIVDB_singlemut_annot(list(PR_mut_freq.keys()), 'PI')
    INI_annot_df['Freq'] = IN_mut_freq.values()
    INI_annot_df['Prot'] = 'IN'
    NNRTI_annot_df['Freq'] = RT_mut_freq.values()
    NNRTI_annot_df['Prot'] = 'RT'
    NRTI_annot_df['Freq'] = RT_mut_freq.values()
    NRTI_annot_df['Prot'] = 'RT'
    PI_annot_df['Freq'] = PR_mut_freq.values()
    PI_annot_df['Prot'] = 'PR'
    
    merged_NNRTI_NRTI = pd.merge(NNRTI_annot_df, NRTI_annot_df, on=['Mutation', 'Freq', 'Prot'], how='outer', suffixes=('_NNRTI', '_NRTI'))##We unify NNRTI and NRTI dataframes
    merged_NNRTI_NRTI['Annotation'] = merged_NNRTI_NRTI.apply(lambda row: row['Annotation_NNRTI'] if row['Annotation_NNRTI'] != 'Unknown' else row['Annotation_NRTI'], axis=1) ##The NRTI and NNRTI commented positions do not overlap
    merged_NNRTI_NRTI['Comment'] = merged_NNRTI_NRTI.apply(lambda row: row['Comment_NNRTI'] if row['Comment_NNRTI'] != 'Unknown' else row['Comment_NRTI'], axis=1)
    merged_NNRTI_NRTI = merged_NNRTI_NRTI[['Mutation', 'Freq', 'Prot', 'Annotation', 'Comment']]
    merged_NNRTI_NRTI["Pos"] = merged_NNRTI_NRTI['Mutation'].apply(lambda x: int(re.sub(r'\D+', '', x)))
    merged_NNRTI_NRTI = merged_NNRTI_NRTI.sort_values(by=['Pos'], ascending = True).drop(columns=['Pos'], axis=1)


    mutations_df = pd.concat([INI_annot_df, merged_NNRTI_NRTI, PI_annot_df], axis=0, ignore_index=True)
    mutations_df = mutations_df[['Mutation', 'Freq', 'Prot', 'Annotation', 'Comment']]

    return mutations_df

def ensemble_table(mut_tsv, higher_cutoff: float = 0.15, HIVDB:bool = True, LSR: bool = True, RF: bool = True):
    '''Generates a table with the ensemble predictions for the mutation list above the high cutoff for all drugs and drug class.
    INPUT:
    mut_tsv: tsv file output from annotate_vcf.py with the observed mutations. Columns are Position, Ref, Mut, Freq and Prot.
    higher_cutoff: minimum frequency to be considered into the majority mutation list.

    OUTPUT:
    ensemble_preds_df: dataframe with the ensemble predictions for all drugs and drug classes.
    '''
    mutations_df = pd.read_csv(mut_tsv, sep='\t')#we read the tsv file
    mutations_filter_freq = mutations_df[mutations_df['Freq'] > higher_cutoff]#we apply the frequency cutoff
    mutations_filter_freq = mutations_filter_freq[mutations_filter_freq['Ref'] != mutations_filter_freq['Mut']]##We filter out the rows where Ref == Mut
    mutations_filter_freq = mutations_filter_freq[mutations_filter_freq['Mut'] != '*']

    ##We check if for a same Prot and Position there are different mutations
    IN_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'IN')
    RT_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'RT')
    PR_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'PR')
    IN_mut_freq_allcombs = get_mutlist_comb(list(IN_mut_freq.keys()))
    RT_mut_freq_allcombs = get_mutlist_comb(list(RT_mut_freq.keys()))
    PR_mut_freq_allcombs = get_mutlist_comb(list(PR_mut_freq.keys()))
    INI_mut_dic, NNRTI_mut_dic, NRTI_mut_dic, PI_mut_dic = {}, {}, {}, {}
    for mut_list in IN_mut_freq_allcombs:
        INI_ensemble_preds_df = ensemble_predictions(mut_list, 'INI', HIVDB = HIVDB, LSR = LSR, RF = RF)
        INI_ensemble_preds_binary = INI_ensemble_preds_df.drop(columns = list(set(['HIVDB_five_labels', 'LSR_RF'])&set(INI_ensemble_preds_df.columns)))
        INI_binary_dic = [df.drop(columns = list(set(['HIVDB_five_labels', 'LSR_RF'])&set(INI_ensemble_preds_df.columns))) for df in INI_mut_dic.values()]
        INI_is_df, INI_df_idx = df_in_list(INI_ensemble_preds_binary, INI_binary_dic)
        if not INI_is_df:
            INI_mut_dic[', '.join(mut_list)] = INI_ensemble_preds_df
        elif calc_HIVDB_score(mut_list, 'INI') > calc_HIVDB_score(list(INI_mut_dic.keys())[INI_df_idx], 'INI'): ##When resistance for all drugs is the same, we keep the mutations with the highest HIVDB score
            INI_mut_dic[', '.join(mut_list)] = INI_ensemble_preds_df
            INI_mut_dic.pop(list(INI_mut_dic.keys())[INI_df_idx])
        elif calc_HIVDB_score(mut_list, 'INI') == calc_HIVDB_score(list(INI_mut_dic.keys())[INI_df_idx], 'INI'):          
            new_mut_freq = 0 ##we get the summed frequency of the mutations
            for mut in mut_list:
                new_mut_freq += IN_mut_freq[mut]
            old_mut_freq = 0
            for mut in list(INI_mut_dic.keys())[INI_df_idx].split(', '):
                old_mut_freq += IN_mut_freq[mut]
            if new_mut_freq > old_mut_freq: ##we keep the most frequent of them
                INI_mut_dic[', '.join(mut_list)] = INI_ensemble_preds_df
                INI_mut_dic.pop(list(INI_mut_dic.keys())[INI_df_idx])

    
    for mut_list in RT_mut_freq_allcombs:
        NRTI_ensemble_preds_df = ensemble_predictions(mut_list, 'NRTI', HIVDB = HIVDB, LSR = LSR, RF = RF)
        NRTI_ensemble_preds_binary = NRTI_ensemble_preds_df.drop(columns = list(set(['HIVDB_five_labels', 'LSR_RF'])&set(NRTI_ensemble_preds_df.columns)))
        NNRTI_ensemble_preds_df = ensemble_predictions(mut_list, 'NNRTI', HIVDB = HIVDB, LSR = LSR, RF = RF)
        NNRTI_ensemble_preds_binary = NNRTI_ensemble_preds_df.drop(columns = list(set(['HIVDB_five_labels', 'LSR_RF'])&set(NNRTI_ensemble_preds_df.columns)))
        NRTI_binary_dic = [df.drop(columns = list(set(['HIVDB_five_labels', 'LSR_RF'])&set(NRTI_ensemble_preds_df.columns))) for df in NRTI_mut_dic.values()]
        NNRTI_binary_dic = [df.drop(columns = list(set(['HIVDB_five_labels', 'LSR_RF'])&set(NNRTI_ensemble_preds_df.columns))) for df in NNRTI_mut_dic.values()]
        is_NRTI_df, NRTI_df_idx = df_in_list(NRTI_ensemble_preds_binary, NRTI_binary_dic)
        is_NNRTI_df, NNRTI_df_idx = df_in_list(NNRTI_ensemble_preds_binary, NNRTI_binary_dic)

        if not is_NRTI_df:
            NRTI_mut_dic[', '.join(mut_list)] = NRTI_ensemble_preds_df
        elif calc_HIVDB_score(mut_list, 'NRTI') > calc_HIVDB_score(list(NRTI_mut_dic.keys())[NRTI_df_idx], 'NRTI'):
            NRTI_mut_dic[', '.join(mut_list)] = NRTI_ensemble_preds_df
            NRTI_mut_dic.pop(list(NRTI_mut_dic.keys())[NRTI_df_idx])
        elif calc_HIVDB_score(mut_list, 'NRTI') == calc_HIVDB_score(list(NRTI_mut_dic.keys())[NRTI_df_idx], 'NRTI'):
            new_mut_freq = 0
            for mut in mut_list:
                new_mut_freq += RT_mut_freq[mut]
            old_mut_freq = 0
            for mut in list(NRTI_mut_dic.keys())[NRTI_df_idx].split(', '):
                old_mut_freq += RT_mut_freq[mut]
            if new_mut_freq > old_mut_freq:
                NRTI_mut_dic[', '.join(mut_list)] = NRTI_ensemble_preds_df
                NRTI_mut_dic.pop(list(NRTI_mut_dic.keys())[NRTI_df_idx])

        if not is_NNRTI_df:
            NNRTI_mut_dic[', '.join(mut_list)] = NNRTI_ensemble_preds_df
        elif calc_HIVDB_score(mut_list, 'NNRTI') > calc_HIVDB_score(list(NNRTI_mut_dic.keys())[NNRTI_df_idx], 'NNRTI'):
            NNRTI_mut_dic[', '.join(mut_list)] = NNRTI_ensemble_preds_df
            NNRTI_mut_dic.pop(list(NNRTI_mut_dic.keys())[NNRTI_df_idx])
        elif calc_HIVDB_score(mut_list, 'NNRTI') == calc_HIVDB_score(list(NNRTI_mut_dic.keys())[NNRTI_df_idx], 'NNRTI'):
            new_mut_freq = 0
            for mut in mut_list:
                new_mut_freq += RT_mut_freq[mut]
            old_mut_freq = 0
            for mut in list(NNRTI_mut_dic.keys())[NNRTI_df_idx].split(', '):
                old_mut_freq += RT_mut_freq[mut]
            if new_mut_freq > old_mut_freq:
                NNRTI_mut_dic[', '.join(mut_list)] = NNRTI_ensemble_preds_df
                NNRTI_mut_dic.pop(list(NNRTI_mut_dic.keys())[NNRTI_df_idx])
        
    
    for mut_list in PR_mut_freq_allcombs:
        PI_ensemble_preds_df = ensemble_predictions(mut_list, 'PI', HIVDB = HIVDB, LSR = LSR, RF = RF)
        PI_ensemble_preds_binary = PI_ensemble_preds_df.drop(['HIVDB_five_labels', 'LSR_RF'], errors='ignore')
        PI_binary_dic = [df.drop(['HIVDB_five_labels', 'LSR_RF'], errors='ignore') for df in PI_mut_dic.values()]
        PI_is_df, PI_df_idx = df_in_list(PI_ensemble_preds_binary, PI_binary_dic)
        if not PI_is_df:
            PI_mut_dic[', '.join(mut_list)] = PI_ensemble_preds_df
        elif calc_HIVDB_score(mut_list, 'PI') > calc_HIVDB_score(list(PI_mut_dic.keys())[PI_df_idx], 'PI'):
            PI_mut_dic.pop(list(PI_mut_dic.keys())[PI_df_idx])
            PI_mut_dic[', '.join(mut_list)] = PI_ensemble_preds_df
        elif calc_HIVDB_score(mut_list, 'PI') == calc_HIVDB_score(list(PI_mut_dic.keys())[PI_df_idx], 'PI'):
            new_mut_freq = 0
            for mut in mut_list:
                new_mut_freq += PR_mut_freq[mut]
            old_mut_freq = 0
            for mut in list(PI_mut_dic.keys())[PI_df_idx].split(', '):
                old_mut_freq += PR_mut_freq[mut]
            if new_mut_freq > old_mut_freq:
                PI_mut_dic.pop(list(PI_mut_dic.keys())[PI_df_idx])
                PI_mut_dic[', '.join(mut_list)] = PI_ensemble_preds_df

    return [INI_mut_dic, NNRTI_mut_dic, NRTI_mut_dic, PI_mut_dic]

def extract_prot_mut_freq(mutations_df, prot = ''):
    '''Outputs a list of dictionaries, one for each alternative list of mutations.
       Each dictionary presents the mutations for a given protein and their observed frequency, coming from a mutation_freq.tsv file.
    INPUT:
    mutations_df: tsv file output from annotate_vcf.py with the observed mutations. Columns are Position, Ref, Mut, Freq and Prot.
    prot: protein to filter the mutations. PR, RT or IN. If not one of those, it will return all mutations.
    
    OUTPUT:
    prot_mut_freq: dictionary with the mutations for the given proteinas keys and the frequencies as values.
    '''
    if prot in ['PR', 'RT', 'IN']:
        mutations_df = mutations_df[mutations_df['Prot'] == prot]
    
    prot_mut_freq = {} #we get a dictionary with the mutations as keys and the frequencies as values
    for i, row in mutations_df.iterrows():
        
        full_mut = mutations_df.loc[i, 'Ref'] + mutations_df.loc[i, 'Position'].astype(str) + mutations_df.loc[i, 'Mut']

        if full_mut not in prot_mut_freq:
            prot_mut_freq[full_mut] = row['Freq']
        else:
            prot_mut_freq[full_mut] += row['Freq']
    
    return prot_mut_freq

def get_mutlist_comb(mut_list: list):
    ''' Generates a list of all possible combinations of mutations in the input list.
    INPUT:
    mut_list: list of mutations i.e [M41L,L90M].
    
    OUTPUT:
    comb_list: list of all possible combinations of mutations. List of lists.
    '''
    pos_list = [re.sub(r'\D+', '', mut) for mut in mut_list] #we get the positions of the mutations
    
    pos_dict = {}    #we make a dictionary with the positions as keys and a list of their index as values
    for i, pos in enumerate(pos_list):
        if pos not in pos_dict:
            pos_dict[pos] = [i]
        else:
            pos_dict[pos].append(i)

    comb_list = []    #we create a list of all possible combinations of mutations
    for i in range(len(pos_dict)):
        pos = list(pos_dict.keys())[i]
        if len(pos_dict[pos]) > 1:
            comb_list.append([mut_list[j] for j in pos_dict[pos]])
        else:
            comb_list.append([mut_list[pos_dict[pos][0]]])
    
    comb_list = list(product(*comb_list))#we create all the combinations of the mutations
    comb_list = [list(comb) for comb in comb_list]

    return comb_list

def unify_mut_list(list_of_mut_list: list):
    '''Unifies a list of lists of mutations into a single list of mutations.
    INPUT:
    list_of_mut_list: list of lists of mutations i.e [[M41L,L90M], [M41K,L90M]].
    
    OUTPUT:
    unified_list: list of unique mutations. i.e [M41LK, L90M].
    '''
    mut_dic = {}
    for mut_list in list_of_mut_list:
        for mut in mut_list:
            pos = re.sub(r'\D+', '', mut)
            aa_ref = re.sub(r'\d+', '', mut)[0]
            aa_mut = re.sub(r'\d+', '', mut)[1:]
            if pos not in mut_dic.keys():
                mut_dic[pos] = aa_ref + aa_mut
            elif mut_dic[pos][1] != aa_mut:
                mut_dic[pos] = mut_dic[pos] + aa_mut

    unified_list = []
    for pos in mut_dic.keys():
        unified_list.append(mut_dic[pos][0] + pos + ''.join(mut_dic[pos][1:]))
    
    return unified_list

def calc_HIVDB_score(mut_list: list, dataset: str):
    '''Calculates the HIVDB score for a list of mutations.
    INPUT:
    mutation_list: list of mutations i.e [M41L,L90M].
    dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.

    OUTPUT:
    score: HIVDB score for the mutations.
    '''
    if dataset == "INI":
        drug_class = INIs
    elif dataset == "NNRTI":
        drug_class = NNRTIs
    elif dataset == "NRTI":
        drug_class = NRTIs
    elif dataset == "PI":
        drug_class = PIs

    mutation_scores = pd.read_csv(f'HIVDB_rules/{dataset}_muts_score_Stanford_HIVDB', sep=',') #we read the mutation scores
    combination_scores = pd.read_csv(f'HIVDB_rules/{dataset}_combinations_score_Stanford_HIVDB', sep=',') #we read the combination scores
        
    mutation_scores.columns = mutation_scores.columns.str.replace('/r','')
    combination_scores.columns = combination_scores.columns.str.replace('/r','')    
    
    HIVDB_score_total = 0
    for drug in drug_class:
        
        drug_score = 0
            
        ###We add the single mutation scores
        for mut in mut_list:
            if mut in mutation_scores['Rule'].values and drug in mutation_scores.columns:
                drug_score += float(mutation_scores[mutation_scores['Rule']==mut][f'{drug}'].iloc[0]) #we add the score

        ##We add the combination scores
        for combination_rule in combination_scores['Combination Rule'].values:
            combination_rule = combination_rule.split(' + ')
            all_in = True
            for mut_c_rule in combination_rule:
                #we separate the mutation from the position
                position = re.sub(r'\D+', '', mut_c_rule)
                mut_aa = re.sub(r'\d+', '', mut_c_rule)
                mut_aa_1 = mut_aa[0]
                mut_aa_2 = mut_aa[1:]
                if len(mut_aa_2) == 1:
                    if mut_c_rule not in mut_list:
                        all_in = False
                        break
                else:
                    posib_mut = False
                    for l in mut_aa_2:
                        if mut_aa_1+position+l in mut_list:
                            posib_mut = True
                            break
                        else:
                            posib_mut = False

                    if not posib_mut:
                        all_in = False
                        break
                
            if all_in: #if all mutations in the combination rule are in the input list
                if drug in combination_scores.columns:
                    drug_score += float(combination_scores[combination_scores['Combination Rule']==' + '.join(combination_rule)][f'{drug}'].iloc[0])
            
        HIVDB_score_total += drug_score ##We add it to the total score

    return HIVDB_score_total


def df_in_list(df, list_of_dfs):
    '''Checks if a dataframe is in a list of dataframes.
    Outputs the index of the dataframe in the list if it is found.'''
    df = df.sort_index(axis=1).reset_index(drop=True)
    for i, other_df in enumerate(list_of_dfs):
        other_df = other_df.sort_index(axis=1).reset_index(drop=True)
        if df.equals(other_df):
            return True, i
    return False, None

def HIVDB_singlemut_annot(mut_list, dataset: str):
    '''Retrieves the HIVDB annotation for each mutation in a list of mutations.
    INPUT:
    mut_list: list of mutations i.e [M41L,L90M].
    dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.

    OUTPUT:
    mut_annot_df: dataframe with the mutations and their HIVDB annotations.
    '''
    #we read the comments files for the dataset
    comments_file = pd.read_csv(f'HIVDB_rules/{dataset}_comments_Stanford_HIVDB', sep=',')
    comments_file.columns = comments_file.columns.str.replace('\n','')

    mut_list = [mut[1:] for mut in mut_list] #we remove the first letter of the mutations
    mut_annot = []
    for mut in mut_list:
        position = re.sub(r'\D+', '', mut)
        mut_aa = re.sub(r'\d+', '', mut)
        comment_in = False
        for i, comm in enumerate(comments_file['Condition'].values):
            if re.match(r'^\d+[A-Z]+$', comm) and position==re.sub(r'\D+', '', comm): #we check if the mutation is in the comment
                comm_ambig = re.sub(r'\d+', '', comm)
                if mut_aa in comm_ambig:
                    mut_annot.append([mut, 
                                     comments_file['Comment/Mutation Type'].values[i], 
                                     comments_file['Comment'].values[i]])
                    comment_in = True
                    break
        
        if not comment_in:
            mut_annot.append([mut, "Unknown", "Unknown"])
                       
    mut_annot_df = pd.DataFrame(mut_annot, columns=['Mutation', 'Annotation', 'Comment'])
    return mut_annot_df


###HIVDB offline implementation for the prediction of the resistance
def HIVDB_pred(mut_list, dataset: str, score = "SR"):
    '''Calculates the HIVDB score given a set of mutations using the HIVDB scores for mutations and combinations extracted from StanfordHIVDB HIVDB program.
        Scores downloaded on 10th April 2025.
        
        INPUT:
        mut_list: list of mutations separated by ', ' as string i.e M41L, L90M
        dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.
        score: binary or HIVDB labeling system. 'SR' for binary Susceptible/Resistant labeling or 'HIVDB' for the HIVDB labeling system (5 labels).

        OUTPUT:
        resistance_pred_df: dataframe with the predicted resistance labels for each drug in the drug class.
    '''

    if dataset == "INI":
        drug_class = INIs
    elif dataset == "NNRTI":
        drug_class = NNRTIs
    elif dataset == "NRTI":
        drug_class = NRTIs
    elif dataset == "PI":
        drug_class = PIs

    mutation_scores = pd.read_csv(f'HIVDB_rules/{dataset}_muts_score_Stanford_HIVDB', sep=',') #we read the mutation scores
    combination_scores = pd.read_csv(f'HIVDB_rules/{dataset}_combinations_score_Stanford_HIVDB', sep=',') #we read the combination scores
        
    mutation_scores.columns = mutation_scores.columns.str.replace('/r','')
    combination_scores.columns = combination_scores.columns.str.replace('/r','')    
    
    resistance_pred_df = pd.DataFrame(columns=drug_class) #we create the dataframe to store the results
    for drug in drug_class:
        
        drug_score = 0
            
        ###We add the single mutation scores
        for mut in mut_list:
            if mut in mutation_scores['Rule'].values and drug in mutation_scores.columns:
                drug_score += float(mutation_scores[mutation_scores['Rule']==mut][f'{drug}'].iloc[0]) #we add the score

        ##We add the combination scores
        for combination_rule in combination_scores['Combination Rule'].values:
            combination_rule = combination_rule.split(' + ')
            all_in = True
            for mut_c_rule in combination_rule:
                #we separate the mutation from the position
                position = re.sub(r'\D+', '', mut_c_rule)
                mut_aa = re.sub(r'\d+', '', mut_c_rule)
                mut_aa_1 = mut_aa[0]
                mut_aa_2 = mut_aa[1:]
                if len(mut_aa_2) == 1:
                    if mut_c_rule not in mut_list:
                        all_in = False
                        break
                else:
                    posib_mut = False
                    for l in mut_aa_2:
                        if mut_aa_1+position+l in mut_list:
                            posib_mut = True
                            break
                        else:
                            posib_mut = False

                    if not posib_mut:
                        all_in = False
                        break
                
            if all_in: #if all mutations in the combination rule are in the input list
                if drug in combination_scores.columns:
                    drug_score += float(combination_scores[combination_scores['Combination Rule']==' + '.join(combination_rule)][f'{drug}'].iloc[0])
            
        resistance_pred_df.loc[0, drug] = drug_score ##We add it to the dataframe


    ###We label the resistance based on the HIVDB scoring method
    if score == "SR": #we label them just as Susceptible or Resistant
        ##Resistant if score > 29, Susceptible from 0 to 29
        for drug in resistance_pred_df.columns:
            if resistance_pred_df.loc[0, drug] > 29:
                resistance_pred_df.loc[0, drug] = "Resistant"
            else:
                resistance_pred_df.loc[0, drug] = "Susceptible"
    elif score == "HIVDB": #we label them following the original HIVDB labeling
        for drug in resistance_pred_df.columns:
            if resistance_pred_df.loc[0, drug] > 59:
                resistance_pred_df.loc[0, drug] = "High-Level Resistance"
            elif resistance_pred_df.loc[0, drug] > 29:
                resistance_pred_df.loc[0, drug] = "Intermediate Resistance"
            elif resistance_pred_df.loc[0, drug] > 14:
                resistance_pred_df.loc[0, drug] = "Low-Level Resistance"
            elif resistance_pred_df.loc[0, drug] > 9:
                resistance_pred_df.loc[0, drug] = "Potential Low-Level Resistance"
            else:
                resistance_pred_df.loc[0, drug] = "Susceptible"
    
    return resistance_pred_df

###Linear regression implementation
def LSR_res_pred(mut_list, dataset: str, report_RF = False):
    '''Calculates the resistance factor given a set of mutations using the coefficients obtained from LSR with TSM features.
        INPUT:
        mut_list: list of mutations i.e [M41L,L90M]. 
        dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.
        report_RF: if True, the resistance factor will be reported

        OUTPUT:
        resistance_pred_df: dataframe with the predicted resistance labels for each drug in the drug class.       
    '''
    cutoff = {'FPV':3, 'ATV':3, 'IDV':3, 'LPV':9, 'NFV':3, 'SQV':3, 'TPV':2, 'DRV':10, #resistance fold cutoff defined to decide drug resistance
              '3TC':3, 'ABC':2, 'AZT':3, 'D4T':1.5, 'DDI':1.5, 'TDF':1.5, 
              'EFV':3, 'ETR':3, 'NVP':3, 'RPV':3,
              'RAL': 4, 'EVG': 4, 'DTG': 3.5, 'BIC':3.5}


    if dataset == "INI":
        drug_class = INIs
    elif dataset == "NNRTI":
        drug_class = NNRTIs
    elif dataset == "NRTI":
        drug_class = NRTIs
    elif dataset == "PI":
        drug_class = PIs
    
    mut_list = [mut[1:] for mut in mut_list] 
    
    resistance_pred_df = pd.DataFrame(columns=drug_class) #we create the dataframe to store the results

    for drug in drug_class:
        #we read the coefficients
        LSR_coefficients = pd.read_csv(f"linear_regression_coefficients/OLS_{drug}_combinations_tsm_all_folds.txt")
        coefficients = {}

        for n in range(LSR_coefficients.shape[0]):
            if n == 0:
                intercept = LSR_coefficients.iloc[n,0].split(' ')[1]
            else:
                coef_line = LSR_coefficients.iloc[n,0].split(' ')
                mutation = coef_line[0][2:]
                coef = coef_line[1]
                coefficients[mutation] = coef

        RF = 0 + float(intercept) #resistance factor

        ##we multiply by 1 the coefficients of the mutations that are present in the mut_list
        if len(mut_list) > 0:
            for mut in mut_list:
                if mut in coefficients:
                    RF += float(coefficients[mut])        
        
            ##and now we parse the mutation combinations
            ###we look for the combinations in the coefficient keys
            for mut_key in coefficients.keys():
                if "." in mut_key: #we look for the combinations coefficients
                    comb_muts = mut_key.split('.')

                    #now we will go through comb_muts and check if all the mutations are in mut_list
                    all_in = True
                    for m in comb_muts:
                        ##we separate the mutation from the position
                        letters = re.sub(r'\d+', '', m)
                        position = re.sub(r'\D+', '', m)

                        if len(letters) == 1:
                            if m not in mut_list:
                                all_in = False
                                break
                        else:
                            posib_mut = False
                            for l in letters:
                                if position+l in mut_list:
                                    posib_mut = True
                                    break
                                else:
                                    posib_mut = False

                            if not posib_mut:
                                all_in = False
                                break
                    
                    if all_in:
                        RF += float(coefficients[mut_key])
    
    
        ##we elevate the RF to the power of 10, as it is in log scale
        RF = 10**RF
        ##we apply the cutoff and predict the resistance
        if RF > float(cutoff[drug]):
            label = "Resistant"
        else:
            label = "Susceptible"

        if report_RF: #if report_RF is True, we report the predicted resistance factor
            label = RF
        
        resistance_pred_df.loc[0, drug] = label
    
    return resistance_pred_df

###Random Forest implementation
def RandomForest_HIV(mut_list, dataset:str):
    '''The function performs anti-HIV drug resistance prediction with a Random Forest classification based on Raposo et al. 2020 implementation.
    INPUT:
    mut_list: list of mutations i.e [M41L,L90M]. 
    dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.

    OUTPUT:
    resistance_pred_df: dataframe with the predicted resistance labels for each drug in the drug class.
    '''

    if dataset == "INI":
        drug_class = INIs
    elif dataset == "NNRTI":
        drug_class = NNRTIs
    elif dataset == "NRTI":
        drug_class = NRTIs
    elif dataset == "PI":
        drug_class = PIs

    mut_list = [mut[1:] for mut in mut_list] 

    resistance_df = pd.DataFrame(columns=drug_class) #we create the dataframe to store the results

    if dataset == "PI":
        ref_seq = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    elif dataset == "NNRTI":
        ref_seq = "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWT"
    elif dataset == "NRTI":
        ref_seq = "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWT"
    elif dataset == "INI":
        ref_seq = "FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"

    X = pd.DataFrame([list(ref_seq)], columns=[f'V{i+1}' for i in range(len(ref_seq))]) #we create a dataframe with a column for each position in the sequence as input
        
    for mut in mut_list: 
        pos = int(mut[:-1])
        aa = mut[-1]
        X.loc[0, f'V{pos}'] = aa

    #we apply the Kyte-Doolittle hydrophobicity scale to the aminoacids in the X dataframe
    #we first define the scale
    kyte_doolittle_scale = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
        'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
        'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
        'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }

    #we apply the scale to the X dataframe
    for col in X.columns:
        X[col] = X[col].apply(lambda x: kyte_doolittle_scale[x])

    for drug in drug_class:
        #we load the .pkl model
        with open(f"random_forest_models/random_forest_python_{drug}_RF_model_allfolds.pkl", "rb") as file:
            rf = pickle.load(file)
    
        #we only keep the columns that are in the model
        X_drug = X[rf.feature_names_in_] #we only keep the columns that are in the model

        #we predict the resistance
        predictions = rf.predict(X_drug)

        #we save it into the dataframe
        resistance_df.loc[0, drug] = predictions[0]

    return resistance_df

###Allow prediction of only one of the methods with options
###Also get the resistance factor for the LSR method

###Ensemble prediction
def ensemble_predictions(mut_list, dataset: str, HIVDB:bool = True, LSR: bool = True, RF: bool = True):
    '''Generates an ensemble prediction by majority voting of the HIVDB, Linear Regression and Random Forest predictions.

    INPUT:
    mut_list: list of mutations i.e [M41L,L90M].
    dataset: drug class. 'INI', 'NNRTI', 'NRTI', 'PI'.
    HIVDB: if True, the HIVDB prediction will be used.
    LSR: if True, the Linear Regression prediction will be used.
    RF: if True, the Random Forest prediction will be used.
    
    OUTPUT:
    resistance_df: dataframe with the predicted resistance labels for each drug in the drug class.
    '''
    if not HIVDB and not LSR and not RF:
        print("At least one method should be selected")
        return
    
    ##we predict the resistance with the different methods
    available_preds, available_index = [], []
    ##HIVDB
    if HIVDB:
        HIVDB_results = HIVDB_pred(mut_list, dataset)
        HIVDB_results_five_labels = HIVDB_pred(mut_list, dataset, score = "HIVDB") #we get the five labels for the HIVDB method
        available_preds.append(HIVDB_results), available_preds.append(HIVDB_results_five_labels)
        available_index.append('HIVDB'), available_index.append('HIVDB_five_labels')
    ##Linear Regression
    if LSR:
        LSR_results = LSR_res_pred(mut_list, dataset)
        LSR_RF = LSR_res_pred(mut_list, dataset, report_RF = True) #we get the resistance factor for the LSR method
        available_preds.append(LSR_results), available_preds.append(LSR_RF)
        available_index.append('LSR'), available_index.append('LSR_RF')
    ##Random Forest
    if RF:
        RF_results = RandomForest_HIV(mut_list, dataset)
        available_preds.append(RF_results)
        available_index.append('RF')

    ##We concatenate the available dataframes
    ensemble_results = pd.concat(available_preds, axis=0, ignore_index=True)
    ensemble_results.index = available_index

    if HIVDB and LSR and RF:
        ensemble_results.loc['Ensemble'] = ensemble_results.drop(['HIVDB_five_labels', 'LSR_RF']).mode().iloc[0] #we get the mode of the predictions

    ensemble_results = ensemble_results.T
    return ensemble_results

def check_mut_input(input_mut: str):
    '''Checks if the input mutations are valid. If not, it returns an Error message.
    
    INPUT:
    input_mut: input string with the list of mutations i.e "M41L, L90M".
    
    OUTPUT:
    valid_mut: list of valid mutations or empty list.
    '''
    # Check if the input is a string
    if not isinstance(input_mut, str):
        print("Input mutations should be a string")
        return

    # Check if the input is empty
    if pd.isnull(input_mut):
        print("No mutations provided")
        return []
    elif len(input_mut) == 0:
        print("No mutations provided")
        return []
    
    # Check if the input is correctly formatted
    mutations = input_mut.split(',')
    mutations = [mut.strip() for mut in mutations]  # Remove leading/trailing spaces
   
    for mut in mutations: #we check if the mutations are in the required format [A-Z][0-9]+[A-Z]
        if not re.match(r'^[A-Z]\d+[A-Z]$', mut):
            print(f"Invalid mutation format: {mut}. Please input a list of comma separated mutations in the correct format. i.e. M41L,L90M")
            return
    
    return mutations

def write_HIVDB_table(HIVDB_table, cutoff, handle, unknown: bool = False, comments: bool = False):
    '''Writes the ensemble results to a markdown file.
    INPUT:
    HIVDB_table: Dataframe with the single mutations, frequency and HIVDB comments. Output from HIVDB_table function.
    cutoff: frequency cutoff to consider the mutations.
    handle: file handle to write the results.
    unknown: if True, the unknown annotations will be written.
    comments: if True, the comments will be written.
    '''

    handle.write("Single mutation annotations obtained from HIVDB program.")
    handle.write("Mutations below " + str(cutoff) + " frequency were not included in the analysis.\n")
    IN_table, RT_table, PR_table = HIVDB_table[HIVDB_table["Prot"]=="IN"], HIVDB_table[HIVDB_table["Prot"]=="RT"], HIVDB_table[HIVDB_table["Prot"]=="PR"]
    if not unknown:
        IN_table, RT_table, PR_table = IN_table[IN_table["Annotation"] != "Unknown"], RT_table[RT_table["Annotation"] != "Unknown"], PR_table[PR_table["Annotation"] != "Unknown"]

    if IN_table.shape[0] > 0:
        handle.write("**INTEGRASE** mutations:\n")
        handle.write("|{: ^14}|{: ^11}|{: ^12}|\n".format('Mutation','Frequency', 'Annotation'))
        handle.write('|:{:-^12}:|:{:-^9}:|:{:-^10}:|'.format('', '', ''))
        handle.write('\n')
        for n in range(IN_table.shape[0]):    
            handle.write('|{: ^14}|{: ^11}|{: ^12}|'.format(IN_table["Mutation"].iloc[n], IN_table["Freq"].iloc[n], IN_table["Annotation"].iloc[n]))
            handle.write("\n")
        
        if comments:
            handle.write("\nComments:\n\n")
            for n in range(IN_table.shape[0]):
                if IN_table["Comment"].iloc[n] != "Unknown":
                    handle.write(f"- **{IN_table['Mutation'].iloc[n]}**: {IN_table['Comment'].iloc[n]}\n") 
    else:
        handle.write("\nNo **INTEGRASE** mutations found.\n")
    
    
    if RT_table.shape[0] > 0:
        handle.write("\n**REVERSE TRANSCRIPTASE** mutations:\n\n")
        handle.write("|{: ^14}|{: ^11}|{: ^12}|\n".format('Mutation','Frequency', 'Annotation'))
        handle.write('|:{:-^12}:|:{:-^9}:|:{:-^10}:|'.format('', '', ''))
        handle.write('\n')
        for n in range(RT_table.shape[0]):
            handle.write('|{: ^14}|{: ^11}|{: ^12}|'.format(RT_table["Mutation"].iloc[n], RT_table["Freq"].iloc[n], RT_table["Annotation"].iloc[n]))
            handle.write("\n")
        
        if comments:
            handle.write("\nComments:\n\n")
            for n in range(RT_table.shape[0]):
                if RT_table["Comment"].iloc[n] != "Unknown":
                    handle.write(f"- **{RT_table['Mutation'].iloc[n]}**: {RT_table['Comment'].iloc[n]}\n")    
    else:
        handle.write("\nNo **REVERSE TRANSCRIPTASE** mutations found.\n")



    if PR_table.shape[0] > 0:
        handle.write("\n**PROTEASE** mutations:\n\n")
        handle.write("|{: ^14}|{: ^11}|{: ^12}|\n".format('Mutation','Frequency', 'Annotation'))
        handle.write('|:{:-^12}:|:{:-^9}:|:{:-^10}:|'.format('', '', ''))
        handle.write('\n')
        for n in range(PR_table.shape[0]):
            handle.write('|{: ^14}|{: ^11}|{: ^12}|'.format(PR_table["Mutation"].iloc[n], PR_table["Freq"].iloc[n], PR_table["Annotation"].iloc[n]))
            handle.write("\n")
        
        if comments:
            handle.write("\nComments:\n\n")
            for n in range(PR_table.shape[0]):
                if PR_table["Comment"].iloc[n] != "Unknown":
                    handle.write(f"- **{PR_table['Mutation'].iloc[n]}**: {PR_table['Comment'].iloc[n]}\n")    
    else:
        handle.write("\nNo **PROTEASE** mutations found.\n")  
    


def write_ensemble_table(ensemble_results, cutoff, handle, HIVDB = True, LSR = True, RF = True):
    '''Writes the ensemble results to a markdown file.
    INPUT:
    ensemble_results: List of dictionaries for each drug class containing ensemble predictions for each mutation list combination. Output from ensemble_predictions function.
    cutoff: cutoff for the mutations.
    handle: file handle to write the results.
    '''
    handle.write("Mutations below " + str(cutoff) + " frequency were not included in the analysis.\n\n")
    
    for i, drug_class_dic in enumerate(ensemble_results):

        if '' in drug_class_dic.keys(): #we skip the empty keys
            continue

        if i == 0:
            handle.write("**INI** predictions:\n")
            dataset = "INI"
        elif i == 1:
            handle.write("\n**NNRTI** predictions:\n")
            dataset = "NNRTI"
        elif i == 2:
            handle.write("\n**NRTI** predictions:\n")
            dataset = "NRTI"
        elif i == 3:
            handle.write("\n**PI** predictions:\n")
            dataset = "PI"
        
        for mut_comb in drug_class_dic.keys():

            handle.write("\n**" + mut_comb.replace(", ", "+") + "**\n\n")

            mut_df = drug_class_dic[mut_comb]
            mut_df["Drug"] = mut_df.index
            print(mut_df)

            keep_cols, keep_color = ["Drug"], ["Drug"]
            if HIVDB:
                keep_cols.append("HIVDB_five_labels")
                keep_color.append("HIVDB")
            if LSR:
                mut_df["LSR_RF"] = mut_df["LSR_RF"].astype(float).round(4)
                keep_cols.append("LSR_RF")
                keep_color.append("LSR")
            if RF:
                keep_cols.append("RF")
                keep_color.append("RF")
            if HIVDB and LSR and RF:
                keep_cols.append("Ensemble")
                keep_color.append("Ensemble")

            mut_df_plot = mut_df[keep_cols]
            mut_df_color = mut_df[keep_color]
            # columns = mut_df_plot.columns
            mut_df_plot.columns = [col.replace("_five_labels", "").replace("LSR_RF", "Linear Regression").replace("RF", "Random Forest") for col in list(mut_df_plot.columns)]
            
            fig, ax = plt.subplots(figsize=(10, mut_df_plot.shape[0] * 0.25))
            ax.axis('off')
            
            table = ax.table(cellText=mut_df_plot.values, colLabels=mut_df_plot.columns, cellLoc='center', loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(12)
            table.scale(1.4, 1.4)
            table.auto_set_column_width([0,1])

            for i in range(mut_df_plot.shape[1]):
                table[(0, i)].set_fontsize(13)
                table[(0, i)].set_text_props(weight="bold")

            for i in range(mut_df_plot.shape[0]):
                table[(i+1, 0)].set_text_props(weight="bold")
                for j in range(mut_df_plot.shape[1]):
                    if mut_df_color.iloc[i, j] == "Resistant":
                        table[(i+1, j)].set_facecolor("red")
                        table[(i+1, j)].set_alpha(0.5)
                    elif mut_df_color.iloc[i, j] == "Susceptible":
                        table[(i+1, j)].set_facecolor("green")
                        table[(i+1, j)].set_alpha(0.3)
                    else:
                        table[(i+1, j)].set_facecolor("white")

            s=f"example_files/{dataset}_{mut_comb.replace(', ', '_')}.pdf"
            plt.tight_layout()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            plt.savefig(s, bbox_inches='tight', pad_inches = 0.05, transparent=True)
            plt.close(fig)
            handle.write(r"\begin{center}")
            handle.write("\n")
            handle.write(r"\includegraphics[width=\textwidth]")
            handle.write("{%s}\n" % s)
            handle.write("\end{center}\n")


def write_coverage_disclaimer(coverage_tsv, handle, read_cutoff:int = 1):
    '''Writes the coverage disclaimer and writes a step plot showing the position of the sample's coverage.
    INPUT:
    coverage_tsv: tsv file with the coverage data, annotated with annotate_vcf.py. Similar to samtools depth.
    read_cutoff: minimum number of reads to consider a position as covered.
    handle: file handle to write the results.
    '''
    handle.write("\n## Coverage disclaimer\n\n")
    handle.write(f"The coverage of the sample was calculated from the FASTQ files. Only Drug Resistance Mutation (DRM) positions are analysed for robustness assessment. Positions were considered as covered above {read_cutoff} read(s).\n")

    coverage_data = pd.read_csv(coverage_tsv, sep='\t')
    coverage_data.columns = ['Reference','Position', 'Reads']#, 'Protein']

    for dataset in ["INI", "NNRTI", "NRTI", "PI"]:
        if dataset == "INI":
            drm_positions = [51, 66, 74, 75, 92, 95, 97, 118, 121, 122, 138, 140, 143, 145, 146, 147, 148, 151, 153, 155, 157, 163, 230, 232, 263]
            prot = "IN"
        elif dataset == "NNRTI":
            drm_positions = [90, 98, 100, 101, 103, 106, 108, 138, 179, 181, 188, 190, 221, 225, 227, 230, 234, 236, 238, 318, 348]
            prot = "RT"
        elif dataset == "NRTI":
            drm_positions = [41, 62, 65, 67, 68, 69, 70, 74, 75, 77, 115, 116, 151, 184, 210, 215, 219]
            prot = "RT"
        elif dataset == "PI":
            drm_positions = [10, 20, 24, 32, 33, 46, 47, 48, 50, 53, 54, 73, 74, 76, 82, 83, 84, 88, 89, 90]
            prot = "PR"

        # coverage_data_prot = coverage_data[coverage_data['Protein'] == prot]
        coverage_data_prot = coverage_data
        coverage_data_drm = coverage_data_prot[coverage_data_prot['Position'].isin(drm_positions)]
        #we calculate the % of covered positions
        sample_drm_coverage = coverage_data_drm[coverage_data_drm['Reads'] > read_cutoff].shape[0] / len(drm_positions) * 100

        handle.write(f"\n* {dataset} coverage\n")
        if sample_drm_coverage == 0:
            handle.write(f"\nNo DRM positions were covered for the {dataset} dataset.\n")
            continue

        robustness_step_plot(dataset, sample_drm_coverage, handle) #we plot the robustness step plot
    
    handle.write("\n\n\*Balanced accuracy takes into account the ability to predict both resistant and susceptible mutations. In contrast with accuracy, it is calculated as the mean of sensitivity and specificity.\n")





def robustness_step_plot(dataset: str , coverage_value: float, handle):
    '''Writes a step plot showing the balanced accuracy on the Y axis and the % of covered DRM positions on the X axis.
    We also plot with a red line the position of the analyzed sample coverage.
    INPUT:
    dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.
    coverage_value: coverage value of the analyzed sample as a %. i.e. 80% coverage = 80 .
    '''
    if coverage_value > 100 or coverage_value < 0:
        raise ValueError(f'Coverage value should be between 0 and 100. {coverage_value} provided.')
    
    if dataset == "INI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    elif dataset == "NNRTI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    elif dataset == "NRTI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    elif dataset == "PI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    else:
        raise ValueError(f'Unknown dataset: {dataset}')

    fig, axs = plt.subplots(1, 1, figsize=(10, 5))
    x = [0.5*m+1 for m in range(11)]
    y = step_data.drop('Drug', axis=1).groupby(["Miss_muts"]).mean()
    max_accuracy = y["Balanced_Accuracy"].max()
    coverage_step = 1-(int(coverage_value/10)+1)/10
    #we get the accuracy for the coverage step
    y_coverage = y.loc[coverage_step]["Balanced_Accuracy"]
    # print(f"Coverage step: {coverage_step}, y_coverage: {y_coverage}, max_accuracy: {max_accuracy}, coverage_value: {coverage_value}")
    handle.write(f"\nThe expected (balanced) accuracy for the {dataset} dataset drug resistance prediction is **{round(y_coverage, 2)}** at {round(coverage_value, 2)}% DRM positions coverage. The reported performance for a 100% coverage is {round(max_accuracy,2)}.\n\n")
    
    axs.axvline(x=((100-coverage_value)/10)*0.5 +1, color='red', linewidth = 2, label=f'Your sample ({round(coverage_value, 2)}% coverage)')#we plot a red line at the coverage value
    axs.step(x, y["Balanced_Accuracy"], where = 'post', alpha = 0.7, linewidth = 2, color='blue')
    axs.axhline(y=y_coverage, color='blue', linestyle='--', alpha = 0.35)
    axs.set_xticks([0.5*i+1 for i in range(11)])
    axs.set_xticklabels([100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0], fontsize=12)
    axs.set_ylim(0, 1)
    axs.set_xlabel('% of covered DRM positions', fontsize=15)
    axs.set_ylabel('Balanced Accuracy', fontsize=15)
    axs.set_title(f'Balanced accuracy for {dataset} class', fontsize=15)
    axs.legend(loc='upper right', fontsize=10)    #and we plot a small legend 

    path_to_plot = f'example_files/robustness_{dataset}.pdf'
    plt.savefig(path_to_plot)

    handle.write(r"\begin{center}")
    handle.write("\n")
    handle.write(r"\includegraphics[width=\textwidth]")
    handle.write("{%s}\n" % path_to_plot)
    handle.write("\end{center}\n")

    return f'example_files/coverage_robustness_{dataset}.pdf'

    
            
def write_report_md(mut_tsv, coverage_tsv, higher_cutoff: float = 0.15, lower_cutoff: float = 0.015, HIVDB:bool = True, LSR: bool = True, RF: bool = True):
    ''' Writes the report in markdown format.'''
    HIVDB_single_table = HIVDB_table(mut_tsv, lower_cutoff = lower_cutoff)

    md = open('example_files/report.md', 'w')
    md.write("# HIV-1 drug resistance report\n")
    md.write("## Drug resistance prediction\n")
    ensemble_predictions = ensemble_table(mut_tsv, higher_cutoff = higher_cutoff, HIVDB = HIVDB, LSR = LSR, RF = RF)
    write_ensemble_table(ensemble_predictions, cutoff = higher_cutoff, handle = md, HIVDB = HIVDB, LSR = LSR, RF = RF)
    write_coverage_disclaimer(coverage_tsv, handle = md)
    md.write("\\newpage\n")    
    md.write("\n## Single mutation annotation\n")
    md.write("### HIVDB single mutation relevant annotations\n\n")
    write_HIVDB_table(HIVDB_single_table, cutoff = lower_cutoff, handle = md, unknown = False, comments = True)
    md.write("\\newpage\n")
    md.write("\n### All single mutations\n")
    write_HIVDB_table(HIVDB_single_table, cutoff = lower_cutoff, handle = md, unknown = True, comments = False)





# ###Main function
# def main(input_mut: Annotated[str ,typer.Argument(help="Input string to process, should be a list of mutations separated by commas")],
#             HIVDB: Annotated[Optional[bool], typer.Option('--HIVDB', '-H', help='Use HIVDB method for prediction')] = None,
#             LSR: Annotated[Optional[bool], typer.Option('--LSR', '-L', help='Use Linear Regression method for prediction')] = None,
#             RF: Annotated[Optional[bool], typer.Option('--RF', '-R', help='Use Random Forest method for prediction')] = None):
#     # Example usage
#     mut_list = check_mut_input(input_mut)

#     annotations = HIVDB_singlemut_annot(mut_list, "NRTI") #we get the single mutation annotations
#     print("Single mutation annotations:")
#     print(annotations)

#     if HIVDB == None and LSR == None and RF == None:
#         ensemble_preds = ensemble_predictions(mut_list, "NRTI", HIVDB = True, LSR = True, RF = True)
#     else:
#         if HIVDB == None:
#             HIVDB = False
#         if LSR == None:
#             LSR = False
#         if RF == None:
#             RF = False

#         ensemble_preds = ensemble_predictions(mut_list, "NRTI", HIVDB = HIVDB, LSR = LSR, RF = RF)
#     print(ensemble_preds)


# input_mut = "K20R, V35I, T39A, M41L, S68G, K103N, I135T, M184V, T200K, Q207E, T215Y" 

# if __name__ == "__main__":
#     typer.run(main)

# table = pd.read_csv('example_files/mutation_freq.tsv', sep='\t')
# table = table[table['Freq'] > 0.05]
# print(extract_prot_mut_freq(table, 'RT'))
# input_list = [["K20R", "V35I", "T39A", "M41L", "S68G", "K103N", "I135T", "M184V", "T200K", "Q207E", "T215Y"], ["K20M", "V35I", "T39A", "M41L", "S68G", "K103N", "I135T", "M184V", "T200K", "Q207E", "T215Y"]]
# print(get_mutlist_comb(input_list))
# print(ensemble_table('example_files/mutation_freq.tsv', higher_cutoff = 0.05, HIVDB = True, LSR = True, RF = True))
# print(unify_mut_list(input_list))
write_report_md('example_files/mutation_freq.tsv', 'example_files/CAP257/week_54/alignments/coverage.tsv', higher_cutoff = 0.15, lower_cutoff = 0.015, HIVDB = True, LSR = True, RF = True)
# robustness_step_plot("NNRTI", 78)