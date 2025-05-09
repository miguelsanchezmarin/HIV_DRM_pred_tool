import pandas as pd
import re
import pickle
from itertools import product

##Drugs per drug class used in this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']


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

def df_in_list(df, list_of_dfs):
    '''Checks if a dataframe is in a list of dataframes.
    Outputs the index of the dataframe in the list if it is found.'''
    df = df.sort_index(axis=1).reset_index(drop=True)
    for i, other_df in enumerate(list_of_dfs):
        other_df = other_df.sort_index(axis=1).reset_index(drop=True)
        if df.equals(other_df):
            return True, i
    return False, None

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