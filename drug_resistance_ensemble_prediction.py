import pandas as pd
import re
from itertools import product

##Drugs per drug class used in this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']

def get_mutlist_comb_dict(mut_dict: dict):
    '''Outputs a list of dicts with all possible dictionaries having only one mutation per position.
    INPUT:
    mut_dict: dictionary with the mutations as keys and the frequencies as values.
    
    OUTPUT:
    comb_list: list of dictionaries with all possible combinations of mutations. List of dicts.
    '''
    comb_list = []
    mutation_list = list(mut_dict.keys())
    pos_list = [re.sub(r'\D+', '', mut) for mut in mutation_list] #we get the positions of the mutations

    pos_dict = {}    #we make a dictionary with the positions as keys and a list of their index as values
    for i, pos in enumerate(pos_list):
        if pos not in pos_dict:
            pos_dict[pos] = [i]
        else:
            pos_dict[pos].append(i)
    #we create a list of all possible combinations of mutations
    for i in range(len(pos_dict)):
        pos = list(pos_dict.keys())[i]
        if len(pos_dict[pos]) > 1:
            comb_list.append([mutation_list[j] for j in pos_dict[pos]])
        else:
            comb_list.append([mutation_list[pos_dict[pos][0]]])
    
    comb_list = list(product(*comb_list))#we create all the combinations of the mutations
    comb_list = [list(comb) for comb in comb_list]

    comb_list_dict = []
    for comb in comb_list:
        comb_dict = {}
        for mut in comb:
            comb_dict[mut] = mut_dict[mut]
        comb_list_dict.append(comb_dict)

    return comb_list_dict


###HIVDB offline implementation for the prediction of the resistance
def HIVDB_pred(mut_dic: dict, hivdb_data_dic: dict, drug: str):
    '''Calculates the HIVDB score given a dictionary of mutations using the HIVDB scores for mutations and combinations extracted from StanfordHIVDB HIVDB program.
        Scores downloaded on 10th April 2025.
        
        INPUT:
        mut_dic: dictionary with the positions as keys and 'Mut' and 'Freq' values.
        drug: drug to be used. 'RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV'.

        OUTPUT:
        output_1: dictionary with the predicted resistance labels for each drug in the drug class. SR gives Susceptible/Resistant labels, HIVDB_five_labels gives the HIVDB labeling system (5 labels).
        output_2: dictionary with the frequencies, annotations and comments from HIVDB.
    '''
    output_1 = dict()
    
    if drug in INIs:
        dataset = "INI"
        protein = 'IN'
    elif drug in NNRTIs:
        dataset = "NNRTI"
        protein = 'RT'
    elif drug in NRTIs:
        dataset = "NRTI"
        protein = 'RT'
    elif drug in PIs:
        dataset = "PI"
        protein = 'PR'

    drug_score = 0

    for mutation in mut_dic.keys():
        reference = mut_dic[mutation]["Ref"]
        #we add the single mutation scores
        if reference + mutation in hivdb_data_dic[dataset]['single_mut'].keys():
            drug_score += float(hivdb_data_dic[dataset]['single_mut'][reference+mutation][drug])
                
    #we add the combination scores
    for combination_rule in hivdb_data_dic[dataset]['combination'].keys():
        combination_rule_list = combination_rule.split('+')
        all_in = True
        for mut_c_rule in combination_rule_list:
            #we separate the mutation from the position
            position = re.sub(r'\D+', '', mut_c_rule)
            mut_aa = re.sub(r'\d+', '', mut_c_rule)
            # mut_aa_1 = mut_aa[0]
            mut_aa_2 = mut_aa[1:]
            if len(mut_aa_2) == 1:
                if mut_c_rule[1:] not in mut_dic:
                    all_in = False
                    break
            else:
                posib_mut = False
                for l in mut_aa_2:
                    if position+l in mut_dic:
                        posib_mut = True
                        break
                    else:
                        posib_mut = False

                if not posib_mut:
                    all_in = False
                    break

        if all_in:
            drug_score += float(hivdb_data_dic[dataset]['combination'][combination_rule][drug])
                    
        output_1['HIVDB_score'] = drug_score
        output_1["SR"] = "Susceptible" if drug_score < 30 else "Resistant"
        output_1['HIVDB_five_labels'] = "Susceptible" if drug_score < 10 else "Potential Low-Level Resistance" if drug_score < 15 else "Low-Level Resistance" if drug_score < 30 else "Intermediate Resistance" if drug_score < 60 else "High-Level Resistance"
        
    return output_1



###Linear regression implementation
def LSR_res_pred(mut_dic: dict, lsr_coefficients: dict, drug:str):
    '''Calculates the resistance factor given a set of mutations using the coefficients obtained from LSR with TSM features.
        INPUT:
        mut_dic: dictionary with the mutations as keys and the frequencies as values.
        lsr_coefficients: dictionary with the coefficients for each drug.
        drug: drug to be used. 'RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV'.

        OUTPUT:
        LSR_predictions: dictionary with the predicted resistance labels for each drug in the drug class. Labels are Susceptible/Resistant, RF is the resistance factor.      
    '''
    cutoff = {'FPV':3, 'ATV':3, 'IDV':3, 'LPV':9, 'NFV':3, 'SQV':3, 'TPV':2, 'DRV':10, #resistance fold cutoff defined to decide drug resistance
              '3TC':3, 'ABC':2, 'AZT':3, 'D4T':1.5, 'DDI':1.5, 'TDF':1.5, 
              'EFV':3, 'ETR':3, 'NVP':3, 'RPV':3,
              'RAL': 4, 'EVG': 4, 'DTG': 3.5, 'BIC':3.5}

    LSR_predictions = dict()
    if drug in INIs:
        drug_class = INIs
        protein = 'IN'
    elif drug in NNRTIs:
        drug_class = NNRTIs
        protein = 'RT'
    elif drug in NRTIs:
        drug_class = NRTIs
        protein = 'RT'
    elif drug in PIs:
        drug_class = PIs
        protein = 'PR'

    coefficients = lsr_coefficients[drug]
    RF = 0 + float(coefficients['intercept']) #resistance factor

    ##we multiply by 1 the coefficients of the mutations that are present in the mut_list
    for mut in mut_dic.keys():
        if mut in coefficients:
            RF += float(coefficients[mut])  
            
    ##we parse the mutation combinations
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
                    if m not in mut_dic:
                        all_in = False
                        break
                else:
                    posib_mut = False
                    for l in letters:
                        if position+l in mut_dic:
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
    LSR_predictions['RF'] = RF
    ##we apply the cutoff and predict the resistance
    if RF > float(cutoff[drug]):
        label = "Resistant"
    else:
        label = "Susceptible"
    LSR_predictions['SR'] = label

    return LSR_predictions


###Random Forest implementation
def RandomForest_pred(mut_dic: dict, rf_models: dict, drug: str):
    '''Predicts Resistance/Susceptible labels with Random Forest models.
        INPUT:
        mut_dic: dictionary with the mutations as keys and the frequencies as values.
        rf_models: dictionary with the models for each drug.
        drug: drug to be used. 'RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV'.

        OUTPUT:
        RF_predictions: dictionary with the predicted resistance labels for each drug in the drug class. Labels are Susceptible/Resistant.    
    '''
    RF_predictions = dict()
    if drug in INIs:
        drug_class = INIs
        protein = 'IN'
    elif drug in NNRTIs:
        drug_class = NNRTIs
        protein = 'RT'
    elif drug in NRTIs:
        drug_class = NRTIs
        protein = 'RT'
    elif drug in PIs:
        drug_class = PIs
        protein = 'PR'
        
            
    if protein == "PR":
        ref_seq = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    elif protein == "RT":
        ref_seq = "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWT"
    elif protein == "IN":
        ref_seq = "FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"

    X = pd.DataFrame([list(ref_seq)], columns=[f'V{i+1}' for i in range(len(ref_seq))]) #we create a dataframe with a column for each position in the sequence as input
            
    for mut in mut_dic.keys():
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
        
    rf_model = rf_models[drug]

    #we only keep the columns that are in the model
    X_drug = X[rf_model.feature_names_in_] #we only keep the columns that are in the model

    RF_label = rf_model.predict(X_drug) #we predict the label
    RF_predictions = RF_label[0] #we save the label

    return RF_predictions


###Ensemble prediction
def ensemble_predictions(mut_dict:dict, hivdb_data, lsr_coefficients, rf_models, cutoff, HIVDB:bool = True, LSR: bool = True, RF: bool = True):
    '''Generates an ensemble prediction by majority voting of the HIVDB, Linear Regression and Random Forest predictions.

    INPUT:
    mut_dict: dictionary with the mutations as keys and the frequencies as values.
    hivdb_data: dictionary with the HIVDB data for each drug class.
    lsr_coefficients: dictionary with the LSR coefficients for each drug.
    rf_models: dictionary with the models for each drug.
    HIVDB: if True, the HIVDB prediction will be used.
    LSR: if True, the Linear Regression prediction will be used.
    RF: if True, the Random Forest prediction will be used.

    OUTPUT:
    ensemble_predictions: dictionary with the predicted resistance labels for each drug in the drug class for all listed methods.
    single_mut_annotations: dictionary with the annotations for each mutation in the input list.
    '''

    if not HIVDB and not LSR and not RF:
        print("At least one method should be selected")
        return
    
    ensemble_predictions = dict()
    ensemble_predictions["cutoff"] = cutoff
    
    for dataset in ['INI', 'NNRTI', 'NRTI', 'PI']:
        if dataset == "INI":
            drug_class = INIs
            protein = 'IN'
        elif dataset == "NNRTI":
            drug_class = NNRTIs
            protein = 'RT'
        elif dataset == "NRTI":
            drug_class = NRTIs
            protein = 'RT'
        elif dataset == "PI":
            drug_class = PIs
            protein = 'PR'

        ensemble_predictions[dataset] = dict()
        mut_dict_protein_unfiltered = mut_dict[protein]
        mut_dict_protein = freq_filter(mut_dict_protein_unfiltered, cutoff) #we filter the mutations by frequency
        mut_dict_combinations = get_mutlist_comb_dict(mut_dict_protein) #we get the combinations of mutations

        for mut_dict_protein in mut_dict_combinations:
            mut_combination_string = get_mutation_string(mut_dict_protein)
            # print(f"Mutations for {dataset}: {mut_combination_string}")
            ensemble_predictions[dataset][mut_combination_string] = dict()

            for drug in drug_class:
                ensemble_predictions[dataset][mut_combination_string][drug] = dict()

                ##HIVDB
                if HIVDB:
                    HIVDB_results = HIVDB_pred(mut_dict_protein, hivdb_data, drug)
                    ensemble_predictions[dataset][mut_combination_string][drug]['HIVDB'] = HIVDB_results['SR']
                    ensemble_predictions[dataset][mut_combination_string][drug]['HIVDB_five_labels'] = HIVDB_results['HIVDB_five_labels']
                
                ##Linear Regression
                if LSR:
                    LSR_results = LSR_res_pred(mut_dict_protein, lsr_coefficients, drug)
                    ensemble_predictions[dataset][mut_combination_string][drug]['LSR'] = LSR_results['SR']
                    ensemble_predictions[dataset][mut_combination_string][drug]['LSR_RF'] = LSR_results['RF']
                
                ##Random Forest
                if RF:
                    RF_results = RandomForest_pred(mut_dict_protein, rf_models, drug)
                    ensemble_predictions[dataset][mut_combination_string][drug]['RandomForest'] = RF_results

                ##Ensemble prediction
                if HIVDB and LSR and RF:
                    ensemble_predictions[dataset][mut_combination_string][drug]['Ensemble'] = pd.Series([ensemble_predictions[dataset][mut_combination_string][drug]['HIVDB'], ensemble_predictions[dataset][mut_combination_string][drug]['LSR'], ensemble_predictions[dataset][mut_combination_string][drug]['RandomForest']]).mode()[0]

    return ensemble_predictions

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

def get_mutation_string(mut_dict: dict):
    '''Generates a string with the mutations in the input dictionary.
    
    INPUT:
    mut_dict: dictionary with the mutations as keys and the frequencies as values.
    
    OUTPUT:
    mutation_string: string with the mutations in the input dictionary.
    '''
    mutation_string = ''
    mutation_list = list(mut_dict.keys())
    mutation_list.sort(key=lambda x: int(re.sub(r'\D+', '', x)))
    for mut in mutation_list:
        mutation_string += mut + '+'
    
    return mutation_string[:-1] #we remove the last comma and space

def get_ensemble_df(ensemble_dict_list: list, sample_ids: list=None):
    '''Creates pandas DataFrame with the ensemble predictions for each sample.
    INPUT:
    ensemble_dict_list: list of dictionaries with the ensemble predictions for each sample.
    sample_ids: list with the sample ids.
    OUTPUT:
    ensemble_df: pandas DataFrame with the ensemble predictions for each sample.
    '''
    ensemble_df = []
    if sample_ids is None:
        cutoff = ensemble_dict_list['cutoff']
        for dataset in ensemble_dict_list.keys():
            if dataset == 'cutoff':
                continue
            for mut_combination in ensemble_dict_list[dataset].keys():
                for drug in ensemble_dict_list[dataset][mut_combination].keys():
                    HIVDB, HIVDB_5_labels = ensemble_dict_list[dataset][mut_combination][drug]['HIVDB'], ensemble_dict_list[dataset][mut_combination][drug]['HIVDB_five_labels']
                    LSR, LSR_RF = ensemble_dict_list[dataset][mut_combination][drug]['LSR'], ensemble_dict_list[dataset][mut_combination][drug]['LSR_RF']
                    RF = ensemble_dict_list[dataset][mut_combination][drug]['RandomForest']
                    ensemble = ensemble_dict_list[dataset][mut_combination][drug]['Ensemble']
                    ensemble_df.append([dataset, mut_combination, drug, HIVDB, HIVDB_5_labels, LSR, LSR_RF, RF, ensemble, cutoff])
        ensemble_df = pd.DataFrame(ensemble_df, columns=['Drug_Class', 'Mutations', 'Drug', 'HIVDB', 'HIVDB_5_labels', 'LSR', 'LSR_RF', 'RF', 'Ensemble', 'Cutoff'])
    else:
        for sample_ID, ensemble_pred in zip(sample_ids, ensemble_dict_list):
            cutoff = ensemble_pred['cutoff']
            for dataset in ensemble_pred.keys():
                if dataset == 'cutoff':
                    continue
                for mut_combination in ensemble_pred[dataset].keys():
                    for drug in ensemble_pred[dataset][mut_combination].keys():
                        HIVDB, HIVDB_5_labels = ensemble_pred[dataset][mut_combination][drug]['HIVDB'], ensemble_pred[dataset][mut_combination][drug]['HIVDB_five_labels']
                        LSR, LSR_RF = ensemble_pred[dataset][mut_combination][drug]['LSR'], ensemble_pred[dataset][mut_combination][drug]['LSR_RF']
                        RF = ensemble_pred[dataset][mut_combination][drug]['RandomForest']
                        ensemble = ensemble_pred[dataset][mut_combination][drug]['Ensemble']
                        ensemble_df.append([sample_ID, dataset, mut_combination, drug, HIVDB, HIVDB_5_labels, LSR, LSR_RF, RF, ensemble, cutoff])

        ensemble_df = pd.DataFrame(ensemble_df, columns=['Sample_ID', 'Drug_Class', 'Mutations', 'Drug', 'HIVDB', 'HIVDB_5_labels', 'LSR', 'LSR_RF', 'RF', 'Ensemble', 'Cutoff'])
    return ensemble_df
        