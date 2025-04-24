###Labels a list of mutations as Resistant or Susceptible for different anti-HIV drugs
import pandas as pd
import re
import pickle
import typer
from typing_extensions import Annotated, Optional

##Drugs per drug class used in this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']

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
            if re.match(r'^\d+[A-Z]+$', comm) and comm.startswith(position):
                comm_ambig = re.sub(r'\d+', '', comm)
                if mut_aa in comm_ambig:
                    mut_annot.append([comments_file['Condition'].values[i], 
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
    elif score == "HIVDB": #we label them follwoing the original HIVDB labeling
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

###Main function
def main(input_mut: Annotated[str ,typer.Argument(help="Input string to process, should be a list of mutations separated by commas")],
            HIVDB: Annotated[Optional[bool], typer.Option('--HIVDB', '-H', help='Use HIVDB method for prediction')] = None,
            LSR: Annotated[Optional[bool], typer.Option('--LSR', '-L', help='Use Linear Regression method for prediction')] = None,
            RF: Annotated[Optional[bool], typer.Option('--RF', '-R', help='Use Random Forest method for prediction')] = None):
    # Example usage
    mut_list = check_mut_input(input_mut)

    annotations = HIVDB_singlemut_annot(mut_list, "NRTI") #we get the single mutation annotations
    print("Single mutation annotations:")
    print(annotations)

    if HIVDB == None and LSR == None and RF == None:
        ensemble_preds = ensemble_predictions(mut_list, "NRTI", HIVDB = True, LSR = True, RF = True)
    else:
        if HIVDB == None:
            HIVDB = False
        if LSR == None:
            LSR = False
        if RF == None:
            RF = False

        ensemble_preds = ensemble_predictions(mut_list, "NRTI", HIVDB = HIVDB, LSR = LSR, RF = RF)
    print(ensemble_preds)


input_mut = "K20R, V35I, T39A, M41L, S68G, K103N, I135T, M184V, T200K, Q207E, T215Y" 

if __name__ == "__main__":
    typer.run(main)