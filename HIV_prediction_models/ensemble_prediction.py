###We will check the predictions from the different methods and use majority voting to decide the final prediction
import pandas as pd
from pathlib import Path


###We load the per drug test sets
test_drugs = ['RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
INIs=['RAL', 'EVG', 'DTG', 'BIC']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV']
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']

cutoff = {'FPV':3, 'ATV':3, 'IDV':3, 'LPV':9, 'NFV':3, 'SQV':3, 'TPV':2, 'DRV':10,
          '3TC':3, 'ABC':2, 'AZT':3, 'D4T':1.5, 'DDI':1.5, 'TDF':1.5, 
          'EFV':3, 'ETR':3, 'NVP':3, 'RPV':3,
          'RAL': 4, 'EVG': 4, 'DTG': 3.5, 'BIC':3.5}


##SMALL ENSEMBLE PREDICTIONS

ensemble_predictions = []
for drug in test_drugs:
    if drug in INIs:
        dataset = "INI"
    elif drug in PIs:
        dataset = "PI"
    elif drug in NNRTIs:
        dataset = "NNRTI"
    elif drug in NRTIs:
        dataset = "NRTI"

    test_set_data = pd.read_csv(f'../datasets/{dataset}/{dataset}_{drug}_5folds.tsv', sep='\t')

    for k in range(5):
        test_set = test_set_data[test_set_data["fold"]==k]
        print(drug, k, test_set_data.shape[0])

        hivdb_predictions = pd.read_csv(f'../method_predictions/HIVDB/HIVDB_{dataset}_predictions.tsv', sep='\t')
        lsr_comb_predictions = pd.read_csv(f'../method_predictions/linear_regression/{dataset}/{drug.replace('3TC','X3TC')}/LSR_I_{drug.replace('3TC','X3TC')}_fold_{k}_predictions.txt', sep=' ', index_col=0)
        rf_predictions = pd.read_csv(f'../method_predictions/random_forest/{dataset}/{drug}/random_forest_{drug}_fold_{k}_predictions.tsv', sep='\t')


        for seq in test_set["SeqID"].tolist():
            if pd.isnull(test_set[test_set["SeqID"]==seq]["CompMutList"].values[0]): ##WT predictions, we count as susceptible per default
                hivdb_pred = "Susceptible"
            else:
                hivdb_pred = hivdb_predictions[hivdb_predictions["seqID"]==seq][f"{drug.replace("3TC", "FTC")}"].values[0]
                hivdb_pred = hivdb_pred.replace('Potential Low-Level Resistance', 'Susceptible')
                hivdb_pred = hivdb_pred.replace('Low-Level Resistance', 'Susceptible')
                hivdb_pred = hivdb_pred.replace('Intermediate Resistance', 'Resistant')
                hivdb_pred = hivdb_pred.replace('High-Level Resistance', 'Resistant')

            lsr_comb_pred = float(lsr_comb_predictions[lsr_comb_predictions["SeqID"]==seq]["RF"].values[0])
            cutoff_value = cutoff[drug]
            if lsr_comb_pred > cutoff_value:
                lsr_comb_pred = "Resistant"
            else:
                lsr_comb_pred = "Susceptible"
            rf_pred = rf_predictions[rf_predictions["SeqID"]==seq]["Predictions"].values[0]

            most_common = max(set([hivdb_pred, lsr_comb_pred, rf_pred]), key = [hivdb_pred, lsr_comb_pred, rf_pred].count)


            real_pred = test_set[test_set["SeqID"]==seq][f"{drug}_res"].values[0]

            ensemble_predictions.append([drug, k, seq, hivdb_pred, lsr_comb_pred, rf_pred, most_common, real_pred])

ensemble_predictions_df = pd.DataFrame(ensemble_predictions, columns=["Drug", "Fold", "SeqID", "HIVDB", "LSR", "RandomForest", "HIVDB+LSR+RF", "Real"])
ensemble_predictions_df.to_csv("../method_predictions/small_ensemble_predictions.tsv", sep='\t', index=False)



###We do the same for the ensemble HIVDB+LSR_comb+LARS_comb+RF+CNN


ensemble_predictions = []
for drug in test_drugs:
    if drug in INIs:
        dataset = "INI"
    elif drug in PIs:
        dataset = "PI"
    elif drug in NNRTIs:
        dataset = "NNRTI"
    elif drug in NRTIs:
        dataset = "NRTI"

    test_set_data = pd.read_csv(f'../datasets/{dataset}/{dataset}_{drug}_5folds.tsv', sep='\t')

    for k in range(5):
        print(drug, k)
        test_set = test_set_data[test_set_data["fold"]==k]
        hivdb_predictions = pd.read_csv(f'../method_predictions/HIVDB/HIVDB_{dataset}_predictions.tsv', sep='\t')
        lsr_comb_predictions = pd.read_csv(f'../method_predictions/linear_regression/{dataset}/{drug.replace('3TC','X3TC')}/LSR_I_{drug.replace('3TC','X3TC')}_fold_{k}_predictions.txt', sep=' ', index_col=0)
        lars_comb_predictions = pd.read_csv(f'../method_predictions/linear_regression/{dataset}/{drug.replace('3TC','X3TC')}/LARS_I_{drug.replace('3TC','X3TC')}_fold_{k}_predictions.txt', sep=' ', index_col=0)
        rf_predictions = pd.read_csv(f'../method_predictions/random_forest/{dataset}/{drug}/random_forest_{drug}_fold_{k}_predictions.tsv', sep='\t')
        try:
            cnn_predictions = pd.read_csv(f'../method_predictions/cnn/cnn_{drug}_fold_{k}_predictions.tsv', sep='\t')
            cnn = True
        except:
            cnn = False

        for seq in test_set["SeqID"].tolist():
    
            if pd.isnull(test_set[test_set["SeqID"]==seq]["CompMutList"].values[0]): ##WT predictions, we count as susceptible per default
                hivdb_pred = "Susceptible"
            else:
                hivdb_pred = hivdb_predictions[hivdb_predictions["seqID"]==seq][f"{drug.replace("3TC", "FTC")}"].values[0]
                hivdb_pred = hivdb_pred.replace('Potential Low-Level Resistance', 'Susceptible')
                hivdb_pred = hivdb_pred.replace('Low-Level Resistance', 'Susceptible')
                hivdb_pred = hivdb_pred.replace('Intermediate Resistance', 'Resistant')
                hivdb_pred = hivdb_pred.replace('High-Level Resistance', 'Resistant')

            lsr_comb_pred = float(lsr_comb_predictions[lsr_comb_predictions["SeqID"]==seq]["RF"].values[0])
            cutoff_value = cutoff[drug]
            if lsr_comb_pred > cutoff_value:
                lsr_comb_pred = "Resistant"
            else:
                lsr_comb_pred = "Susceptible"
        
            lars_comb_pred = float(lars_comb_predictions[lars_comb_predictions["SeqID"]==seq]["lambda.min"].values[0])
            if lars_comb_pred > cutoff_value:
                lars_comb_pred = "Resistant"
            else:
                lars_comb_pred = "Susceptible"
            rf_pred = rf_predictions[rf_predictions["SeqID"]==seq]["Predictions"].values[0]
            try:
                cnn_pred = cnn_predictions[cnn_predictions["SeqID"]==seq]["CNN_Prediction"].values[0]
            except:
                cnn_pred = "NaN"
            #1 is Resistant and 0 is Susceptible
            if cnn_pred == 1:
                cnn_pred = "Resistant"
            elif cnn_pred == 0:
                cnn_pred = "Susceptible"


            #if CNN is not available we use the other 4 methods, but if we have the same counts for two labels, we do not use lars_comb prediction
            if cnn:
                most_common = max(set([hivdb_pred, lsr_comb_pred, lars_comb_pred, rf_pred, cnn_pred]), key = [hivdb_pred, lsr_comb_pred, lars_comb_pred, rf_pred, cnn_pred].count)
            else:
                if [hivdb_pred, lsr_comb_pred, lars_comb_pred, rf_pred].count("Resistant") == [hivdb_pred, lsr_comb_pred, lars_comb_pred, rf_pred].count("Susceptible"):
                    most_common = max(set([hivdb_pred, lsr_comb_pred, rf_pred]), key = [hivdb_pred, lsr_comb_pred, rf_pred].count)
                else:
                    most_common = max(set([hivdb_pred, lsr_comb_pred, lars_comb_pred, rf_pred]), key = [hivdb_pred, lsr_comb_pred, lars_comb_pred, rf_pred].count)

            real_pred = test_set[test_set["SeqID"]==seq][f"{drug}_res"].values[0]

            ensemble_predictions.append([drug, k, seq, hivdb_pred, lsr_comb_pred, lars_comb_pred, rf_pred, cnn_pred, most_common, real_pred])

ensemble_predictions_df = pd.DataFrame(ensemble_predictions, columns=["Drug", "Fold", "SeqID", "HIVDB", "LSR", "LARS", "RandomForest", "CNN","HIVDB+LSR+LARS+RF+CNN", "Real"])
ensemble_predictions_df.to_csv("../method_predictions/large_ensemble_predictions.tsv", sep='\t', index=False)
