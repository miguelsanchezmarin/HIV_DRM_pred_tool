###We implement the trained random forest models in Python
import pandas as pd
from pathlib import Path
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import confusion_matrix



test_drugs = ['RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
INIs=['RAL', 'EVG', 'DTG', 'BIC']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV']
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']

for drug in test_drugs:
    if drug in PIs:
        dataset = "PI"
    elif drug in NNRTIs:
        dataset = "NNRTI"
    elif drug in NRTIs:
        dataset = "NRTI"
    elif drug in INIs:
        dataset = "INI"

    # we load the data
    data_set = pd.read_csv(f'../datasets/{dataset}/{dataset}_{drug}_5folds.tsv', sep='\t')
    #we only keep SeqID, drug_res and CompMutList columns
    data_set = data_set[['SeqID', f'{drug}_res', 'CompMutList', 'fold']]

    #we get the res values as the Y
    Y = data_set[f'{drug}_res']

    if dataset == "PI":
        ref_seq = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    elif dataset == "NNRTI" or dataset == "NRTI":
        ref_seq = "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWT"
    else:
        ref_seq = "FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"

    #we build the X, a datafram with a column per position (P1, P2, P3, ...) and a row per sequence
    #we first set the reference sequence in the data frame
    X = pd.DataFrame([list(ref_seq)], columns=[f'V{i+1}' for i in range(len(ref_seq))])
    #we multiply it by the number of sequences
    X = pd.concat([X]*len(data_set), ignore_index=True)

    for seq in range(len(data_set)):
        mut_list = data_set.loc[seq, 'CompMutList']
        try:
            mut_list = mut_list.split(', ')
            mut_list = [m[1:] for m in mut_list]
        except:
            mut_list = []
        
        for mut in mut_list:
            pos = int(mut[:-1])
            aa = mut[-1]
            X.loc[seq, f'V{pos}'] = aa
    

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

    #we bind the Y to the X dataframe and the SeqIDs
    XY = pd.concat([data_set[['SeqID', f'{drug}_res']], X], axis=1)

    for fold in range(5):

        #we split in train and test sets
        test_set = data_set[data_set['fold'] == fold]
        train_set = data_set[data_set['fold'] != fold]      

        trainXY = XY[XY['SeqID'].isin(train_set['SeqID'])]
        testXY = XY[XY['SeqID'].isin(test_set['SeqID'])]
        trainXY_filtered = trainXY.drop(['SeqID', f'{drug}_res'], axis=1)
        testXY_filtered = testXY.drop(['SeqID', f'{drug}_res'], axis=1)

        #now we do the VarianceThreshold
        selector = VarianceThreshold(threshold=0.01)
        train_selector = selector.fit(trainXY_filtered)
        trainXY_filtered = trainXY_filtered[trainXY_filtered.columns[selector.get_support(indices=True)]]
        testXY_filtered = testXY_filtered[testXY_filtered.columns[selector.get_support(indices=True)]]
        #we add back the SeqID and the drug_res columns
        trainXY_filtered = pd.concat([trainXY[['SeqID', f'{drug}_res']], trainXY_filtered], axis=1)
        testXY_filtered = pd.concat([testXY[['SeqID', f'{drug}_res']], testXY_filtered], axis=1)

        X_train = trainXY_filtered.drop(['SeqID', f'{drug}_res'], axis=1)
        Y_train = trainXY_filtered[f'{drug}_res']
        X_test = testXY_filtered.drop(['SeqID', f'{drug}_res'], axis=1)
        Y_test = testXY_filtered[f'{drug}_res']

        #we train the model
        n_trees = 500 #extracted from Raposo et al.
        mtry = 4
        rf = RandomForestClassifier(n_estimators=n_trees, max_features=mtry)
        rf.fit(X_train, Y_train)

        #we predict the test set
        predictions = rf.predict(X_test)

        #we save a .tsv with the predictions and their SeqIDs
        predictions_df = pd.DataFrame(testXY['SeqID'])
        predictions_df['Predictions'] = predictions
        predictions_df.to_csv(f'../method_predictions/random_forest/random_forest_{drug}_fold_{fold}_predictions.tsv', sep='\t', index=False)
        print(f"Fold {fold} for {drug} done")

        #we calculate the confusion matrix
        print(f"Confusion matrix for {drug}:")
        print(confusion_matrix(Y_test, predictions))

        # #we save the probabilities
        # probabilities = rf.predict_proba(X_test)
        # probabilities_df = pd.DataFrame(probabilities, columns=['Prob_Susceptible', 'Prob_Resistant'])
        # probabilities_df['SeqID'] = testXY['SeqID'].tolist()
        # probabilities_df['Predictions'] = predictions
        # probabilities_df.to_csv(Path(DATADIR, f'random_forest/RF_{drug}_fold_{fold}_probabilities.tsv'), sep='\t', index=False)

        # #we import the model as a pkl file
        # with open(Path(DATADIR,f"random_forest/models_RF/random_forest_python_{drug}_RF_model_fold_{fold}.pkl"), "wb") as file:
        #     pickle.dump(rf, file)
