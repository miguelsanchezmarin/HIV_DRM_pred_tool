###We will plot the ROC and the precision-recall curve for each model
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc

DATADIR = '/cluster/home/msanche/MiguelSM_2025/HIV_DRM_Prediction'

test_drugs = ['RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
test_drugs_no_INIs = ['EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
INIs=['RAL', 'EVG', 'DTG', 'BIC']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] #RPV test set is empty
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']

# ensemble_predictions = pd.read_csv(Path(DATADIR, f'ensemble_learning/ensemble_predictions_5folds.tsv'), sep='\t')
# ensemble_5_predictions = pd.read_csv(Path(DATADIR, f'ensemble_learning/ensemble_predictions_hivdb_lsrcomb_larscomb_rf_cnn_5folds.tsv'), sep='\t')

# ###We will create a dataframe with all the predictions and real value for the different drugs
# predictions_list = []
# for drug in test_drugs: #we loop over the drugs
#     #we define the dataset
#     if drug in INIs:
#         dataset = 'INI'
#     elif drug in PIs:
#         dataset = 'PI'
#     elif drug in NNRTIs:
#         dataset = 'NNRTI'
#     elif drug in NRTIs:
#         dataset = 'NRTI'

#     drug_test_set = pd.read_csv(Path(DATADIR, f'{dataset}_{drug}_5folds.tsv'), sep='\t')
#     hivdb_pred = pd.read_csv(Path(DATADIR, f'HIVDB/{dataset.lower()}_HIVDB_fullpreds.tsv'), sep='\t')


#     for fold in range(5): #we loop over the folds
#         ##first we read the test set
#         drug_test_set_fold = drug_test_set[drug_test_set['fold']==fold]

#         #and the predictions
#         lsr_preds_fold = pd.read_csv(Path(DATADIR, f'linear_regression_HIVDB/SeqID_preds_{drug.replace('3TC','X3TC')}_tsm_fold_{fold}.txt'), sep=' ', index_col=0)
#         lsr_comb_preds_fold = pd.read_csv(Path(DATADIR, f'linear_regression_HIVDB/SeqID_preds_{drug.replace('3TC','X3TC')}_combinations_tsm_fold_{fold}.txt'), sep=' ', index_col=0)
#         lars_preds_fold = pd.read_csv(Path(DATADIR, f'linear_regression_HIVDB/SeqID_preds_LARS_{drug.replace('3TC','X3TC')}_tsm_fold_{fold}.txt'), sep=' ', index_col=0)
#         lars_comb_preds_fold = pd.read_csv(Path(DATADIR, f'linear_regression_HIVDB/SeqID_preds_LARS_{drug.replace('3TC','X3TC')}_combinations_tsm_fold_{fold}.txt'), sep=' ', index_col=0)

#         for n in range(drug_test_set_fold.shape[0]):
#             seq_id, true_label = drug_test_set_fold.iloc[n]['SeqID'], drug_test_set_fold.iloc[n][f'{drug}_res'] #we read SeqID and true label

#             #first we read the HIVDB predictions   
#             try:         
#                 hivdb_pred_label = hivdb_pred[hivdb_pred['seqID']==seq_id][f'{drug.replace('3TC', 'FTC')}'].values[0]
#             except:
#                 hivdb_pred_label = 'NA'
#             hivdb_pred_label = hivdb_pred_label.replace('Potential Low-Level Resistance', 'Susceptible')
#             hivdb_pred_label = hivdb_pred_label.replace('Low-Level Resistance', 'Susceptible')
#             hivdb_pred_label = hivdb_pred_label.replace('Intermediate Resistance', 'Resistant')
#             hivdb_pred_label = hivdb_pred_label.replace('High-Level Resistance', 'Resistant')
#             print("HIVDB:", hivdb_pred_label)

#             #now LSR, LSR_comb, LARS, LARS_comb predictions
#             lsr_pred_label = lsr_preds_fold[lsr_preds_fold['SeqID']==seq_id]['RF'].values[0]
#             lsr_comb_pred_label = lsr_comb_preds_fold[lsr_comb_preds_fold['SeqID']==seq_id]['RF'].values[0]
#             lars_pred_label = lars_preds_fold[lars_preds_fold['SeqID']==seq_id]['lambda.min'].values[0]
#             lars_comb_pred_label = lars_comb_preds_fold[lars_comb_preds_fold['SeqID']==seq_id]['lambda.min'].values[0]
#             print("LSR:", lsr_pred_label, "LSR_Comb:", lsr_comb_pred_label, "LARS:", lars_pred_label, "LARS_Comb:", lars_comb_pred_label)

#             #RandomForest
#             rf_predictions = pd.read_csv(Path(DATADIR, f"random_forest/RF_{drug}_fold_{fold}_predictions.tsv"), sep='\t')
#             rf_pred_label = rf_predictions[rf_predictions['SeqID']==seq_id]['Predictions'].values[0]
#             print("RF:", rf_pred_label)

#             ##CNN
#             if drug in test_drugs_no_INIs:
#                 cnn_predictions = pd.read_csv(Path(DATADIR, f"deep_learning_HIV/{drug}_preds_test_fold_{fold}.tsv"), sep='\t')
#                 cnn_pred_label = cnn_predictions[cnn_predictions['SeqID']==seq_id]['CNN_Prediction'].values[0]
#             else:
#                 cnn_pred_label = 'NA'
#             print("CNN:", cnn_pred_label)

#             ###Ensemble predictions
#             ensemble_pred_label = ensemble_predictions[(ensemble_predictions['SeqID']==seq_id)&(ensemble_predictions['Drug']==drug)][f'HIVDB+LSR+RF'].values[0]
#             ensemble_5_pred_label = ensemble_5_predictions[(ensemble_5_predictions['SeqID']==seq_id)&(ensemble_5_predictions['Drug']==drug)][f'HIVDB+LSR+LARS+RF+CNN'].values[0]
#             print("Ensemble:", ensemble_pred_label, "Ensemble_5:", ensemble_5_pred_label)

#             predictions_list.append([drug, fold, seq_id, true_label, hivdb_pred_label, lsr_pred_label, lsr_comb_pred_label, lars_pred_label, lars_comb_pred_label, rf_pred_label, cnn_pred_label, ensemble_pred_label, ensemble_5_pred_label])

# predictions_df = pd.DataFrame(predictions_list, columns=['Drug', 'Fold', 'SeqID', 'True_Label', 'HIVDB', 'LSR', 'LSR_Comb', 'LARS', 'LARS_Comb', 'RF', 'CNN', 'HIVDB+LSR_C+RF', 'HIVDB+LSR_C+LARS_C+RF+CNN'])
# predictions_df.to_csv(Path(DATADIR, 'all_predictions_df.tsv'), sep='\t', index=False)

# #we print how many NAs in HIVDB
# print("HIVDB NAs:", predictions_df[predictions_df['HIVDB']=='NA'].shape[0])

# predictions_df = pd.read_csv(Path(DATADIR, 'all_predictions_df.tsv'), sep='\t')
# # #in the CNN column, we change 0 by 'Susceptible' and 1 by 'Resistant'
# # predictions_df['CNN'] = predictions_df['CNN'].replace(0, 'Susceptible')
# # predictions_df['CNN'] = predictions_df['CNN'].replace(1, 'Resistant')

# cutoff_dic = {'FPV':3, 'ATV':3, 'IDV':3, 'LPV':9, 'NFV':3, 'SQV':3, 'TPV':2, 'DRV':10,
#               '3TC':3, 'ABC':2, 'AZT':3, 'D4T':1.5, 'DDI':1.5, 'TDF':1.5, 
#               'EFV':3, 'ETR':3, 'NVP':3, 'RPV':3,
#               'RAL': 4, 'EVG': 4, 'DTG': 3.5, 'BIC':3.5}

# for n in range(predictions_df.shape[0]):
#     drug = predictions_df.iloc[n]['Drug']
#     seqid = predictions_df.iloc[n]['SeqID']
#     fold = predictions_df.iloc[n]['Fold']
#     cutoff = cutoff_dic[drug]

#     rf_probabilities_df = pd.read_csv(Path(DATADIR, f"random_forest/RF_{drug}_fold_{fold}_probabilities.tsv"), sep='\t')
#     res_prob = rf_probabilities_df[rf_probabilities_df['SeqID']==seqid]['Prob_Susceptible'].values[0] #it is really the probability of being resistant but the labels are wrong
#     predictions_df.at[n, 'RF'] = res_prob
#     if drug in test_drugs_no_INIs:
#         cnn_probabilities_df = pd.read_csv(Path(DATADIR, f"deep_learning_HIV/{drug}_preds_test_fold_{fold}.tsv"), sep='\t') 
#         cnn_prob = cnn_probabilities_df[cnn_probabilities_df['SeqID']==seqid]['CNN_Prob'].values[0]
#         predictions_df.at[n, 'CNN'] = cnn_prob
    
#     # if float(predictions_df.iloc[n]['LSR'])<cutoff:
#     #     # predictions_df.at[n, 'LSR'] = 'Susceptible'
#     #     predictions_df.at[n, 'LSR'] = 0
#     # else:
#     #     # predictions_df.at[n, 'LSR'] = 'Resistant'
#     #     predictions_df.at[n, 'LSR'] = 1
#     # if float(predictions_df.iloc[n]['LSR_Comb'])<cutoff:
#     #     # predictions_df.at[n, 'LSR_Comb'] = 'Susceptible'
#     #     predictions_df.at[n, 'LSR_Comb'] = 0
#     # else:
#     #     # predictions_df.at[n, 'LSR_Comb'] = 'Resistant'
#     #     predictions_df.at[n, 'LSR_Comb'] = 1
#     # if float(predictions_df.iloc[n]['LARS'])<cutoff:
#     #     # predictions_df.at[n, 'LARS'] = 'Susceptible'
#     #     predictions_df.at[n, 'LARS'] = 0
#     # else:
#     #     # predictions_df.at[n, 'LARS'] = 'Resistant'
#     #     predictions_df.at[n, 'LARS'] = 1
#     # if float(predictions_df.iloc[n]['LARS_Comb'])<cutoff:
#     #     # predictions_df.at[n, 'LARS_Comb'] = 'Susceptible'
#     #     predictions_df.at[n, 'LARS_Comb'] = 0
#     # else:
#     #     # predictions_df.at[n, 'LARS_Comb'] = 'Resistant'
#     #     predictions_df.at[n, 'LARS_Comb'] = 1

#     predictions_df.at[n, 'LSR'] = float(predictions_df.iloc[n]['LSR'])
#     predictions_df.at[n, 'LSR_Comb'] = float(predictions_df.iloc[n]['LSR_Comb'])
#     predictions_df.at[n, 'LARS'] = float(predictions_df.iloc[n]['LARS'])
#     predictions_df.at[n, 'LARS_Comb'] = float(predictions_df.iloc[n]['LARS_Comb'])



#     #we change for the rest of the models Resistant by 1 and Susceptible by 0
#     for model in ['True_Label', 'HIVDB', 'HIVDB+LSR_C+RF', 'HIVDB+LSR_C+LARS_C+RF+CNN']:
#         if predictions_df.iloc[n][model]=='Resistant':
#             predictions_df.at[n, model] = 1
#         elif predictions_df.iloc[n][model]=='Susceptible':
#             predictions_df.at[n, model] = 0
    
# #we fill the missing values as NA
# predictions_df = predictions_df.fillna('NA')

# predictions_df.to_csv(Path(DATADIR, 'all_predictions_df_cutoff.tsv'), sep='\t', index=False)

predictions_df = pd.read_csv(Path(DATADIR, 'all_predictions_df_cutoff.tsv'), sep='\t')


#we convert to integers columns LSR, LSR_Comb, LARS and LARS_Comb and CNN when not NA
predictions_df['LSR'] = predictions_df['LSR'].astype(float)
predictions_df['LSR-I'] = predictions_df['LSR_Comb'].astype(float)
predictions_df['LARS'] = predictions_df['LARS'].astype(float)
predictions_df['LARS-I'] = predictions_df['LARS_Comb'].astype(float)
predictions_df['RF'] = predictions_df['RF'].astype(float)
#We replace the non NA values in CNN by integers
for n in range(predictions_df.shape[0]):
    if not pd.isnull(predictions_df.iloc[n]['CNN']):
        predictions_df.at[n, 'CNN'] = float(predictions_df.iloc[n]['CNN'])
    if pd.isnull(predictions_df.iloc[n]['HIVDB']):
        predictions_df.at[n, 'HIVDB'] = 0


# #we calculate the ROC and precision-recall AUC for each drug, fold and model
# roc_auc_dic = {}
# pr_auc_dic = {}
# roc_auc_list = []
# pr_auc_list = []
# for drug in test_drugs:

#     roc_auc_dic[drug] = {}
#     pr_auc_dic[drug] = {}
#     # for model in ['HIVDB', 'LSR', 'LSR_Comb', 'LARS', 'LARS_Comb', 'RF', 'CNN', 'HIVDB+LSR_C+RF', 'HIVDB+LSR_C+LARS_C+RF+CNN']:
#     for model in ['LSR', 'LSR_Comb', 'LARS', 'LARS_Comb', 'RF', 'CNN']:
#         roc_auc_dic[drug][model] = []
#         pr_auc_dic[drug][model] = []
        
#         if model == 'CNN' and drug in INIs:
#             continue
        
#         #we remove the rows with NA in the model
#         predictions_df_noNA = predictions_df[predictions_df[model].notna()]
#         for fold in range(5):
#             print(drug, model, fold)
            
            
#             y_true = predictions_df_noNA[(predictions_df_noNA['Drug']==drug)&(predictions_df_noNA['Fold']==fold)]['True_Label']
#             y_pred = predictions_df_noNA[(predictions_df_noNA['Drug']==drug)&(predictions_df_noNA['Fold']==fold)][model]
#             print(y_true)
#             print(y_pred)
            
#             # TP, TN, FP, FN = 0, 0, 0, 0
#             # for n in range(len(y_true)):
#             #     if y_true.iloc[n]==1 and y_pred.iloc[n]==1:
#             #         TP += 1
#             #     elif y_true.iloc[n]==1 and y_pred.iloc[n]==0:
#             #         FN += 1
#             #     elif y_true.iloc[n]==0 and y_pred.iloc[n]==1:
#             #         FP += 1
#             #     elif y_true.iloc[n]==0 and y_pred.iloc[n]==0:
#             #         TN += 1

#             # fpr, tpr = FP/(FP+TN), TP/(TP+FN)
#             # recall = TP/(TP+FN)
#             # precision = TP/(TP+FP)

#             fpr, tpr, _ = roc_curve(y_true, y_pred)
#             roc_auc = auc(fpr, tpr)
#             roc_auc_dic[drug][model].append(roc_auc)
#             roc_auc_list.append([drug, model, fold, roc_auc])

#             precision, recall, _ = precision_recall_curve(y_true, y_pred)
#             pr_auc = auc(recall, precision)
#             pr_auc_dic[drug][model].append(pr_auc)
#             pr_auc_list.append([drug, model, fold, pr_auc])

# #we store the results in a dataframe
# # roc_auc_df = pd.DataFrame(roc_auc_dic)
# # pr_auc_df = pd.DataFrame(pr_auc_dic)
# roc_auc_df = pd.DataFrame(roc_auc_list, columns=['Drug', 'Model', 'Fold', 'ROC_AUC'])
# pr_auc_df = pd.DataFrame(pr_auc_list, columns=['Drug', 'Model', 'Fold', 'PR_AUC'])
# roc_auc_df.to_csv(Path(DATADIR, 'roc_auc_df_5folds.tsv'), sep='\t', index=False)
# pr_auc_df.to_csv(Path(DATADIR, 'pr_auc_df_5folds.tsv'), sep='\t', index=False)

###We plot the ROC and precision-recall curves for each drug (summing the 5 folds)
test_drugs = INIs
letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
fig, ax = plt.subplots(int(len(test_drugs)/2), 4, figsize = (20,int(5*(len(test_drugs)/2))))

for i, drug in enumerate(test_drugs):

    # for model in ['HIVDB', 'LSR', 'LSR_Comb', 'LARS', 'LARS_Comb', 'RF', 'CNN', 'HIVDB+LSR_C+RF', 'HIVDB+LSR_C+LARS_C+RF+CNN']:
    for model in ['LSR', 'LSR-I', 'LARS', 'LARS-I', 'RF', 'CNN']:
        if model == 'CNN' and drug in INIs:
            continue
        # we remove the rows with NA in the model
        predictions_df_noNA = predictions_df[predictions_df[model].notna()]
        y_true = predictions_df_noNA[predictions_df_noNA['Drug']==drug]['True_Label']
        y_pred = predictions_df_noNA[predictions_df_noNA['Drug']==drug][model]
        print(drug, model)
        print(y_true.shape, y_pred.shape)
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = roc_auc_score(y_true, y_pred)
        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        pr_auc = auc(recall, precision)

        if i % 2 == 0:
            j,k = 0,1
        else:
            j,k = 2,3

        l = int((i+0.01)//2)
        ax[l,j].plot([0, 1], [0, 1], 'k--')
        ax[l,j].set_xlim([0.0, 1.0])
        ax[l,j].set_ylim([0.0, 1.05])
        ax[l,j].set_xlabel('False Positive Rate')
        ax[l,j].set_ylabel('True Positive Rate')
        ax[l,j].set_title(f'{drug} ROC curve')
        ax[l,j].plot(fpr, tpr, label=f'{model} (AUC = {roc_auc:.3f})')
        ax[l,j].legend(loc='lower right', fontsize='small')
        ax[l,j].text(-0.25,1,letters[i], fontsize = 15)
        ax[l,k].plot([0, 1], [1, 0], 'k--')
        ax[l,k].set_xlim([0.0, 1.0])
        ax[l,k].set_ylim([0.0, 1.05])
        ax[l,k].set_xlabel('Recall')
        ax[l,k].set_ylabel('Precision')
        ax[l,k].set_title(f'{drug} Precision-Recall curve')
        ax[l,k].plot(recall, precision, label=f'{model} (AUC = {pr_auc:.3f})')
        ax[l,k].legend(loc='lower left', fontsize='small')

fig.suptitle("AUROC and AUPRC for the different INI drugs", y=0.99, fontsize = 15)
fig.tight_layout(pad=2.0)
plt.savefig(Path(DATADIR, f'roc_pr_curves/INI_roc_pr_curves_5folds_combined.png'))
plt.close()
