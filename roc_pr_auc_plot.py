###We will plot the ROC and the precision-recall curve for each model
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc

test_drugs = ['RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
test_drugs_no_INIs = ['EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
INIs=['RAL', 'EVG', 'DTG', 'BIC']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] #RPV test set is empty
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
dataset_names = ["INI", "PI", "NNRTI", "NRTI"]


predictions_df = pd.read_csv("method_predictions/roc_pr_auc_data.tsv", sep='\t')


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


###We plot the ROC and precision-recall curves for each drug (summing the 5 folds)
letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']

for m,test_drugs in enumerate([INIs, PIs, NNRTIs, NRTIs]):
    fig, ax = plt.subplots(int(len(test_drugs)/2), 4, figsize = (20,int(5*(len(test_drugs)/2))))
    dataset_name = dataset_names[m]

    for i, drug in enumerate(test_drugs):

        for model in ['LSR', 'LSR-I', 'LARS', 'LARS-I', 'RF', 'CNN']:
            if model == 'CNN' and drug in INIs:
                continue
            # we remove the rows with NA in the model
            predictions_df_noNA = predictions_df[predictions_df[model].notna()]
            y_true = predictions_df_noNA[predictions_df_noNA['Drug']==drug]['True_Label']
            y_pred = predictions_df_noNA[predictions_df_noNA['Drug']==drug][model]
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
    plt.savefig(f'figures/roc_pr_curves/{dataset_name}_roc_pr_curves_5folds_combined.png', dpi=300)
    plt.close()
