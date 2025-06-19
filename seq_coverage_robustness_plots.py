###We will plot bargraphs showing accuracy, balanced accuracy, F1-score and MCC for each drug class and for different number of missing major DRM positions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def perf_metrics(drug, dataset, mutations_df, mut_class = "all"):
    '''Calculates the performance metrics for the ensemble predictions.
    Based on % of uncovered positions.

    drug: drug name i.e AZT
    dataset: drug class i.e NRTI
    mutations_df: dataframe with SeqID, CompMutList_no_minor and fold columns
    '''
    if dataset == 'NNRTI':
        major_positions = [100, 101, 103, 106, 181, 188, 190, 230]
        all_positions = [90, 98, 100, 101, 103, 106, 108, 138, 179, 181, 188, 190, 221, 225, 227, 230, 234, 236, 238, 318, 348]
    elif dataset == 'NRTI':
        major_positions = [41, 65, 70, 74, 75, 151, 184, 210, 215]
        all_positions = [41, 62, 65, 67, 68, 69, 70, 74, 75, 77, 115, 116, 151, 184, 210, 215, 219]

    elif dataset == 'INI':
        major_positions = [66, 92, 118, 143, 148, 155, 263]
        all_positions = [51, 66, 74, 75, 92, 95, 97, 118, 121, 122, 138, 140, 143, 145, 146, 147, 148, 151, 153, 155, 157, 163, 230, 232, 263]
    elif dataset == 'PI':
        major_positions = [30, 32, 47, 48, 50, 54, 76, 82, 84, 88]
        all_positions = [10, 20, 24, 32, 33, 46, 47, 48, 50, 53, 54, 73, 74, 76, 82, 83, 84, 88, 89, 90]
    else:
        raise ValueError(f'Unknown dataset: {dataset}')
    
    #we define the minor positions as the ones that are not in the major positions
    minor_positions = [i for i in all_positions if i not in major_positions]

    metrics = []
    real_df = pd.read_csv(f"datasets/subsampled_datasets/{dataset}_subsampled_dataset.tsv", sep='\t')
    if mut_class == "all":
        print('Calculating metrics for all mutations')
        real_df['n_missing_muts'] =  real_df['keep_seq_positions'].apply(lambda x: [] if pd.isnull(x) else [y for y in x.split(', ') if int(y) in all_positions])
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: len(x)/len(all_positions))
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: 1-x)
    elif mut_class == "major":
        print('Calculating metrics for major mutations')
        real_df['n_missing_muts'] =  real_df['keep_seq_positions'].apply(lambda x: [] if pd.isnull(x) else [y for y in x.split(', ') if int(y) in major_positions])
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: len(x)/len(major_positions))
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: 1-x)
    elif mut_class == "minor":
        print('Calculating metrics for minor mutations')
        real_df['n_missing_muts'] =  real_df['keep_seq_positions'].apply(lambda x: [] if pd.isnull(x) else [y for y in x.split(', ') if int(y) in minor_positions])
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: len(x)/len(minor_positions))
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: 1-x)
    else:
        raise ValueError(f'Unknown mutation class: {mut_class}')
    

    for m in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
        # print(f'Calculating metrics for {m} missing mutations')
        real_m_df = real_df[real_df['n_missing_muts'].between(m, m+0.1, inclusive='left')]
        #we keep only rows where SeqID and prc equal than the ones in real_m_df
        mutations_prc = pd.DataFrame()
        for prc in percentage_sequence:
            real_m_df_prc = real_m_df[real_m_df['seq_percent'] == prc]
            mutations_prc_n = mutations_df[(mutations_df['SeqID'].isin(real_m_df_prc['SeqID'])) & (mutations_df['seq_percent'] == prc)]
            mutations_prc = pd.concat([mutations_prc, mutations_prc_n], axis=0)

        mutations_m = mutations_prc
            
        for fold in range(5):
            # print(f'Calculating metrics for fold {fold}')
            TP, TN, FP, FN = 0, 0, 0, 0
            mutations_df_fold = mutations_m[mutations_m['fold'] == fold]
            length = mutations_df_fold.shape[0]
            for n in range(mutations_df_fold.shape[0]):
                seqID, prc = mutations_df_fold.iloc[n]['SeqID'], mutations_df_fold.iloc[n]['seq_percent']
                true = real_df[real_df['SeqID'] == seqID][f'{drug}_res'].values[0]
                prediction = mutations_df_fold[(mutations_df_fold['SeqID'] == seqID)&(mutations_df_fold['seq_percent']==prc)]['Ensemble'].values[0]
            
                if true == 'Resistant' and prediction == 'Resistant':
                    TP += 1
                elif true == 'Susceptible' and prediction == 'Susceptible':
                    TN += 1
                elif true == 'Resistant' and prediction == 'Susceptible':
                    FN += 1
                elif true == 'Susceptible' and prediction == 'Resistant':
                    FP += 1

            if TP + TN + FP + FN == 0:
                print(f'No predictions for fold {fold} and missing mutations {m}')
                continue

            accuracy = (TP + TN) / (TP + TN + FP + FN)
            try:
                sensitivity = TP / (TP + FN)
            except ZeroDivisionError:
                sensitivity = 0
            try:
                specificity = TN / (TN + FP)
            except ZeroDivisionError:
                specificity = 0
            balanced_accuracy = (sensitivity + specificity) / 2
            try:
                precision = TP / (TP + FP)
            except ZeroDivisionError:
                precision = 0
            recall = sensitivity
            try:
                F1 = 2 * ((precision * recall) / (precision + recall))
            except ZeroDivisionError:
                F1 = 0
            try:
                MCC = ((TP * TN) - (FP * FN)) / ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))**(1/2)
            except ZeroDivisionError:
                MCC = 0
            metrics.append([m, fold, accuracy, balanced_accuracy, F1, MCC, length])
    
    metrics_df = pd.DataFrame(metrics, columns=['Miss_muts', 'Fold', 'Accuracy', 'Balanced_Accuracy', 'F1', 'MCC', 'Length'])
    return metrics_df

def perf_metrics_number(drug, dataset, mutations_df, mut_class = "all"):
    '''Calculates the performance metrics for the ensemble predictions.
    Based on total number of uncovered positions.

    drug: drug name i.e AZT
    dataset: drug class i.e NRTI
    mutations_df: dataframe with SeqID, CompMutList_no_minor and fold columns
    '''
    if dataset == 'NNRTI':
        major_positions = [100, 101, 103, 106, 181, 188, 190, 230]
        all_positions = [90, 98, 100, 101, 103, 106, 108, 138, 179, 181, 188, 190, 221, 225, 227, 230, 234, 236, 238, 318, 348]
    elif dataset == 'NRTI':
        major_positions = [41, 65, 70, 74, 75, 151, 184, 210, 215]
        all_positions = [41, 62, 65, 67, 68, 69, 70, 74, 75, 77, 115, 116, 151, 184, 210, 215, 219]

    elif dataset == 'INI':
        major_positions = [66, 92, 118, 143, 148, 155, 263]
        all_positions = [51, 66, 74, 75, 92, 95, 97, 118, 121, 122, 138, 140, 143, 145, 146, 147, 148, 151, 153, 155, 157, 163, 230, 232, 263]
    elif dataset == 'PI':
        major_positions = [30, 32, 47, 48, 50, 54, 76, 82, 84, 88]
        all_positions = [10, 20, 24, 32, 33, 46, 47, 48, 50, 53, 54, 73, 74, 76, 82, 83, 84, 88, 89, 90]
    else:
        raise ValueError(f'Unknown dataset: {dataset}')
    
    #we define the minor positions as the ones that are not in the major positions
    minor_positions = [i for i in all_positions if i not in major_positions]

    metrics = []
    real_df = pd.read_csv(f"datasets/subsampled_datasets/{dataset}_subsampled_dataset.tsv", sep='\t')
    if mut_class == "all":
        print('Calculating metrics for all mutations')
        real_df['n_missing_muts'] =  real_df['keep_seq_positions'].apply(lambda x: [] if pd.isnull(x) else [y for y in x.split(', ') if int(y) in all_positions])
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: len(all_positions)-len(x))

    elif mut_class == "major":
        print('Calculating metrics for major mutations')
        real_df['n_missing_muts'] =  real_df['keep_seq_positions'].apply(lambda x: [] if pd.isnull(x) else [y for y in x.split(', ') if int(y) in major_positions])
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: len(major_positions)-len(x))
    elif mut_class == "minor":
        print('Calculating metrics for minor mutations')
        real_df['n_missing_muts'] =  real_df['keep_seq_positions'].apply(lambda x: [] if pd.isnull(x) else [y for y in x.split(', ') if int(y) in minor_positions])
        real_df['n_missing_muts'] = real_df['n_missing_muts'].apply(lambda x: len(minor_positions)-len(x))
    else:
        raise ValueError(f'Unknown mutation class: {mut_class}')
    

    for m in range(len(all_positions) + 1):
        # print(f'Calculating metrics for {m} missing mutations')
        real_m_df = real_df[real_df['n_missing_muts']==m]
        #we keep only rows where SeqID and prc equal than the ones in real_m_df
        mutations_prc = pd.DataFrame()
        for prc in percentage_sequence:
            real_m_df_prc = real_m_df[real_m_df['seq_percent'] == prc]
            mutations_prc_n = mutations_df[(mutations_df['SeqID'].isin(real_m_df_prc['SeqID'])) & (mutations_df['seq_percent'] == prc)]
            mutations_prc = pd.concat([mutations_prc, mutations_prc_n], axis=0)

        mutations_m = mutations_prc
            
        for fold in range(5):
            # print(f'Calculating metrics for fold {fold}')
            TP, TN, FP, FN = 0, 0, 0, 0
            mutations_df_fold = mutations_m[mutations_m['fold'] == fold]
            length = mutations_df_fold.shape[0]
            for n in range(mutations_df_fold.shape[0]):
                seqID, prc = mutations_df_fold.iloc[n]['SeqID'], mutations_df_fold.iloc[n]['seq_percent']
                true = real_df[real_df['SeqID'] == seqID][f'{drug}_res'].values[0]
                prediction = mutations_df_fold[(mutations_df_fold['SeqID'] == seqID)&(mutations_df_fold['seq_percent']==prc)]['Ensemble'].values[0]
            
                if true == 'Resistant' and prediction == 'Resistant':
                    TP += 1
                elif true == 'Susceptible' and prediction == 'Susceptible':
                    TN += 1
                elif true == 'Resistant' and prediction == 'Susceptible':
                    FN += 1
                elif true == 'Susceptible' and prediction == 'Resistant':
                    FP += 1

            if TP + TN + FP + FN == 0:
                print(f'No predictions for fold {fold} and missing mutations {m}')
                continue

            accuracy = (TP + TN) / (TP + TN + FP + FN)
            try:
                sensitivity = TP / (TP + FN)
            except ZeroDivisionError:
                sensitivity = 0
            try:
                specificity = TN / (TN + FP)
            except ZeroDivisionError:
                specificity = 0
            balanced_accuracy = (sensitivity + specificity) / 2
            try:
                precision = TP / (TP + FP)
            except ZeroDivisionError:
                precision = 0
            recall = sensitivity
            try:
                F1 = 2 * ((precision * recall) / (precision + recall))
            except ZeroDivisionError:
                F1 = 0
            try:
                MCC = ((TP * TN) - (FP * FN)) / ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))**(1/2)
            except ZeroDivisionError:
                MCC = 0
            metrics.append([m, fold, accuracy, balanced_accuracy, F1, MCC, length])
    
    metrics_df = pd.DataFrame(metrics, columns=['Miss_muts', 'Fold', 'Accuracy', 'Balanced_Accuracy', 'F1', 'MCC', 'Length'])
    return metrics_df



test_drugs = ['RAL', 'EVG', 'DTG', 'BIC', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
INIs=['RAL', 'EVG', 'DTG', 'BIC']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
drug_classes = [INIs, NNRTIs, NRTIs, PIs]
drug_classes_names = ['INI', 'NNRTI', 'NRTI', 'PI']

pi_major_positions = [30, 32, 47, 48, 50, 54, 76, 82, 84, 88]
nrti_major_positions = [41, 65, 70, 74, 75, 151, 184, 210, 215]
nnrti_major_positions = [100, 101, 103, 106, 181, 188, 190, 230]
ini_major_positions = [66, 92, 118, 143, 148, 155, 263]

percentage_sequence = [0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.0]


ensemble_results_all, ensemble_results_major, ensemble_results_minor = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
ensemble_results_total_all, ensemble_results_total_major, ensemble_results_total_minor = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
for drug in test_drugs:
    if drug in PIs:
        dataset = "PI"
    elif drug in NNRTIs:
        dataset = "NNRTI"
    elif drug in NRTIs:
        dataset = "NRTI"
    elif drug in INIs:
        dataset = "INI"

    ensemble_results = pd.read_csv(f'method_predictions/ensemble_subsampled_predictions/{dataset}/{drug}_ensemble_predictions.tsv', sep='\t')
    
    #% and total of all positions
    metrics_total_all = perf_metrics_number(drug, dataset, ensemble_results)
    # metrics_all, metrics_total_all = perf_metrics(drug, dataset, ensemble_results), perf_metrics_number(drug, dataset, ensemble_results)
    # metrics_all['Drug'] = drug
    metrics_total_all['Drug'] = drug
    # ensemble_results_all = pd.concat([ensemble_results_all, metrics_all], ignore_index=True)
    ensemble_results_total_all = pd.concat([ensemble_results_total_all, metrics_total_all], ignore_index=True)
    # print(metrics_all, metrics_total_all)

    #% of major positions
    metrics_total_major = perf_metrics_number(drug, dataset, ensemble_results, mut_class="major")
    # metrics_major, metrics_total_major = perf_metrics(drug, dataset, ensemble_results, mut_class="major"), perf_metrics_number(drug, dataset, ensemble_results, mut_class="major")
    # metrics_major['Drug'] = drug
    metrics_total_major['Drug'] = drug
    # ensemble_results_major = pd.concat([ensemble_results_major, metrics_major], ignore_index=True)
    ensemble_results_total_major = pd.concat([ensemble_results_total_major, metrics_total_major], ignore_index=True)
    # print(metrics_major, metrics_total_major)

    #% of minor positions
    metrics_total_minor = perf_metrics_number(drug, dataset, ensemble_results, mut_class="minor")
    # metrics_minor, metrics_total_minor = perf_metrics(drug, dataset, ensemble_results, mut_class="minor"), perf_metrics_number(drug, dataset, ensemble_results, mut_class="minor")
    # metrics_minor['Drug'] = drug
    metrics_total_minor['Drug'] = drug
    # ensemble_results_minor = pd.concat([ensemble_results_minor, metrics_minor], ignore_index=True)
    ensemble_results_total_minor = pd.concat([ensemble_results_total_minor, metrics_total_minor], ignore_index=True)
    # print(metrics_minor, metrics_total_minor)

    print(f'Drug {drug} metrics done')

# ensemble_performance_seq_percent_data = ensemble_results_all
# ensemble_performance_seq_percent_data_major = ensemble_results_major
# ensemble_performance_seq_percent_data_minor = ensemble_results_minor

# #we do not consider BIC fold 1
# ensemble_performance_seq_percent_data = ensemble_performance_seq_percent_data.drop(ensemble_performance_seq_percent_data[(ensemble_performance_seq_percent_data['Drug'] == 'BIC') & (ensemble_performance_seq_percent_data['Fold'] == 1)].index)
# ensemble_performance_seq_percent_data_major = ensemble_performance_seq_percent_data_major.drop(ensemble_performance_seq_percent_data_major[(ensemble_performance_seq_percent_data_major['Drug'] == 'BIC') & (ensemble_performance_seq_percent_data_major['Fold'] == 1)].index)
# ensemble_performance_seq_percent_data_minor = ensemble_performance_seq_percent_data_minor.drop(ensemble_performance_seq_percent_data_minor[(ensemble_performance_seq_percent_data_minor['Drug'] == 'BIC') & (ensemble_performance_seq_percent_data_minor['Fold'] == 1)].index)


# ##Step-line graph for all four measurements related to sequence coverage % (Figure 3 and Supplemetary Figures 9-11)
# for metric in ['Accuracy', 'Balanced_Accuracy', 'F1', 'MCC']:
#     fig, axs = plt.subplots(4, 1, figsize=(20, 20))
    
#     for i, cl in enumerate(drug_classes):

#         ensemble_performance_seq_percent_data_drug = ensemble_performance_seq_percent_data[ensemble_performance_seq_percent_data['Drug'].isin(cl)]
#         ensemble_performance_seq_percent_data_major_drug = ensemble_performance_seq_percent_data_major[ensemble_performance_seq_percent_data_major['Drug'].isin(cl)]
#         ensemble_performance_seq_percent_data_minor_drug = ensemble_performance_seq_percent_data_minor[ensemble_performance_seq_percent_data_minor['Drug'].isin(cl)]

#         if metric == "Balanced_Accuracy":
#             axs[i].axhline(y=0.7, color='black', linestyle='--', alpha= 0.3)

#         if cl == NNRTIs:
#             major_positions = nnrti_major_positions
#         elif cl == NRTIs:
#             major_positions = nrti_major_positions
#         elif cl == INIs:
#             major_positions = ini_major_positions
#         elif cl == PIs:
#             major_positions = pi_major_positions
        
#         x = [0.5*m+1 for m in range(11)]
#         y, y_major, y_minor = [], [], []

#         for m,p in enumerate([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]):
#             all_drms, major, minor = [], [], []

#             ensemble_m_mutations_df = ensemble_performance_seq_percent_data_drug[ensemble_performance_seq_percent_data_drug['Miss_muts'] == p]
#             ensemble_m_mutations_major_df = ensemble_performance_seq_percent_data_major_drug[ensemble_performance_seq_percent_data_major_drug['Miss_muts'] == p]
#             ensemble_m_mutations_minor_df = ensemble_performance_seq_percent_data_minor_drug[ensemble_performance_seq_percent_data_minor_drug['Miss_muts'] == p]

#             for k in range(5):
#                 all_drms.append(ensemble_m_mutations_df[ensemble_m_mutations_df['Fold'] == k][metric].mean())
#                 major.append(ensemble_m_mutations_major_df[ensemble_m_mutations_major_df['Fold'] == k][metric].mean())
#                 minor.append(ensemble_m_mutations_minor_df[ensemble_m_mutations_minor_df['Fold'] == k][metric].mean())

#             #we plot the error bars for all positions
#             if m < 10:
#                 axs[i].errorbar(0.5*m+1.25, np.mean(all_drms), yerr=np.std(all_drms), color='blue', alpha = 0.7)
#                 axs[i].errorbar(0.5*m+1.25, np.mean(major), yerr=np.std(major), color='orange', alpha = 0.7)
#                 axs[i].errorbar(0.5*m+1.25, np.mean(minor), yerr=np.std(minor), color='green', alpha = 0.7)

#             y.append(np.mean(all_drms)), y_major.append(np.mean(major)), y_minor.append(np.mean(minor))
           
#         axs[i].step(x, y, where = 'post', label='All DRM', color='blue', alpha = 0.7)
#         axs[i].step(x, y_major, where = 'post', label='Major DRM', color='orange', alpha = 0.7)
#         axs[i].step(x, y_minor, where = 'post', label='Minor DRM', color='green', alpha = 0.7)
        
#     axs[3].legend(fontsize=15, loc='lower center', bbox_to_anchor=(0.5, -0.35), ncol=3)

#     axs[0].set_title(f'INIs', fontsize=20)
#     axs[1].set_title(f'NNRTIs', fontsize=20)
#     axs[2].set_title(f'NRTIs', fontsize=20)
#     axs[3].set_title(f'PIs', fontsize=20)
    
#     for n in range(4):
#         axs[n].set_xticks([0.5*i+1 for i in range(11)])
#         axs[n].set_xticklabels([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
#         axs[n].set_ylim(0, 1)
#         axs[n].set_ylabel(metric.replace("_"," ").replace("F1", "F1 score") + " (5 folds)", fontsize=15)
#         if metric == "Accuracy" or metric == "Balanced_Accuracy":
#             axs[n].set_ylim(0.4,1)
#         axs[n].tick_params(axis='x', labelsize=15)
#         axs[n].tick_params(axis='y', labelsize=15)

#     axs[3].set_xlabel('% of uncovered DRM positions', fontsize=20)
#     fig.suptitle(f"{metric.replace("_"," ").replace("F1", "F1 score")} for different DRM positions coverage percentages", size=25, y=0.95)

#     plt.savefig(f"figures/{metric.replace("F1", "f1score").lower()}_robustness_percent_plot.png")

##Now we look at the total number of uncovered positions
ensemble_performance_seq_percent_data = ensemble_results_total_all
ensemble_performance_seq_percent_data_major = ensemble_results_total_major
ensemble_performance_seq_percent_data_minor = ensemble_results_total_minor

#we do not consider BIC fold 1
ensemble_performance_seq_percent_data = ensemble_performance_seq_percent_data.drop(ensemble_performance_seq_percent_data[(ensemble_performance_seq_percent_data['Drug'] == 'BIC') & (ensemble_performance_seq_percent_data['Fold'] == 1)].index)
ensemble_performance_seq_percent_data_major = ensemble_performance_seq_percent_data_major.drop(ensemble_performance_seq_percent_data_major[(ensemble_performance_seq_percent_data_major['Drug'] == 'BIC') & (ensemble_performance_seq_percent_data_major['Fold'] == 1)].index)
ensemble_performance_seq_percent_data_minor = ensemble_performance_seq_percent_data_minor.drop(ensemble_performance_seq_percent_data_minor[(ensemble_performance_seq_percent_data_minor['Drug'] == 'BIC') & (ensemble_performance_seq_percent_data_minor['Fold'] == 1)].index)


##Line plot for all four measurements related to sequence coverage in total number of positions
##(Supplemetary Figures 12-15)

for metric in ['Accuracy', 'Balanced_Accuracy', 'F1', 'MCC']:
    fig, axs = plt.subplots(4, 1, figsize=(20, 20))

    for i, cl in enumerate(drug_classes):

        ensemble_performance_seq_percent_data_drug = ensemble_performance_seq_percent_data[ensemble_performance_seq_percent_data['Drug'].isin(cl)]
        ensemble_performance_seq_percent_data_major_drug = ensemble_performance_seq_percent_data_major[ensemble_performance_seq_percent_data_major['Drug'].isin(cl)]
        ensemble_performance_seq_percent_data_minor_drug = ensemble_performance_seq_percent_data_minor[ensemble_performance_seq_percent_data_minor['Drug'].isin(cl)]
        if metric == 'Balanced_Accuracy':
            axs[i].axhline(y=0.7, color='black', linestyle='--', alpha= 0.3)

        if cl == NNRTIs:
            major_positions = nnrti_major_positions
            max_len = 21
        elif cl == NRTIs:
            major_positions = nrti_major_positions
            max_len = 17
        elif cl == INIs:
            major_positions = ini_major_positions
            max_len = 25
        elif cl == PIs:
            major_positions = pi_major_positions
            max_len = 20
        
        x = [0.5*m+1 for m in range(max_len+1)]
        y, y_major, y_minor = [], [], []
        for m in range(max_len+1):
            all_drms, major, minor = [], [], []
            ensemble_m_mutations_df = ensemble_performance_seq_percent_data_drug[ensemble_performance_seq_percent_data_drug['Miss_muts'] == m]
            ensemble_m_mutations_major_df = ensemble_performance_seq_percent_data_major_drug[ensemble_performance_seq_percent_data_major_drug['Miss_muts'] == m]
            ensemble_m_mutations_minor_df = ensemble_performance_seq_percent_data_minor_drug[ensemble_performance_seq_percent_data_minor_drug['Miss_muts'] == m]

            for k in range(5):
                all_drms.append(ensemble_m_mutations_df[ensemble_m_mutations_df['Fold'] == k][metric].mean())
                major.append(ensemble_m_mutations_major_df[ensemble_m_mutations_major_df['Fold'] == k][metric].mean())
                minor.append(ensemble_m_mutations_minor_df[ensemble_m_mutations_minor_df['Fold'] == k][metric].mean())

            # #we plot the error bars for all positions
            axs[i].errorbar(0.5*m+1, np.mean(all_drms), yerr=np.std(all_drms), color='blue', alpha = 0.7)
            axs[i].errorbar(0.5*m+1, np.mean(major), yerr=np.std(major), color='orange', alpha = 0.7)
            axs[i].errorbar(0.5*m+1, np.mean(minor), yerr=np.std(minor), color='green', alpha = 0.7)

            y.append(np.mean(all_drms)), y_major.append(np.mean(major)), y_minor.append(np.mean(minor))

        axs[i].plot(x, y, label='All DRM', color='blue', alpha = 0.7, marker='o')
        axs[i].plot(x, y_major, label='Major DRM', color='orange', alpha = 0.7, marker='o')
        axs[i].plot(x, y_minor, label='Minor DRM', color='green', alpha = 0.7, marker='o')

        axs[i].set_xticks([0.5*i+1 for i in range(max_len+1)])
        axs[i].set_xticklabels(range(0, max_len+1))
        
    axs[3].legend(fontsize=15, loc='lower center', bbox_to_anchor=(0.5, -0.35), ncol=3)

    axs[0].set_title(f'INIs', fontsize=20)
    axs[1].set_title(f'NNRTIs', fontsize=20)
    axs[2].set_title(f'NRTIs', fontsize=20)
    axs[3].set_title(f'PIs', fontsize=20)
    
    for n in range(4):
        axs[n].set_ylim(0, 1)
        if metric == "Accuracy" or metric == "Balanced_Accuracy":
            axs[n].set_ylim(0.4, 1)
        axs[n].set_ylabel(metric.replace("_"," ").replace("F1", "F1 score") + " (5 folds)", fontsize=15)
        axs[n].tick_params(axis='x', labelsize=15)
        axs[n].tick_params(axis='y', labelsize=15)

    axs[3].set_xlabel('Number of uncovered DRM positions', fontsize=20)
    fig.suptitle(f"{metric.replace("_"," ").replace("F1", "F1 score")} for different numbers of uncovered DRM positions", size=25, y=0.95)

    plt.savefig(f"figures/{metric.replace("F1", "f1score").lower()}_robustness_total_number_plot.png")
