###We create barplots with error bars to compare the different methods performance
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

#data for 5 folds
hivdb = pd.read_csv("method_predictions/HIVDB_performance_5folds.tsv", sep="\t")
lsr = pd.read_csv("method_predictions/LSR_performance_5folds.tsv", sep="\t")
lsr_i = pd.read_csv("method_predictions/LSR_I_performance_5folds.tsv", sep="\t")
lars = pd.read_csv("method_predictions/LARS_performance_5folds.tsv", sep="\t")
lars_i = pd.read_csv("method_predictions/LARS_I_performance_5folds.tsv", sep="\t")
random_forest = pd.read_csv("method_predictions/random_forest_performance_5folds.tsv", sep="\t")
cnn = pd.read_csv("method_predictions/cnn_performance_5folds.tsv", sep="\t")
small_ensemble = pd.read_csv("method_predictions/small_ensemble_performance_5folds.tsv", sep="\t")
large_ensemble = pd.read_csv("method_predictions/large_ensemble_performance_5folds.tsv", sep="\t")

##We remove BIC drug fold 1 row from the result calculations as it is completely unbalanced
hivdb = hivdb.drop(hivdb[(hivdb["Drug"] == "BIC") & (hivdb["Fold"] == 1)].index)
lsr = lsr.drop(lsr[(lsr["Drug"] == "BIC") & (lsr["Fold"] == 1)].index)
lsr_i = lsr_i.drop(lsr_i[(lsr_i["Drug"] == "BIC") & (lsr_i["Fold"] == 1)].index)
lars = lars.drop(lars[(lars["Drug"] == "BIC") & (lars["Fold"] == 1)].index)
lars_i = lars_i.drop(lars_i[(lars_i["Drug"] == "BIC") & (lars_i["Fold"] == 1)].index)
random_forest = random_forest.drop(random_forest[(random_forest["Drug"] == "BIC") & (random_forest["Fold"] == 1)].index)
cnn = cnn.drop(cnn[(cnn["Drug"] == "BIC") & (cnn["Fold"] == 1)].index)
small_ensemble = small_ensemble.drop(small_ensemble[(small_ensemble["Drug"] == "BIC") & (small_ensemble["Fold"] == 1)].index)
large_ensemble = large_ensemble.drop(large_ensemble[(large_ensemble["Drug"] == "BIC") & (large_ensemble["Fold"] == 1)].index)

#we define lists with drugs and drug classes labels
test_drugs = ["RAL", "EVG", "DTG", "BIC", "EFV", "NVP", "ETR", "RPV", "3TC", "ABC", "AZT", "D4T", "DDI", "TDF", "FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"]
test_drugs_no_INIs = ["EFV", "NVP", "ETR", "RPV", "3TC", "ABC", "AZT", "D4T", "DDI", "TDF", "FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"]
INIs=["RAL", "EVG", "DTG", "BIC"]
PIs=["FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"]
NNRTIs=["EFV", "NVP", "ETR", "RPV"]
NRTIs=["3TC", "ABC", "AZT", "D4T", "DDI", "TDF"]
classes_drugs = [INIs, NNRTIs, NRTIs, PIs]
classes_names = ["INIs", "NNRTIs", "NRTIs", "PIs"]

methods = ["HIVDB", "LSR", "LSR-I", "LARS", "LARS-I", "RandomForest",  "CNN", "HIVDB+LSR-I+RF", "HIVDB+LSR-I+LARS-I+RF+CNN"]

performance_list = []
for cl in classes_drugs:
    #we select the drugs in the defined class
    hivdb_class = hivdb[hivdb["Drug"].isin(cl)]
    lsr_class = lsr[lsr["Drug"].isin(cl)]
    lsr_i_class = lsr_i[lsr_i["Drug"].isin(cl)]
    lars_class = lars[lars["Drug"].isin(cl)]
    lars_i_class = lars_i[lars_i["Drug"].isin(cl)]
    random_forest_class = random_forest[random_forest["Drug"].isin(cl)]
    cnn_class = cnn[cnn["Drug"].isin(cl)]
    small_ensemble_class = small_ensemble[small_ensemble["Drug"].isin(cl)]
    large_ensemble_class = large_ensemble[large_ensemble["Drug"].isin(cl)]

    #we retrieve each metric for all folds
    for metric in ["Accuracy", "Balanced_accuracy", "F1_score", "MCC"]:
        HIVDB_folds, lsr_folds, lsr_i_folds, lars_folds, lars_i_folds, rf_folds, cnn_folds, ensemble_folds, ensemble_5_folds = [], [], [], [], [], [], [], [], []
        for k in range(5):
            HIVDB_folds.append(hivdb_class[hivdb_class["Fold"] == k][metric].mean())
            lsr_folds.append(lsr_class[lsr_class["Fold"] == k][metric].mean())
            lsr_i_folds.append(lsr_i_class[lsr_i_class["Fold"] == k][metric].mean())
            lars_folds.append(lars_class[lars_class["Fold"] == k][metric].mean())
            lars_i_folds.append(lars_i_class[lars_i_class["Fold"] == k][metric].mean())
            rf_folds.append(random_forest_class[random_forest_class["Fold"] == k][metric].mean())
            cnn_folds.append(cnn_class[cnn_class["Fold"] == k][metric].mean())
            ensemble_folds.append(small_ensemble_class[small_ensemble_class["Fold"] == k][metric].mean())
            ensemble_5_folds.append(large_ensemble_class[large_ensemble_class["Fold"] == k][metric].mean())

        if cl == INIs: #There are no data available for INIs in the CNN implementation
            CNN_mean_class, CNN_std_class = 0,0
        else:
            CNN_mean_class, CNN_std_class = np.mean(cnn_folds), np.std(cnn_folds)
        
        #we calculate the mean and the standard deviation for the five folds
        HIVDB_mean_class, HIVDB_std_class = np.mean(HIVDB_folds), np.std(HIVDB_folds)
        LSR_mean_class, LSR_std_class = np.mean(lsr_folds), np.std(lsr_folds)
        LSR_i_mean_class, LSR_i_std_class = np.mean(lsr_i_folds), np.std(lsr_i_folds)
        LARS_mean_class, LARS_std_class = np.mean(lars_folds), np.std(lars_folds)
        LARS_i_mean_class, LARS_i_std_class = np.mean(lars_i_folds), np.std(lars_i_folds)
        rf_mean_class, rf_std_class = np.mean(rf_folds), np.std(rf_folds)
        small_ensemble_mean_class, small_ensemble_std_class = np.mean(ensemble_folds), np.std(ensemble_folds)
        large_ensemble_mean_class, large_ensemble_std_class = np.mean(ensemble_5_folds), np.std(ensemble_5_folds)

        class_name = classes_names[classes_drugs.index(cl)]
        performance_list.append([metric, class_name, HIVDB_mean_class, HIVDB_std_class, LSR_mean_class, LSR_std_class, LSR_i_mean_class, LSR_i_std_class, LARS_mean_class, LARS_std_class, LARS_i_mean_class, LARS_i_std_class, rf_mean_class, rf_std_class, CNN_mean_class, CNN_std_class, small_ensemble_mean_class, small_ensemble_std_class, large_ensemble_mean_class, large_ensemble_std_class])
                             
performance_table = pd.DataFrame(performance_list, columns=["Metric", "Drug_Class", "HIVDB_mean", "HIVDB_std", "LSR_mean", "LSR_std", "LSR_I_mean", "LSR_I_std", "LARS_mean", "LARS_std", "LARS_I_mean", "LARS_I_std", "rf_mean", "rf_std", "CNN_mean", "CNN_std", "small_ensemble_mean", "small_ensemble_std", "large_ensemble_mean", "large_ensemble_std"])

##PERFORMANCE BARPLOTS (Figure 2)

for metric in ["Accuracy", "Balanced_accuracy", "F1_score", "MCC"]:

    #we filter the performance table for the metric we want to plot
    performance_metric_table = performance_table[performance_table["Metric"] == metric]
    ini_table = performance_metric_table[performance_metric_table["Drug_Class"] == "INIs"]
    nnrti_table = performance_metric_table[performance_metric_table["Drug_Class"] == "NNRTIs"]
    nrti_table = performance_metric_table[performance_metric_table["Drug_Class"] == "NRTIs"]
    pi_table = performance_metric_table[performance_metric_table["Drug_Class"] == "PIs"]

    fig, ax = plt.subplots(figsize=(10, 6))

    barWidth = 0.15

    ###INIs (this is the only one without CNN and large ensemble then)
    ax.bar(0.15, ini_table["HIVDB_mean"], yerr=ini_table["HIVDB_std"], color="b", alpha = 0.7,  width=barWidth, edgecolor="black", label="HIVDB")
    ax.bar(0.3, ini_table["LSR_mean"], yerr=ini_table["LSR_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR")
    ax.bar(0.45, ini_table["LSR_I_mean"], yerr=ini_table["LSR_I_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR_I", hatch="//")
    ax.bar(0.6, ini_table["LARS_mean"], yerr=ini_table["LARS_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS")
    ax.bar(0.75, ini_table["LARS_I_mean"], yerr=ini_table["LARS_I_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS_I", hatch="//")
    ax.bar(0.9, ini_table["rf_mean"], yerr=ini_table["rf_std"], color="purple", alpha = 0.7, width=barWidth, edgecolor="black", label="rf")
    ax.bar(1.05, ini_table["small_ensemble_mean"], yerr=ini_table["small_ensemble_std"], color="orange", alpha = 0.7, width=barWidth, edgecolor="black", label="small_ensemble")

    ###NNRTIs
    ax.bar(1.85, nnrti_table["HIVDB_mean"], yerr=nnrti_table["HIVDB_std"], color="b", alpha = 0.7,  width=barWidth, edgecolor="black", label="HIVDB")
    ax.bar(2.00, nnrti_table["LSR_mean"], yerr=nnrti_table["LSR_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR")
    ax.bar(2.15, nnrti_table["LSR_I_mean"], yerr=nnrti_table["LSR_I_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR_I", hatch="//")
    ax.bar(2.30, nnrti_table["LARS_mean"], yerr=nnrti_table["LARS_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS")
    ax.bar(2.45, nnrti_table["LARS_I_mean"], yerr=nnrti_table["LARS_I_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS_I", hatch="//")
    ax.bar(2.60, nnrti_table["rf_mean"], yerr=nnrti_table["rf_std"], color="purple", alpha = 0.7, width=barWidth, edgecolor="black", label="rf")
    ax.bar(2.75, nnrti_table["CNN_mean"], yerr=nnrti_table["CNN_std"], color="cyan", alpha = 0.7, width=barWidth, edgecolor="black", label="CNN")
    ax.bar(2.90, nnrti_table["small_ensemble_mean"], yerr=nnrti_table["small_ensemble_std"], color="orange", alpha = 0.7, width=barWidth, edgecolor="black", label="small_ensemble")
    ax.bar(3.05, nnrti_table["large_ensemble_mean"], yerr=nnrti_table["large_ensemble_std"], color="pink", alpha = 0.7, width=barWidth, edgecolor="black", label="large_ensemble")

    ###NRTIs
    ax.bar(3.70, nrti_table["HIVDB_mean"], yerr=nrti_table["HIVDB_std"], color="b", alpha = 0.7,  width=barWidth, edgecolor="black", label="HIVDB")
    ax.bar(3.85, nrti_table["LSR_mean"], yerr=nrti_table["LSR_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR")
    ax.bar(4.00, nrti_table["LSR_I_mean"], yerr=nrti_table["LSR_I_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR_I", hatch="//")
    ax.bar(4.15, nrti_table["LARS_mean"], yerr=nrti_table["LARS_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS")
    ax.bar(4.30, nrti_table["LARS_I_mean"], yerr=nrti_table["LARS_I_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS_I", hatch="//")
    ax.bar(4.45, nrti_table["rf_mean"], yerr=nrti_table["rf_std"], color="purple", alpha = 0.7, width=barWidth, edgecolor="black", label="rf")
    ax.bar(4.60, nrti_table["CNN_mean"], yerr=nrti_table["CNN_std"], color="cyan", alpha = 0.7, width=barWidth, edgecolor="black", label="CNN")
    ax.bar(4.75, nrti_table["small_ensemble_mean"], yerr=nrti_table["small_ensemble_std"], color="orange", alpha = 0.7, width=barWidth, edgecolor="black", label="small_ensemble")
    ax.bar(4.90, nrti_table["large_ensemble_mean"], yerr=nrti_table["large_ensemble_std"], color="pink", alpha = 0.7, width=barWidth, edgecolor="black", label="large_ensemble")

    ###PIs
    ax.bar(5.55, pi_table["HIVDB_mean"], yerr=pi_table["HIVDB_std"], color="b", alpha = 0.7,  width=barWidth, edgecolor="black", label="HIVDB")
    ax.bar(5.70, pi_table["LSR_mean"], yerr=pi_table["LSR_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR")
    ax.bar(5.85, pi_table["LSR_I_mean"], yerr=pi_table["LSR_I_std"], color="r", alpha = 0.7, width=barWidth, edgecolor="black", label="LSR_I", hatch="//")
    ax.bar(6.00, pi_table["LARS_mean"], yerr=pi_table["LARS_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS")
    ax.bar(6.15, pi_table["LARS_I_mean"], yerr=pi_table["LARS_I_std"], color="y", alpha = 0.7, width=barWidth, edgecolor="black", label="LARS_I", hatch="//")
    ax.bar(6.30, pi_table["rf_mean"], yerr=pi_table["rf_std"], color="purple", alpha = 0.7, width=barWidth, edgecolor="black", label="rf")
    ax.bar(6.45, pi_table["CNN_mean"], yerr=pi_table["CNN_std"], color="cyan", alpha = 0.7, width=barWidth, edgecolor="black", label="CNN")
    ax.bar(6.60, pi_table["small_ensemble_mean"], yerr=pi_table["small_ensemble_std"], color="orange", alpha = 0.7, width=barWidth, edgecolor="black", label="small_ensemble")
    ax.bar(6.75, pi_table["large_ensemble_mean"], yerr=pi_table["large_ensemble_std"], color="pink", alpha = 0.7, width=barWidth, edgecolor="black", label="large_ensemble")

    ax.set_xticks([0.675, 2.45, 4.30, 6.15])
    ax.set_xticklabels(performance_metric_table["Drug_Class"])
    ax.set_ylabel(f"{metric.replace("_"," ")} (5 folds)")
    ax.set_xlabel("Drug Class")
    ax.set_title(f"{metric.replace("_"," ")} per Drug Class")

    ax.legend(loc="upper left", bbox_to_anchor=(1,1), ncol=1, title="Method", labels=methods)
    legend_elements = [Patch(facecolor="b", edgecolor="black", alpha = 0.7, label="HIVDB"),
                    Patch(facecolor="r", edgecolor="black", alpha = 0.7, label="LSR"),
                    Patch(facecolor="r", edgecolor="black", alpha = 0.7, label="LSR-I", hatch="//"),
                    Patch(facecolor="y", edgecolor="black", alpha = 0.7, label="LARS"),
                    Patch(facecolor="y", edgecolor="black", alpha = 0.7, label="LARS-I", hatch="//"),
                    Patch(facecolor="purple", edgecolor="black", alpha = 0.7, label="RandomForest"),
                    Patch(facecolor="cyan", edgecolor="black", alpha = 0.7, label="CNN"),
                    Patch(facecolor="orange", edgecolor="black", alpha = 0.7, label="HIVDB+LSR-I+RF"),
                    Patch(facecolor="pink", edgecolor="black", alpha = 0.7, label="HIVDB+LSR-I+LARS-I\n+RF+CNN")]

    ax.legend(handles = legend_elements, loc="upper left", bbox_to_anchor=(1,1), ncol=1, title="Method", fontsize="small")
    plt.tight_layout()
    plt.savefig(f"figures/method_performance_{metric.lower()}_barplot.png")
    print(f"{metric.replace('_',' ')} barplot saved as figures/method_performance_{metric.lower()}_barplot.png")



#PERFORMANCE HEATMAPS (Supplementary Figures 1-4)
for metric in ["Accuracy", "Balanced_accuracy", "F1_score", "MCC"]:
    heatmap_table = pd.DataFrame(index=test_drugs, columns=methods)
    heatmap_table = heatmap_table.fillna(float(0)).infer_objects(copy=False)
    for drug in test_drugs:
        heatmap_table.loc[drug, 'HIVDB'] = float(hivdb[hivdb['Drug'] == drug][metric].mean())
        heatmap_table.loc[drug, 'LSR'] = float(lsr[lsr['Drug'] == drug][metric].mean())
        heatmap_table.loc[drug, 'LSR-I'] = float(lsr_i[lsr_i['Drug'] == drug][metric].mean())
        heatmap_table.loc[drug, 'LARS'] = float(lars[lars['Drug'] == drug][metric].mean())
        heatmap_table.loc[drug, 'LARS-I'] = float(lars_i[lars_i['Drug'] == drug][metric].mean())
        heatmap_table.loc[drug, 'RandomForest'] = float(random_forest[random_forest['Drug'] == drug][metric].mean())
        if pd.isnull(cnn[cnn['Drug'] == drug][metric].mean()): #for INI drugs there is no CNN predictions
            heatmap_table.loc[drug, 'CNN'] = '-'
        else:
            heatmap_table.loc[drug, 'CNN'] = float(cnn[cnn['Drug'] == drug][metric].mean())

        heatmap_table.loc[drug, 'HIVDB+LSR-I+RF'] = float(small_ensemble[small_ensemble['Drug'] == drug][metric].mean())
        heatmap_table.loc[drug, 'HIVDB+LSR-I+LARS-I+RF+CNN'] = float(large_ensemble[large_ensemble['Drug'] == drug][metric].mean())

    #we add a Mean row of the tested drugs
    for j in range(heatmap_table.shape[1]):
        values = heatmap_table.loc[test_drugs, heatmap_table.columns[j]].values
        values = [x for x in values if x != '-']
        mean = sum(values)/len(values)
        heatmap_table.loc['Mean', heatmap_table.columns[j]] = mean

    #We add a Mean without INIs
    for j in range(heatmap_table.shape[1]):
        values = heatmap_table.loc[test_drugs_no_INIs, heatmap_table.columns[j]].values
        values = [x for x in values if x != '-']# and x != 0]
        mean = sum(values)/len(values)
        heatmap_table.loc['Mean_WO_INIs', heatmap_table.columns[j]] = mean

    #we add a row for the mean of drugs in each class
    heatmap_table.loc['INI'] = heatmap_table.loc[INIs].mean(numeric_only=True)
    heatmap_table.loc['PI'] = heatmap_table.loc[PIs].mean()
    heatmap_table.loc['NNRTI'] = heatmap_table.loc[NNRTIs].mean()
    heatmap_table.loc['NRTI'] = heatmap_table.loc[NRTIs].mean()
    
    heatmap_table = heatmap_table.reindex(['RAL', 'EVG', 'DTG', 'BIC', 'INI', 'EFV', 'NVP', 'ETR', 'RPV', 'NNRTI', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'NRTI', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV', 'PI', 'Mean_WO_INIs','Mean'])
    heatmap_table = heatmap_table.round(4)

    #We plot the heatmaps
    heatmap_table = heatmap_table.apply(pd.to_numeric, errors='coerce')
    
    null_coordinates = []   #we get the coordinates of the NaN values
    min_val, max_val = 1, 0
    for i in range(heatmap_table.shape[0]):
        for j in range(heatmap_table.shape[1]):
            if pd.isnull(heatmap_table.iloc[i, j]):
                null_coordinates.append((i, j))
            else:
                if heatmap_table.iloc[i, j] < min_val:
                    min_val = heatmap_table.iloc[i, j]
                if heatmap_table.iloc[i, j] > max_val:
                    max_val = heatmap_table.iloc[i, j]

    sns.set(font_scale=1.5)
    fig, ax = plt.subplots(figsize=(20, 20))
    sns.heatmap(heatmap_table, cmap='crest', ax=ax, vmin=min_val, vmax=max_val)
    ax.hlines([4, 5, 9, 10, 16, 17, 25, 26,27], *ax.get_xlim(), color='black', lw=2, edgecolor='black')
    ax.yaxis.set_tick_params(rotation=0)
    for i in [4, 9, 16, 25, 26,27]:
        ax.get_yticklabels()[i].set_weight('bold')
        ax.get_yticklabels()[i].set_size(15)

    for i, j in null_coordinates:
        ax.add_patch(plt.Rectangle((j, i), 1, 1, color = 'white',  lw=2))

    #we plot the values in the cells
    for i in range(heatmap_table.shape[0]):
        for j in range(heatmap_table.shape[1]):
            if pd.notnull(heatmap_table.iloc[i, j]) and not i in [4, 9, 16, 25, 26,27]:
                ax.text(j+0.5, i+0.5, round(heatmap_table.iloc[i, j], 3), ha='center', va='center', color='white', fontsize=18)
            elif i in [4, 9, 16, 25, 26,27]:
                ax.text(j+0.5, i+0.5, round(heatmap_table.iloc[i, j], 3), ha='center', va='center', color='white', fontsize=18, fontweight='bold')

    ax.set_title(f'{metric.replace("_", " ")} of the different methods\n', fontsize=20, fontweight='bold')
    ax.set_xlabel('Methods', fontsize=20, fontweight='bold')
    ax.xaxis.set_label_coords(0.5, -0.1)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=25, horizontalalignment='right')
    ax.set_ylabel('Drugs', fontsize=20, fontweight='bold')

    labels = ax.get_yticklabels()
    labels[26] = 'Mean (without INIs)'
    ax.set_yticklabels(labels)
    ax.set_xticklabels(["HIVDB", "LSR", "LSR-I", "LARS", "LARS-I", "RandomForest", "CNN", "HIVDB+LSR-I+RF", "HIVDB+LSR-I+LARS-I+RF+CNN"], rotation=25, horizontalalignment='right')

    plt.savefig(f"figures/method_performance_{metric.lower()}_heatmap.png")
    print(f"{metric.replace('_',' ')} heatmap saved as figures/method_performance_{metric.lower()}_heatmap.png")