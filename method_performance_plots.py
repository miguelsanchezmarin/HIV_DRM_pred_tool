###We create barplots with error bars to compare the different methods performance
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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

methods = ["HIVDB", "LSR", "LSR_I", "LARS", "LARS_comb", "G2P_RF", "G2P_Prob", "rf",  "CNN", "HIVDB+LSR_C+RF_P", "HIVDB+LSR_C+LARS_C+RF_P+CNN"]

##ACCURACY
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

###Now we do a bar plot where we will plot the mean as bars and +-std in the errror bars.

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
