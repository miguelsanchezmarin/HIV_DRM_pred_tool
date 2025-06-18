###Analysis of the prefiltered high quality Stanford HIV dataset

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


ini_dataset = pd.read_csv("datasets/INI_dataset.tsv", sep='\t')
nnrti_dataset = pd.read_csv("datasets/NNRTI_dataset.tsv", sep='\t')
nrti_dataset = pd.read_csv("datasets/NRTI_dataset.tsv", sep='\t')
pi_dataset = pd.read_csv("datasets/PI_dataset.tsv", sep='\t')

used_drugs = ['ABC', 'FTC', '3TC', 'TDF', 'AZT', 'DOR', 'EFV', 'ETR', 'NVP', 'RPV', 'ATV', 'DRV', 'FPV', 'RTV', 'TPV', 'LPV', 'CAB', 'DTG', 'RAL']


###Barplot showing the different drugs dataset size and distribution (Figure 1)

fig, ax = plt.subplots()
ax.bar(['RAL', 'EVG', 'DTG', 'BIC', 'CAB'], [ini_dataset[ini_dataset['RAL_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['EVG_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['DTG_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['BIC_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['CAB_res'] == 'Resistant'].shape[0]], color='white', label='INIs (R)', alpha=1, edgecolor='blue', hatch='//')
ax.bar(['RAL', 'EVG', 'DTG', 'BIC', 'CAB'], [ini_dataset[ini_dataset['RAL_res'] == 'Susceptible'].shape[0], ini_dataset[ini_dataset['EVG_res'] == 'Susceptible'].shape[0], ini_dataset[ini_dataset['DTG_res'] == 'Susceptible'].shape[0], ini_dataset[ini_dataset['BIC_res'] == 'Susceptible'].shape[0], ini_dataset[ini_dataset['CAB_res'] == 'Susceptible'].shape[0]],
       bottom=[ini_dataset[ini_dataset['RAL_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['EVG_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['DTG_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['BIC_res'] == 'Resistant'].shape[0], ini_dataset[ini_dataset['CAB_res'] == 'Resistant'].shape[0]], color='blue', label='INIs (S)', alpha=1, edgecolor='black')
ax.bar(['EFV', 'NVP', 'ETR', 'RPV'], [nnrti_dataset[nnrti_dataset['EFV_res'] == 'Resistant'].shape[0], nnrti_dataset[nnrti_dataset['NVP_res'] == 'Resistant'].shape[0], nnrti_dataset[nnrti_dataset['ETR_res'] == 'Resistant'].shape[0], nnrti_dataset[nnrti_dataset['RPV_res'] == 'Resistant'].shape[0]], color='white', label='NNRTIs (R)', alpha=1, edgecolor='red', hatch='//')
ax.bar(['EFV', 'NVP', 'ETR', 'RPV'], [nnrti_dataset[nnrti_dataset['EFV_res'] == 'Susceptible'].shape[0], nnrti_dataset[nnrti_dataset['NVP_res'] == 'Susceptible'].shape[0], nnrti_dataset[nnrti_dataset['ETR_res'] == 'Susceptible'].shape[0], nnrti_dataset[nnrti_dataset['RPV_res'] == 'Susceptible'].shape[0]],
       bottom=[nnrti_dataset[nnrti_dataset['EFV_res'] == 'Resistant'].shape[0], nnrti_dataset[nnrti_dataset['NVP_res'] == 'Resistant'].shape[0], nnrti_dataset[nnrti_dataset['ETR_res'] == 'Resistant'].shape[0], nnrti_dataset[nnrti_dataset['RPV_res'] == 'Resistant'].shape[0]], color='red', label='NNRTIs (S)', alpha=1, edgecolor='black')
ax.bar(['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF'], [nrti_dataset[nrti_dataset['3TC_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['ABC_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['AZT_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['D4T_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['DDI_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['TDF_res'] == 'Resistant'].shape[0]], color='white', label='NRTIs (R)', alpha=1, edgecolor='green', hatch='//')
ax.bar(['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF'], [nrti_dataset[nrti_dataset['3TC_res'] == 'Susceptible'].shape[0], nrti_dataset[nrti_dataset['ABC_res'] == 'Susceptible'].shape[0], nrti_dataset[nrti_dataset['AZT_res'] == 'Susceptible'].shape[0], nrti_dataset[nrti_dataset['D4T_res'] == 'Susceptible'].shape[0], nrti_dataset[nrti_dataset['DDI_res'] == 'Susceptible'].shape[0], nrti_dataset[nrti_dataset['TDF_res'] == 'Susceptible'].shape[0]],
       bottom=[nrti_dataset[nrti_dataset['3TC_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['ABC_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['AZT_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['D4T_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['DDI_res'] == 'Resistant'].shape[0], nrti_dataset[nrti_dataset['TDF_res'] == 'Resistant'].shape[0]], color='green', label='NRTIs (S)', alpha=1, edgecolor='black')
ax.bar(['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV'], [pi_dataset[pi_dataset['FPV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['ATV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['IDV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['LPV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['NFV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['SQV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['TPV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['DRV_res'] == 'Resistant'].shape[0]], color='white', label='PIs (R)', alpha=1, edgecolor='purple', hatch='//')
ax.bar(['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV'], [pi_dataset[pi_dataset['FPV_res'] == 'Susceptible'].shape[0], pi_dataset[pi_dataset['ATV_res'] == 'Susceptible'].shape[0], pi_dataset[pi_dataset['IDV_res'] == 'Susceptible'].shape[0], pi_dataset[pi_dataset['LPV_res'] == 'Susceptible'].shape[0], pi_dataset[pi_dataset['NFV_res'] == 'Susceptible'].shape[0], pi_dataset[pi_dataset['SQV_res'] == 'Susceptible'].shape[0], pi_dataset[pi_dataset['TPV_res'] == 'Susceptible'].shape[0], pi_dataset[pi_dataset['DRV_res'] == 'Susceptible'].shape[0]],
       bottom=[pi_dataset[pi_dataset['FPV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['ATV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['IDV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['LPV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['NFV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['SQV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['TPV_res'] == 'Resistant'].shape[0], pi_dataset[pi_dataset['DRV_res'] == 'Resistant'].shape[0]], color='purple', label='PIs (S)', alpha=1, edgecolor='black')
                                                                  
ax.set_ylabel('Number of sequences')
ax.set_xticks(['RAL', 'EVG', 'DTG', 'BIC', 'CAB', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV'])
ax.set_xticklabels(['RAL', 'EVG', 'DTG', 'BIC', 'CAB', 'EFV', 'NVP', 'ETR', 'RPV', '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV'], rotation=45)
ax.set_ylim(0, 2500)
ax.set_title(f'Number of sequences per drug\n(StanfordHIVDB prefiltered high-quality dataset)')

#we plot the number of susceptible sequences on top of the bars and resistant sequences on top of the resistant bars
ax.text('RAL', ini_dataset[ini_dataset['RAL_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['RAL_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('RAL', ini_dataset[ini_dataset['RAL_res'] == 'Susceptible'].shape[0] + ini_dataset[ini_dataset['RAL_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['RAL_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('EVG', ini_dataset[ini_dataset['EVG_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['EVG_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('EVG', ini_dataset[ini_dataset['EVG_res'] == 'Susceptible'].shape[0] + ini_dataset[ini_dataset['EVG_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['EVG_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('DTG', ini_dataset[ini_dataset['DTG_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['DTG_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('DTG', ini_dataset[ini_dataset['DTG_res'] == 'Susceptible'].shape[0] + ini_dataset[ini_dataset['DTG_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['DTG_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('BIC', ini_dataset[ini_dataset['BIC_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['BIC_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('BIC', ini_dataset[ini_dataset['BIC_res'] == 'Susceptible'].shape[0] + ini_dataset[ini_dataset['BIC_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['BIC_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('CAB', ini_dataset[ini_dataset['CAB_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['CAB_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('CAB', ini_dataset[ini_dataset['CAB_res'] == 'Susceptible'].shape[0] + ini_dataset[ini_dataset['CAB_res'] == 'Resistant'].shape[0] + 10, str(ini_dataset[ini_dataset['CAB_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)

ax.text('EFV', nnrti_dataset[nnrti_dataset['EFV_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['EFV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('EFV', nnrti_dataset[nnrti_dataset['EFV_res'] == 'Susceptible'].shape[0] + nnrti_dataset[nnrti_dataset['EFV_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['EFV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('NVP', nnrti_dataset[nnrti_dataset['NVP_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['NVP_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('NVP', nnrti_dataset[nnrti_dataset['NVP_res'] == 'Susceptible'].shape[0] + nnrti_dataset[nnrti_dataset['NVP_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['NVP_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('ETR', nnrti_dataset[nnrti_dataset['ETR_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['ETR_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('ETR', nnrti_dataset[nnrti_dataset['ETR_res'] == 'Susceptible'].shape[0] + nnrti_dataset[nnrti_dataset['ETR_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['ETR_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('RPV', nnrti_dataset[nnrti_dataset['RPV_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['RPV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('RPV', nnrti_dataset[nnrti_dataset['RPV_res'] == 'Susceptible'].shape[0] + nnrti_dataset[nnrti_dataset['RPV_res'] == 'Resistant'].shape[0] + 10, str(nnrti_dataset[nnrti_dataset['RPV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)

ax.text('3TC', nrti_dataset[nrti_dataset['3TC_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['3TC_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('3TC', nrti_dataset[nrti_dataset['3TC_res'] == 'Susceptible'].shape[0] + nrti_dataset[nrti_dataset['3TC_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['3TC_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('ABC', nrti_dataset[nrti_dataset['ABC_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['ABC_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('ABC', nrti_dataset[nrti_dataset['ABC_res'] == 'Susceptible'].shape[0] + nrti_dataset[nrti_dataset['ABC_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['ABC_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('AZT', nrti_dataset[nrti_dataset['AZT_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['AZT_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('AZT', nrti_dataset[nrti_dataset['AZT_res'] == 'Susceptible'].shape[0] + nrti_dataset[nrti_dataset['AZT_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['AZT_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('D4T', nrti_dataset[nrti_dataset['D4T_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['D4T_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('D4T', nrti_dataset[nrti_dataset['D4T_res'] == 'Susceptible'].shape[0] + nrti_dataset[nrti_dataset['D4T_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['D4T_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('DDI', nrti_dataset[nrti_dataset['DDI_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['DDI_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('DDI', nrti_dataset[nrti_dataset['DDI_res'] == 'Susceptible'].shape[0] + nrti_dataset[nrti_dataset['DDI_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['DDI_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('TDF', nrti_dataset[nrti_dataset['TDF_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['TDF_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('TDF', nrti_dataset[nrti_dataset['TDF_res'] == 'Susceptible'].shape[0] + nrti_dataset[nrti_dataset['TDF_res'] == 'Resistant'].shape[0] + 10, str(nrti_dataset[nrti_dataset['TDF_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)

ax.text('FPV', pi_dataset[pi_dataset['FPV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['FPV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('FPV', pi_dataset[pi_dataset['FPV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['FPV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['FPV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('ATV', pi_dataset[pi_dataset['ATV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['ATV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('ATV', pi_dataset[pi_dataset['ATV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['ATV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['ATV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('IDV', pi_dataset[pi_dataset['IDV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['IDV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('IDV', pi_dataset[pi_dataset['IDV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['IDV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['IDV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('LPV', pi_dataset[pi_dataset['LPV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['LPV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('LPV', pi_dataset[pi_dataset['LPV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['LPV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['LPV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('NFV', pi_dataset[pi_dataset['NFV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['NFV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('NFV', pi_dataset[pi_dataset['NFV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['NFV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['NFV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('SQV', pi_dataset[pi_dataset['SQV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['SQV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('SQV', pi_dataset[pi_dataset['SQV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['SQV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['SQV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('TPV', pi_dataset[pi_dataset['TPV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['TPV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('TPV', pi_dataset[pi_dataset['TPV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['TPV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['TPV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('DRV', pi_dataset[pi_dataset['DRV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['DRV_res'] == 'Resistant'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)
ax.text('DRV', pi_dataset[pi_dataset['DRV_res'] == 'Susceptible'].shape[0] + pi_dataset[pi_dataset['DRV_res'] == 'Resistant'].shape[0] + 10, str(pi_dataset[pi_dataset['DRV_res'] == 'Susceptible'].shape[0]), ha='center', va='bottom', rotation=45, fontsize=7.5)

ax.legend(ncol=4, loc='upper right')

plt.savefig('figures/dataset_distribution_plot.png', dpi = 300)


###Dataset folds sizes and resistant/susceptible distributions (Supplementary Table 1)
fold_distribution = []
for drug_class in ["INI", "NNRTI", "NRTI", "PI"]:
       for dataset in os.listdir(f"datasets/{drug_class}"):
              drug = dataset.split("_")[1]
              dataset_table = pd.read_csv(f"datasets/{drug_class}/{dataset}", sep = "\t")
              res_table = dataset_table[dataset_table[f"{drug}_res"]=="Resistant"]
              sus_table = dataset_table[dataset_table[f"{drug}_res"]=="Susceptible"]
              fold_distribution.append([drug, "Total", dataset_table.shape[0], res_table.shape[0], sus_table.shape[0]])
              for fold in range(5):
                     dataset_table_fold = dataset_table[dataset_table["fold"]==fold]
                     res_table_fold = res_table[res_table["fold"]==fold]
                     sus_table_fold = sus_table[sus_table["fold"]==fold]
                     fold_distribution.append([drug, fold+1, dataset_table_fold.shape[0], res_table_fold.shape[0], sus_table_fold.shape[0]])

fold_distribution_df = pd.DataFrame(fold_distribution, columns = ["Drug", "Fold", "Total size", "Resistant size", "Susceptible size"])
fold_distribution_df.to_csv("datasets/dataset_fold_distribution.tsv", sep = "\t", index = False)
              
