###Drug dataset and fold division

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold

# Load the prefiltered datasets
ini_dataset = pd.read_csv("datasets/INI_dataset.tsv", sep='\t')
nnrti_dataset = pd.read_csv("datasets/NNRTI_dataset.tsv", sep='\t')
nrti_dataset = pd.read_csv("datasets/NRTI_dataset.tsv", sep='\t')
pi_dataset = pd.read_csv("datasets/PI_dataset.tsv", sep='\t')

##we separate the datasets by drug class
ini_dataset_RAL = ini_dataset[ini_dataset['RAL'].notna()].copy()
ini_dataset_EVG = ini_dataset[ini_dataset['EVG'].notna()].copy()
ini_dataset_DTG = ini_dataset[ini_dataset['DTG'].notna()].copy()
ini_dataset_BIC = ini_dataset[ini_dataset['BIC'].notna()].copy()

nnrti_dataset_EFV = nnrti_dataset[nnrti_dataset['EFV'].notna()].copy()
nnrti_dataset_NVP = nnrti_dataset[nnrti_dataset['NVP'].notna()].copy()
nnrti_dataset_ETR = nnrti_dataset[nnrti_dataset['ETR'].notna()].copy()
nnrti_dataset_RPV = nnrti_dataset[nnrti_dataset['RPV'].notna()].copy()

nrti_dataset_3TC = nrti_dataset[nrti_dataset['3TC'].notna()].copy()
nrti_dataset_ABC = nrti_dataset[nrti_dataset['ABC'].notna()].copy()
nrti_dataset_AZT = nrti_dataset[nrti_dataset['AZT'].notna()].copy()
nrti_dataset_D4T = nrti_dataset[nrti_dataset['D4T'].notna()].copy()
nrti_dataset_DDI = nrti_dataset[nrti_dataset['DDI'].notna()].copy()
nrti_dataset_TDF = nrti_dataset[nrti_dataset['TDF'].notna()].copy()

pi_dataset_FPV = pi_dataset[pi_dataset['FPV'].notna()].copy()
pi_dataset_ATV = pi_dataset[pi_dataset['ATV'].notna()].copy()
pi_dataset_IDV = pi_dataset[pi_dataset['IDV'].notna()].copy()
pi_dataset_LPV = pi_dataset[pi_dataset['LPV'].notna()].copy()
pi_dataset_NFV = pi_dataset[pi_dataset['NFV'].notna()].copy()
pi_dataset_SQV = pi_dataset[pi_dataset['SQV'].notna()].copy()
pi_dataset_TPV = pi_dataset[pi_dataset['TPV'].notna()].copy()
pi_dataset_DRV = pi_dataset[pi_dataset['DRV'].notna()].copy()

dataset_list = [ini_dataset_RAL, ini_dataset_EVG, ini_dataset_DTG, ini_dataset_BIC,
                nnrti_dataset_EFV, nnrti_dataset_NVP, nnrti_dataset_ETR, nnrti_dataset_RPV,
                nrti_dataset_3TC, nrti_dataset_ABC, nrti_dataset_AZT, nrti_dataset_D4T, nrti_dataset_DDI, nrti_dataset_TDF,
                pi_dataset_FPV, pi_dataset_ATV, pi_dataset_IDV, pi_dataset_LPV, pi_dataset_NFV, pi_dataset_SQV, pi_dataset_TPV, pi_dataset_DRV]

###We consider as different patients the PtID 11630
###We will change the PtID 11630 to a different random number which is not present in the PtID column already
for dataset in dataset_list:
    index_list = dataset[dataset['PtID'] == 11630].index
    for index in index_list:
        for i in range(dataset.shape[0]):
            if i not in dataset['PtID'].values:
                dataset.at[index, 'PtID'] = i
                break

###WE CREATE 5 FOLDS OUT OF ALL THE DATA, STRATIFIED AND WITHOUT PATIENT DATA LEAKAGE
#we define the StratifiedKFoldGroup object
sgkf = StratifiedGroupKFold(n_splits=5, random_state=42, shuffle = True)

##INI
##RAL
labels = np.array(ini_dataset_RAL['RAL_res'])
seq_id = np.array(ini_dataset_RAL['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = ini_dataset_RAL, y = ini_dataset_RAL["RAL_res"], groups = ini_dataset_RAL["PtID"])):
    seq_index = ini_dataset_RAL.index[test_index]
    ini_dataset_RAL.loc[seq_index, 'fold'] = i

ini_dataset_RAL.to_csv("datasets/INI/INI_RAL_5folds.tsv", sep = "\t", index = False)

#we do the same for the other datasets
##EVG
labels = np.array(ini_dataset_EVG['EVG_res'])
seq_id = np.array(ini_dataset_EVG['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = ini_dataset_EVG, y = ini_dataset_EVG["EVG_res"], groups = ini_dataset_EVG["PtID"])):
    seq_index = ini_dataset_EVG.index[test_index]
    ini_dataset_EVG.loc[seq_index, 'fold'] = i

ini_dataset_EVG.to_csv("datasets/INI/INI_EVG_5folds.tsv", sep = "\t", index = False)
        
##DTG
labels = np.array(ini_dataset_DTG['DTG_res'])
seq_id = np.array(ini_dataset_DTG['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = ini_dataset_DTG, y = ini_dataset_DTG["DTG_res"], groups = ini_dataset_DTG["PtID"])):
    seq_index = ini_dataset_DTG.index[test_index]
    ini_dataset_DTG.loc[seq_index, 'fold'] = i

ini_dataset_DTG.to_csv("datasets/INI/INI_DTG_5folds.tsv", sep = "\t", index = False)

##BIC
labels = np.array(ini_dataset_BIC['BIC_res'])
seq_id = np.array(ini_dataset_BIC['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = ini_dataset_BIC, y = ini_dataset_BIC["BIC_res"], groups = ini_dataset_BIC["PtID"])):
    seq_index = ini_dataset_BIC.index[test_index]
    ini_dataset_BIC.loc[seq_index, 'fold'] = i

ini_dataset_BIC.to_csv("datasets/INI/INI_BIC_5folds.tsv", sep = "\t", index = False)

##NNRTI
##EFV
labels = np.array(nnrti_dataset_EFV['EFV_res'])
seq_id = np.array(nnrti_dataset_EFV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nnrti_dataset_EFV, y = nnrti_dataset_EFV["EFV_res"], groups = nnrti_dataset_EFV["PtID"])):
    seq_index = nnrti_dataset_EFV.index[test_index]
    nnrti_dataset_EFV.loc[seq_index, 'fold'] = i

nnrti_dataset_EFV.to_csv("datasets/NNRTI/NNRTI_EFV_5folds.tsv", sep = "\t", index = False)


##NVP
labels = np.array(nnrti_dataset_NVP['NVP_res'])
seq_id = np.array(nnrti_dataset_NVP['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nnrti_dataset_NVP, y = nnrti_dataset_NVP["NVP_res"], groups = nnrti_dataset_NVP["PtID"])):
    seq_index = nnrti_dataset_NVP.index[test_index]
    nnrti_dataset_NVP.loc[seq_index, 'fold'] = i

nnrti_dataset_NVP.to_csv("datasets/NNRTI/NNRTI_NVP_5folds.tsv", sep = "\t", index = False)

##ETR
labels = np.array(nnrti_dataset_ETR['ETR_res'])
seq_id = np.array(nnrti_dataset_ETR['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nnrti_dataset_ETR, y = nnrti_dataset_ETR["ETR_res"], groups = nnrti_dataset_ETR["PtID"])):
    seq_index = nnrti_dataset_ETR.index[test_index]
    nnrti_dataset_ETR.loc[seq_index, 'fold'] = i

nnrti_dataset_ETR.to_csv("datasets/NNRTI/NNRTI_ETR_5folds.tsv", sep = "\t", index = False)

##RPV
labels = np.array(nnrti_dataset_RPV['RPV_res'])
seq_id = np.array(nnrti_dataset_RPV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nnrti_dataset_RPV, y = nnrti_dataset_RPV["RPV_res"], groups = nnrti_dataset_RPV["PtID"])):
    seq_index = nnrti_dataset_RPV.index[test_index]
    nnrti_dataset_RPV.loc[seq_index, 'fold'] = i

nnrti_dataset_RPV.to_csv("datasets/NNRTI/NNRTI_RPV_5folds.tsv", sep = "\t", index = False)

###NRTI
##3TC
labels = np.array(nrti_dataset_3TC['3TC_res'])
seq_id = np.array(nrti_dataset_3TC['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nrti_dataset_3TC, y = nrti_dataset_3TC["3TC_res"], groups = nrti_dataset_3TC["PtID"])):
    seq_index = nrti_dataset_3TC.index[test_index]
    nrti_dataset_3TC.loc[seq_index, 'fold'] = i

nrti_dataset_3TC.to_csv("datasets/NRTI/NRTI_3TC_5folds.tsv", sep = "\t", index = False)

##ABC
labels = np.array(nrti_dataset_ABC['ABC_res'])
seq_id = np.array(nrti_dataset_ABC['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nrti_dataset_ABC, y = nrti_dataset_ABC["ABC_res"], groups = nrti_dataset_ABC["PtID"])):
    seq_index = nrti_dataset_ABC.index[test_index]
    nrti_dataset_ABC.loc[seq_index, 'fold'] = i

nrti_dataset_ABC.to_csv("datasets/NRTI/NRTI_ABC_5folds.tsv", sep = "\t", index = False)

##AZT
labels = np.array(nrti_dataset_AZT['AZT_res'])
seq_id = np.array(nrti_dataset_AZT['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nrti_dataset_AZT, y = nrti_dataset_AZT["AZT_res"], groups = nrti_dataset_AZT["PtID"])):
    seq_index = nrti_dataset_AZT.index[test_index]
    nrti_dataset_AZT.loc[seq_index, 'fold'] = i

nrti_dataset_AZT.to_csv("datasets/NRTI/NRTI_AZT_5folds.tsv", sep = "\t", index = False)

##D4T
labels = np.array(nrti_dataset_D4T['D4T_res'])
seq_id = np.array(nrti_dataset_D4T['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nrti_dataset_D4T, y = nrti_dataset_D4T["D4T_res"], groups = nrti_dataset_D4T["PtID"])):
    seq_index = nrti_dataset_D4T.index[test_index]
    nrti_dataset_D4T.loc[seq_index, 'fold'] = i

nrti_dataset_D4T.to_csv("datasets/NRTI/NRTI_D4T_5folds.tsv", sep = "\t", index = False)

##DDI
labels = np.array(nrti_dataset_DDI['DDI_res'])
seq_id = np.array(nrti_dataset_DDI['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nrti_dataset_DDI, y = nrti_dataset_DDI["DDI_res"], groups = nrti_dataset_DDI["PtID"])):
    seq_index = nrti_dataset_DDI.index[test_index]
    nrti_dataset_DDI.loc[seq_index, 'fold'] = i

nrti_dataset_DDI.to_csv("datasets/NRTI/NRTI_DDI_5folds.tsv", sep = "\t", index = False)

##TDF
labels = np.array(nrti_dataset_TDF['TDF_res'])
seq_id = np.array(nrti_dataset_TDF['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = nrti_dataset_TDF, y = nrti_dataset_TDF["TDF_res"], groups = nrti_dataset_TDF["PtID"])):
    seq_index = nrti_dataset_TDF.index[test_index]
    nrti_dataset_TDF.loc[seq_index, 'fold'] = i

nrti_dataset_TDF.to_csv("datasets/NRTI/NRTI_TDF_5folds.tsv", sep = "\t", index = False)

###PI
##FPV
labels = np.array(pi_dataset_FPV['FPV_res'])
seq_id = np.array(pi_dataset_FPV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_FPV, y = pi_dataset_FPV["FPV_res"], groups = pi_dataset_FPV["PtID"])):
    seq_index = pi_dataset_FPV.index[test_index]
    pi_dataset_FPV.loc[seq_index, 'fold'] = i

pi_dataset_FPV.to_csv("datasets/PI/PI_FPV_5folds.tsv", sep = "\t", index = False)

##ATV
labels = np.array(pi_dataset_ATV['ATV_res'])
seq_id = np.array(pi_dataset_ATV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_ATV, y = pi_dataset_ATV["ATV_res"], groups = pi_dataset_ATV["PtID"])):
    seq_index = pi_dataset_ATV.index[test_index]
    pi_dataset_ATV.loc[seq_index, 'fold'] = i

pi_dataset_ATV.to_csv("datasets/PI/PI_ATV_5folds.tsv", sep = "\t", index = False)

##IDV
labels = np.array(pi_dataset_IDV['IDV_res'])
seq_id = np.array(pi_dataset_IDV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_IDV, y = pi_dataset_IDV["IDV_res"], groups = pi_dataset_IDV["PtID"])):
    seq_index = pi_dataset_IDV.index[test_index]
    pi_dataset_IDV.loc[seq_index, 'fold'] = i

pi_dataset_IDV.to_csv("datasets/PI/PI_IDV_5folds.tsv", sep = "\t", index = False)

##LPV
labels = np.array(pi_dataset_LPV['LPV_res'])
seq_id = np.array(pi_dataset_LPV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_LPV, y = pi_dataset_LPV["LPV_res"], groups = pi_dataset_LPV["PtID"])):
    seq_index = pi_dataset_LPV.index[test_index]
    pi_dataset_LPV.loc[seq_index, 'fold'] = i

pi_dataset_LPV.to_csv("datasets/PI/PI_LPV_5folds.tsv", sep = "\t", index = False)

##NFV
labels = np.array(pi_dataset_NFV['NFV_res'])
seq_id = np.array(pi_dataset_NFV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_NFV, y = pi_dataset_NFV["NFV_res"], groups = pi_dataset_NFV["PtID"])):
    seq_index = pi_dataset_NFV.index[test_index]
    pi_dataset_NFV.loc[seq_index, 'fold'] = i

pi_dataset_NFV.to_csv("datasets/PI/PI_NFV_5folds.tsv", sep = "\t", index = False)

##SQV
labels = np.array(pi_dataset_SQV['SQV_res'])
seq_id = np.array(pi_dataset_SQV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_SQV, y = pi_dataset_SQV["SQV_res"], groups = pi_dataset_SQV["PtID"])):
    seq_index = pi_dataset_SQV.index[test_index]
    pi_dataset_SQV.loc[seq_index, 'fold'] = i

pi_dataset_SQV.to_csv("datasets/PI/PI_SQV_5folds.tsv", sep = "\t", index = False)

##TPV
labels = np.array(pi_dataset_TPV['TPV_res'])
seq_id = np.array(pi_dataset_TPV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_TPV, y = pi_dataset_TPV["TPV_res"], groups = pi_dataset_TPV["PtID"])):
    seq_index = pi_dataset_TPV.index[test_index]
    pi_dataset_TPV.loc[seq_index, 'fold'] = i

pi_dataset_TPV.to_csv("datasets/PI/PI_TPV_5folds.tsv", sep = "\t", index = False)

##DRV
labels = np.array(pi_dataset_DRV['DRV_res'])
seq_id = np.array(pi_dataset_DRV['SeqID'])
for i, (train_index, test_index) in enumerate(sgkf.split(X = pi_dataset_DRV, y = pi_dataset_DRV["DRV_res"], groups = pi_dataset_DRV["PtID"])):
    seq_index = pi_dataset_DRV.index[test_index]
    pi_dataset_DRV.loc[seq_index, 'fold'] = i

pi_dataset_DRV.to_csv("datasets/PI/PI_DRV_5folds.tsv", sep = "\t", index = False)
