import pandas as pd

# import matplotlib.pyplot as plt

##Drugs per drug class used in this study
INIs=['RAL', 'EVG', 'DTG', 'BIC']
NNRTIs=['EFV', 'NVP', 'ETR', 'RPV'] 
NRTIs=['3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF']
PIs=['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']


def get_coverage_dict(mut_dic: dict, robustness_data: dict, nreads:int = 1):
    '''Reads the coverage data and returns a dictionary with the information of % of DRM positions covered and estimated balanced accuracy
    INPUT:
    mut_dic: dictionary with the mutations and their frequencies and their coverage per position.
    robustness_data: dictionary with the robustness analysis data.
    nreads: number of reads to consider a position covered. Default is 1.
    
    OUTPUT:
    coverage_dict: dictionary with the % of covered DRM positions and the estimated balanced accuracy per drug class.
    '''
    coverage_dict = dict()
    coverage_dict['read_cutoff'] = nreads

    for dataset in ["INI", "NNRTI", "NRTI", "PI"]:
        
        coverage_dict[dataset] = dict()

        if dataset == "INI":
            drug_class = INIs
            drm_positions = [51, 66, 74, 75, 92, 95, 97, 118, 121, 122, 138, 140, 143, 145, 146, 147, 148, 151, 153, 155, 157, 163, 230, 232, 263]
            prot = "IN"
        elif dataset == "NNRTI":
            drug_class = NNRTIs
            drm_positions = [90, 98, 100, 101, 103, 106, 108, 138, 179, 181, 188, 190, 221, 225, 227, 230, 234, 236, 238, 318, 348]
            prot = "RT"
        elif dataset == "NRTI":
            drug_class = NRTIs
            drm_positions = [41, 62, 65, 67, 68, 69, 70, 74, 75, 77, 115, 116, 151, 184, 210, 215, 219]
            prot = "RT"
        elif dataset == "PI":
            drug_class = PIs
            drm_positions = [10, 20, 24, 32, 33, 46, 47, 48, 50, 53, 54, 73, 74, 76, 82, 83, 84, 88, 89, 90]
            prot = "PR"

        mut_dic_prot = mut_dic[prot]
        robustness_data_prot = robustness_data[dataset]
        covered_positions = []
        for pos in mut_dic_prot.keys():
            if pos in drm_positions:
                if mut_dic_prot[pos] >= nreads: #we check if the coverage is above the threshold
                    covered_positions.append(pos)
        
        drm_coverage = len(covered_positions)/len(drm_positions)*100 #we get the coverage as a % of the total number of DRM positions
        coverage_step = round(1-((int(drm_coverage/10)+1)/10),5)
        if drm_coverage == 100:
            coverage_step = 0
        if drm_coverage == 0:
            coverage_step = 1

        accuracies, max_accuracies = [], []
        for drug in drug_class:
            accuracies.append(robustness_data_prot[coverage_step][drug])
            max_accuracies.append(robustness_data_prot[0][drug])
        
        estimated_accuracy = sum(accuracies)/len(accuracies) #we get the average balanced accuracy for the drug class

        #we save the results in a dictionary
        coverage_dict[dataset]['coverage'] = drm_coverage
        coverage_dict[dataset]['balanced_accuracy'] = estimated_accuracy
        coverage_dict[dataset]['max_balanced_accuracy'] = sum(max_accuracies)/len(max_accuracies)

    return coverage_dict

def get_coverage_df(coverage_disclaimer_dict_list: list, sample_ids:list = None):
    '''Reads the coverage data list of dictionaries and returns a dataframe with the information of % of DRM positions covered and estimated balanced accuracy.
    INPUT:
    coverage_disclaimer_dict_list: list of dictionaries with the coverage data per sample.
    sample_ids: list of sample ids.
    
    OUTPUT:
    coverage_df: dataframe with the % of covered DRM positions and the estimated balanced accuracy per drug class.
    '''
    coverage_df = []
    if sample_ids is None:
        sample_coverage = coverage_disclaimer_dict_list
        read_cutoff = sample_coverage['read_cutoff']
        for dataset in sample_coverage.keys():
            if dataset == "read_cutoff":
                    continue
            coverage = sample_coverage[dataset]['coverage']
            balanced_accuracy = sample_coverage[dataset]['balanced_accuracy']
            max_balanced_accuracy = sample_coverage[dataset]['max_balanced_accuracy']
            coverage_df.append([dataset, coverage, balanced_accuracy, max_balanced_accuracy, read_cutoff])
        coverage_df = pd.DataFrame(coverage_df, columns=['Drug_Class', 'DRM_Coverage', 'Estimated_Balanced_Accuracy', 'Max_Balanced_Accuracy', 'Read_Cutoff'])

    else:
        for sample_id, sample_coverage in zip(sample_ids, coverage_disclaimer_dict_list):
            read_cutoff = sample_coverage['read_cutoff']
            for dataset in sample_coverage.keys():
                if dataset == "read_cutoff":
                    continue
                coverage = sample_coverage[dataset]['coverage']
                balanced_accuracy = sample_coverage[dataset]['balanced_accuracy']
                max_balanced_accuracy = sample_coverage[dataset]['max_balanced_accuracy']
                coverage_df.append([sample_id, dataset, coverage, balanced_accuracy, max_balanced_accuracy, read_cutoff])

        coverage_df = pd.DataFrame(coverage_df, columns=['Sample_ID', 'Drug_Class', 'DRM_Coverage', 'Estimated_Balanced_Accuracy', 'Max_Balanced_Accuracy', 'Read_Cutoff'])

    return coverage_df
        

