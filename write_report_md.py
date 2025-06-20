import pandas as pd
import os
import matplotlib.pyplot as plt


from load_data import load_robustness_data

def write_HIVDB_table(HIVDB_table, handle, unknown: bool = False, comments: bool = False):
    '''Writes the ensemble results to a markdown file.
    INPUT:
    HIVDB_table: Dataframe with the single mutations, frequency and HIVDB comments. Output from HIVDB_table function.
    handle: file handle to write the results.
    unknown: if True, the unknown annotations will be written.
    comments: if True, the comments will be written.
    '''
    cutoff = HIVDB_table.iloc[0]['Cutoff']

    handle.write("Single mutation annotations obtained from HIVDB program.")
    handle.write(" Mutations below " + str(cutoff) + " frequency were not included in the analysis.\n")
    IN_table, RT_table, PR_table = HIVDB_table[HIVDB_table["Protein"]=="IN"], HIVDB_table[HIVDB_table["Protein"]=="RT"], HIVDB_table[HIVDB_table["Protein"]=="PR"]
    if not unknown:
        IN_table, RT_table, PR_table = IN_table[IN_table["Annotation"] != "Unknown"], RT_table[RT_table["Annotation"] != "Unknown"], PR_table[PR_table["Annotation"] != "Unknown"]

    if IN_table.shape[0] > 0:
        handle.write("**INTEGRASE** mutations:\n\n")
        handle.write("|{: ^14}|{: ^11}|{: ^12}|\n".format('Mutation','Frequency', 'Annotation'))
        handle.write('|:{:-^12}:|:{:-^9}:|:{:-^10}:|'.format('', '', ''))
        handle.write('\n')
        for n in range(IN_table.shape[0]):    
            handle.write('|{: ^14}|{: ^11}|{: ^12}|'.format(IN_table["Mutation"].iloc[n], IN_table["Freq"].iloc[n], IN_table["Annotation"].iloc[n]))
            handle.write("\n")
        
        if comments:
            handle.write("\nComments:\n\n")
            for n in range(IN_table.shape[0]):
                if IN_table["Comment"].iloc[n] != "Unknown":
                    handle.write(f"- **{IN_table['Mutation'].iloc[n]}**: {IN_table['Comment'].iloc[n]}\n") 
    else:
        handle.write("\nNo **INTEGRASE** mutations found.\n")
    
    
    if RT_table.shape[0] > 0:
        handle.write("\n**REVERSE TRANSCRIPTASE** mutations:\n\n")
        handle.write("|{: ^14}|{: ^11}|{: ^12}|\n".format('Mutation','Frequency', 'Annotation'))
        handle.write('|:{:-^12}:|:{:-^9}:|:{:-^10}:|'.format('', '', ''))
        handle.write('\n')
        for n in range(RT_table.shape[0]):
            handle.write('|{: ^14}|{: ^11}|{: ^12}|'.format(RT_table["Mutation"].iloc[n], RT_table["Freq"].iloc[n], RT_table["Annotation"].iloc[n]))
            handle.write("\n")
        
        if comments:
            handle.write("\nComments:\n\n")
            for n in range(RT_table.shape[0]):
                if RT_table["Comment"].iloc[n] != "Unknown":
                    handle.write(f"- **{RT_table['Mutation'].iloc[n]}**: {RT_table['Comment'].iloc[n]}\n")    
    else:
        handle.write("\nNo **REVERSE TRANSCRIPTASE** mutations found.\n")


    if PR_table.shape[0] > 0:
        handle.write("\n**PROTEASE** mutations:\n\n")
        handle.write("|{: ^14}|{: ^11}|{: ^12}|\n".format('Mutation','Frequency', 'Annotation'))
        handle.write('|:{:-^12}:|:{:-^9}:|:{:-^10}:|'.format('', '', ''))
        handle.write('\n')
        for n in range(PR_table.shape[0]):
            handle.write('|{: ^14}|{: ^11}|{: ^12}|'.format(PR_table["Mutation"].iloc[n], PR_table["Freq"].iloc[n], PR_table["Annotation"].iloc[n]))
            handle.write("\n")
        
        if comments:
            handle.write("\nComments:\n\n")
            for n in range(PR_table.shape[0]):
                if PR_table["Comment"].iloc[n] != "Unknown":
                    handle.write(f"- **{PR_table['Mutation'].iloc[n]}**: {PR_table['Comment'].iloc[n]}\n")    
    else:
        handle.write("\nNo **PROTEASE** mutations found.\n")  
    
    


def write_ensemble_table(ensemble_results, sample_id, handle):
    '''Writes the ensemble results to a markdown file.
    INPUT:
    ensemble_results: List of dictionaries for each drug class containing ensemble predictions for each mutation list combination. Output from ensemble_predictions function.
    handle: file handle to write the results.
    '''
    cutoff = ensemble_results.iloc[0]['Cutoff']

    HIVDB, LSR, RF = False, False, False

    if "HIVDB" in ensemble_results.columns:
        HIVDB = True
    
    if "LSR" in ensemble_results.columns:
        LSR = True

    for cols in ensemble_results.columns:
        if cols == "RF":
            RF = True

    handle.write(" Mutations below " + str(cutoff) + " frequency were not included in the analysis. Susceptible results are shown in \n")
    handle.write(r"\textcolor{green}{green}")
    handle.write("\nand resistant results are shown in \n")
    handle.write(r"\textcolor{red}{red}.")
    handle.write("\nLinear regression and random forest models were trained on phenotypic data from PhenoSense assays.\n\n")
    
    for dataset in ["INI", "NNRTI", "NRTI", "PI"]:
        ensemble_results_dataset = ensemble_results[ensemble_results["Drug_Class"] == dataset]
        # print(ensemble_results_dataset["Mutations"])
        if 'No mutations' in ensemble_results_dataset["Mutations"].tolist(): #we skip the empty keys
            handle.write("\nNo mutations found for the " + dataset + " drug class.\n")
            continue

        handle.write("\n**" + dataset + "** predictions:\n\n")
        
        for mut_comb in ensemble_results_dataset["Mutations"].unique():

            handle.write("\n**" + mut_comb + "**\n\n")

            mut_df = ensemble_results_dataset[ensemble_results_dataset["Mutations"] == mut_comb]

            keep_cols, keep_color = ["Drug"], ["Drug"]
            if HIVDB:
                keep_cols.append("HIVDB_5_labels")
                keep_color.append("HIVDB")
            if LSR:
                mut_df["LSR_RF"] = mut_df["LSR_RF"].astype(float).round(4)
                keep_cols.append("LSR_RF")
                keep_color.append("LSR")
            if RF:
                keep_cols.append("RF")
                keep_color.append("RF")
            if HIVDB and LSR and RF:
                keep_cols.append("Ensemble")
                keep_color.append("Ensemble")

            mut_df_plot = mut_df[keep_cols]
            mut_df_color = mut_df[keep_color]
            mut_df_plot.columns = [col.replace("_5_labels", "").replace("LSR_RF", "Linear Regression").replace("RF", "Random Forest") for col in list(mut_df_plot.columns)]
            
            fig, ax = plt.subplots(figsize=(10, mut_df_plot.shape[0] * 0.25))
            ax.axis('off')
            
            table = ax.table(cellText=mut_df_plot.values, colLabels=mut_df_plot.columns, cellLoc='center', loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(12)
            table.scale(1.4, 1.4)
            table.auto_set_column_width([0,1])

            for i in range(mut_df_plot.shape[1]):
                table[(0, i)].set_fontsize(13)
                table[(0, i)].set_text_props(weight="bold")

            for i in range(mut_df_plot.shape[0]):
                table[(i+1, 0)].set_text_props(weight="bold")
                for j in range(mut_df_plot.shape[1]):
                    if mut_df_color.iloc[i, j] == "Resistant":
                        table[(i+1, j)].set_facecolor("red")
                        table[(i+1, j)].set_alpha(0.5)
                    elif mut_df_color.iloc[i, j] == "Susceptible":
                        table[(i+1, j)].set_facecolor("green")
                        table[(i+1, j)].set_alpha(0.3)
                    else:
                        table[(i+1, j)].set_facecolor("white")

            # s=f"/cluster/scratch/msanche/{sample_id}_{dataset}_{mut_comb.replace('+', '_')}.pdf"
            s=f"tmpdir_hiv/{sample_id}_{dataset}_{mut_comb.replace('+', '_')}.pdf"
            os.makedirs('tmpdir_hiv', exist_ok=True) #We save the plots in a tmpdir_hiv folder
            plt.tight_layout()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            plt.savefig(s, bbox_inches='tight', pad_inches = 0.05, transparent=True)
            plt.close(fig)
            handle.write(r"\begin{center}")
            handle.write("\n")
            handle.write(r"\includegraphics[width=\textwidth]")
            handle.write("{%s}\n" % s)
            handle.write("\end{center}\n")
    


def write_coverage_disclaimer(coverage_tsv, robustness_data, sample_id, handle):
    '''Writes the coverage disclaimer and writes a step plot showing the position of the sample's coverage.
    INPUT:
    coverage_tsv: tsv file with the coverage data, annotated with annotate_vcf.py. Similar to samtools depth.
    robustness_data: dictionary with the robustness data for each drug class needed for rhe step graph.
    handle: file handle to write the results.
    '''

    read_cutoff = coverage_tsv.iloc[0]['Read_Cutoff']

    handle.write("\n## Ensemble prediction coverage disclaimer\n\n")
    handle.write(f"The coverage of the sample was calculated from the FASTQ files. Only Drug Resistance Mutation (DRM) positions are analysed for robustness assessment. Positions were considered as covered above {read_cutoff} read(s).\n")

    for dataset in ["INI", "NNRTI", "NRTI", "PI"]:

        # if f"No mutations found for the {dataset} dataset." in handle: #We skip the coverage disclaimers from datasets without mutations
        #     handle.write(f"\nNo mutations found for the {dataset} dataset.\n")
        #     continue

        coverage_dataset = coverage_tsv[coverage_tsv["Drug_Class"] == dataset]
        sample_drm_coverage = coverage_dataset["DRM_Coverage"].values[0]
        handle.write(f"\n* {dataset} coverage\n")
        if sample_drm_coverage == 0:
            handle.write(f"\nNo DRM positions were covered for the {dataset} drug class.\n")
            continue

        robustness_step_plot(robustness_data, dataset, coverage_dataset, sample_id, handle) #we plot the robustness step plot
    
    handle.write("\n\n\*Balanced accuracy takes into account the ability to predict both resistant and susceptible mutations. In contrast with accuracy, it is calculated as the mean of sensitivity and specificity.\n")

def robustness_step_plot(robustness_data, dataset: str , coverage_tsv_dataset, sample_id, handle):
    '''Writes a step plot showing the balanced accuracy on the Y axis and the % of covered DRM positions on the X axis.
    We also plot with a red line the position of the analyzed sample coverage.
    INPUT:
    dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.
    coverage_value: coverage value of the analyzed sample as a %. i.e. 80% coverage = 80 .
    '''
    coverage_value = coverage_tsv_dataset["DRM_Coverage"].values[0]
    max_accuracy = coverage_tsv_dataset["Max_Balanced_Accuracy"].values[0]
    y_coverage = coverage_tsv_dataset["Estimated_Balanced_Accuracy"].values[0]
    if coverage_value > 100 or coverage_value < 0:
        raise ValueError(f'Coverage value should be between 0 and 100. {coverage_value} provided.')
    
    step_data = robustness_data[dataset]

    fig, axs = plt.subplots(1, 1, figsize=(10, 5))
    x = [0.5*m+1 for m in range(11)]
    # y = step_data.drop('Drug', axis=1).groupby(["Miss_muts"]).mean()
    y = []
    for mutations in step_data.keys():
        step_acc = []
        for drug in step_data[mutations].keys():
            step_acc.append(step_data[mutations][drug])
        y.append(sum(step_acc)/len(step_acc))
    
    #we get the accuracy for the coverage step
    # y_coverage = y.loc[coverage_step]["Balanced_Accuracy"]
    handle.write(f"\nThe estimated balanced accuracy for the {dataset} ensemble drug resistance prediction is **{round(y_coverage, 2)}** at {round(coverage_value, 2)}% DRM positions coverage. The reported balanced accuracy for a 100% coverage is {round(max_accuracy,2)}.\n\n")
    
    axs.axvline(x=((100-coverage_value)/10)*0.5 +1, color='red', linewidth = 2, label=f'Your sample ({round(coverage_value, 2)}% coverage)')#we plot a red line at the coverage value
    axs.step(x, y, where = 'post', alpha = 0.7, linewidth = 2, color='blue')
    axs.axhline(y=y_coverage, color='blue', linestyle='--', alpha = 0.35)
    axs.set_xticks([0.5*i+1 for i in range(11)])
    axs.set_xticklabels([100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0], fontsize=12)
    axs.set_ylim(0, 1)
    axs.set_xlabel('% of covered DRM positions', fontsize=15)
    axs.set_ylabel('Balanced Accuracy', fontsize=15)
    axs.set_title(f'Balanced accuracy for {dataset} class', fontsize=15)
    axs.legend(loc='upper right', fontsize=10)    #and we plot a small legend 

    # path_to_plot = f'/cluster/scratch/msanche/{sample_id}_robustness_{dataset}.pdf'
    os.makedirs('tmpdir_hiv', exist_ok=True) #We save the plots in a tmpdir_hiv folder
    path_to_plot = f'tmpdir_hiv/{sample_id}_robustness_{dataset}.pdf'
    plt.savefig(path_to_plot)

    handle.write(r"\begin{center}")
    handle.write("\n")
    handle.write(r"\includegraphics[width=\textwidth]")
    handle.write("{%s}\n" % path_to_plot)
    handle.write("\end{center}\n")

    return path_to_plot


def write_report_md(ensemble_tsv, coverage_tsv, single_mut_annotations_tsv, output_file, robustness_data, sample_id):
    ''' Writes the report in markdown format from the ensemble predictions, coverage data and single mutation annotations .tsv tables.'''
    # if not os.path.exists('HIV_resistance_reports'):
    #     os.makedirs('HIV_resistance_reports')

    md = open(output_file, 'w')
    md.write(f"# HIV-1 drug resistance report. Sample ID: {sample_id}\n")
    md.write("## Drug resistance prediction\n")
    # ensemble_predictions = ensemble_table(mut_tsv, higher_cutoff = higher_cutoff, HIVDB = HIVDB, LSR = LSR, RF = RF)
    write_ensemble_table(ensemble_tsv, sample_id, handle = md)
    write_coverage_disclaimer(coverage_tsv, robustness_data, sample_id, handle = md)
    md.write("\\newpage\n")    
    md.write("\n## Single mutation annotation\n")
    md.write("### HIVDB single mutation relevant annotations\n\n")
    write_HIVDB_table(single_mut_annotations_tsv, handle = md, unknown = False, comments = True)
    md.write("\\newpage\n")
    md.write("\n### All single mutations\n")
    write_HIVDB_table(single_mut_annotations_tsv, handle = md, unknown = True, comments = False)

def main(ensemble_tsv, coverage_tsv, single_mut_annotations_tsv, output_file, robustness_data, sample_id):
    ''' Writes the report in markdown format from the ensemble predictions, coverage data and single mutation annotations .tsv tables to .md from samples.'''

    ensemble_preds_table = pd.read_csv(ensemble_tsv, sep="\t")#, dtype={"Sample_ID":str})
    ensemble_preds_table.fillna('No mutations', inplace=True)
    coverage_data = pd.read_csv(coverage_tsv, sep="\t")#, dtype={"Sample_ID":str})
    single_mut_annotations = pd.read_csv(single_mut_annotations_tsv, sep="\t")#, dtype={"Sample_ID":str})
    
    write_report_md(ensemble_preds_table, coverage_data, single_mut_annotations, output_file, robustness_data, sample_id)