import matplotlib.pyplot as plt

from drug_resistance_ensemble_prediction import *
from hivdb_single_mut_annotation import *
from coverage_analysis import *

def write_HIVDB_table(HIVDB_table, cutoff, handle, unknown: bool = False, comments: bool = False):
    '''Writes the ensemble results to a markdown file.
    INPUT:
    HIVDB_table: Dataframe with the single mutations, frequency and HIVDB comments. Output from HIVDB_table function.
    cutoff: frequency cutoff to consider the mutations.
    handle: file handle to write the results.
    unknown: if True, the unknown annotations will be written.
    comments: if True, the comments will be written.
    '''

    handle.write("Single mutation annotations obtained from HIVDB program.")
    handle.write("Mutations below " + str(cutoff) + " frequency were not included in the analysis.\n")
    IN_table, RT_table, PR_table = HIVDB_table[HIVDB_table["Prot"]=="IN"], HIVDB_table[HIVDB_table["Prot"]=="RT"], HIVDB_table[HIVDB_table["Prot"]=="PR"]
    if not unknown:
        IN_table, RT_table, PR_table = IN_table[IN_table["Annotation"] != "Unknown"], RT_table[RT_table["Annotation"] != "Unknown"], PR_table[PR_table["Annotation"] != "Unknown"]

    if IN_table.shape[0] > 0:
        handle.write("**INTEGRASE** mutations:\n")
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
    
    


def write_ensemble_table(ensemble_results, cutoff, handle, HIVDB = True, LSR = True, RF = True):
    '''Writes the ensemble results to a markdown file.
    INPUT:
    ensemble_results: List of dictionaries for each drug class containing ensemble predictions for each mutation list combination. Output from ensemble_predictions function.
    cutoff: cutoff for the mutations.
    handle: file handle to write the results.
    '''
    handle.write("Mutations below " + str(cutoff) + " frequency were not included in the analysis. Susceptible results are shown in \n")
    handle.write(r"\textcolor{green}{green}")
    handle.write("\nand resistant results are shown in \n")
    handle.write(r"\textcolor{red}{red}.")
    handle.write("\n\n")
    
    for i, drug_class_dic in enumerate(ensemble_results):

        if '' in drug_class_dic.keys(): #we skip the empty keys
            continue

        if i == 0:
            handle.write("**INI** predictions:\n")
            dataset = "INI"
        elif i == 1:
            handle.write("\n**NNRTI** predictions:\n")
            dataset = "NNRTI"
        elif i == 2:
            handle.write("\n**NRTI** predictions:\n")
            dataset = "NRTI"
        elif i == 3:
            handle.write("\n**PI** predictions:\n")
            dataset = "PI"
        
        for mut_comb in drug_class_dic.keys():

            handle.write("\n**" + mut_comb.replace(", ", "+") + "**\n\n")

            mut_df = drug_class_dic[mut_comb]
            mut_df["Drug"] = mut_df.index
            # print(mut_df)

            keep_cols, keep_color = ["Drug"], ["Drug"]
            if HIVDB:
                keep_cols.append("HIVDB_five_labels")
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
            # columns = mut_df_plot.columns
            mut_df_plot.columns = [col.replace("_five_labels", "").replace("LSR_RF", "Linear Regression").replace("RF", "Random Forest") for col in list(mut_df_plot.columns)]
            
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

            s=f"example_files/{dataset}_{mut_comb.replace(', ', '_')}.pdf"
            plt.tight_layout()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            plt.savefig(s, bbox_inches='tight', pad_inches = 0.05, transparent=True)
            plt.close(fig)
            handle.write(r"\begin{center}")
            handle.write("\n")
            handle.write(r"\includegraphics[width=\textwidth]")
            handle.write("{%s}\n" % s)
            handle.write("\end{center}\n")
    


def write_coverage_disclaimer(coverage_tsv, handle, read_cutoff:int = 1):
    '''Writes the coverage disclaimer and writes a step plot showing the position of the sample's coverage.
    INPUT:
    coverage_tsv: tsv file with the coverage data, annotated with annotate_vcf.py. Similar to samtools depth.
    read_cutoff: minimum number of reads to consider a position as covered.
    handle: file handle to write the results.
    '''
    handle.write("\n## Coverage disclaimer\n\n")
    handle.write(f"The coverage of the sample was calculated from the FASTQ files. Only Drug Resistance Mutation (DRM) positions are analysed for robustness assessment. Positions were considered as covered above {read_cutoff} read(s).\n")

    coverage_data = pd.read_csv(coverage_tsv, sep='\t')
    coverage_data.columns = ['Subtype', 'Position', 'Reads', 'Protein']

    for dataset in ["INI", "NNRTI", "NRTI", "PI"]:
        if dataset == "INI":
            drm_positions = [51, 66, 74, 75, 92, 95, 97, 118, 121, 122, 138, 140, 143, 145, 146, 147, 148, 151, 153, 155, 157, 163, 230, 232, 263]
            prot = "IN"
        elif dataset == "NNRTI":
            drm_positions = [90, 98, 100, 101, 103, 106, 108, 138, 179, 181, 188, 190, 221, 225, 227, 230, 234, 236, 238, 318, 348]
            prot = "RT"
        elif dataset == "NRTI":
            drm_positions = [41, 62, 65, 67, 68, 69, 70, 74, 75, 77, 115, 116, 151, 184, 210, 215, 219]
            prot = "RT"
        elif dataset == "PI":
            drm_positions = [10, 20, 24, 32, 33, 46, 47, 48, 50, 53, 54, 73, 74, 76, 82, 83, 84, 88, 89, 90]
            prot = "PR"

        coverage_data_prot = coverage_data[coverage_data['Protein'] == prot]
        coverage_data_drm = coverage_data_prot[coverage_data_prot['Position'].isin(drm_positions)]
        #we calculate the % of covered positions
        sample_drm_coverage = coverage_data_drm[coverage_data_drm['Reads'] > read_cutoff].shape[0] / len(drm_positions) * 100

        handle.write(f"\n* {dataset} coverage\n")
        if sample_drm_coverage == 0:
            handle.write(f"\nNo DRM positions were covered for the {dataset} dataset.\n")
            continue

        robustness_step_plot(dataset, sample_drm_coverage, handle) #we plot the robustness step plot
    
    handle.write("\n\n\*Balanced accuracy takes into account the ability to predict both resistant and susceptible mutations. In contrast with accuracy, it is calculated as the mean of sensitivity and specificity.\n")

def write_report_md(mut_tsv, coverage_tsv, higher_cutoff: float = 0.15, lower_cutoff: float = 0.015, HIVDB:bool = True, LSR: bool = True, RF: bool = True):
    ''' Writes the report in markdown format.'''
    HIVDB_single_table = HIVDB_table(mut_tsv, lower_cutoff = lower_cutoff)

    md = open('example_files/report.md', 'w')
    md.write("# HIV-1 drug resistance report\n")
    md.write("## Drug resistance prediction\n")
    print("Starting the HIV drug resistance ensemble predictions...")
    ensemble_predictions = ensemble_table(mut_tsv, higher_cutoff = higher_cutoff, HIVDB = HIVDB, LSR = LSR, RF = RF)
    print("HIV drug resistance ensemble predictions done.")
    write_ensemble_table(ensemble_predictions, cutoff = higher_cutoff, handle = md, HIVDB = HIVDB, LSR = LSR, RF = RF)
    print("HIV drug resistance ensemble predictions written to file.")
    write_coverage_disclaimer(coverage_tsv, handle = md)
    print("Coverage disclaimer written to file.")
    md.write("\\newpage\n")    
    md.write("\n## Single mutation annotation\n")
    md.write("### HIVDB single mutation relevant annotations\n\n")
    write_HIVDB_table(HIVDB_single_table, cutoff = lower_cutoff, handle = md, unknown = False, comments = True)
    md.write("\\newpage\n")
    md.write("\n### All single mutations\n")
    write_HIVDB_table(HIVDB_single_table, cutoff = lower_cutoff, handle = md, unknown = True, comments = False)
    print("HIVDB single mutations annotations written to file.")