###Labels a list of mutations as Resistant or Susceptible for different anti-HIV drugs

import os
import typer
from typing_extensions import Annotated, Optional
from multiprocessing import Pool
from tqdm import tqdm

from load_data import *
from drug_resistance_ensemble_prediction import *
from hivdb_single_mut_annotation import *
from coverage_analysis import *
from write_report_md import *


def predict_annotate_and_coverage(mut_dict:dict, coverage_dict: dict, hivdb_data: dict, lsr_coeff: dict, rf_models: dict, robustness_data: dict, higher_cutoff: float, lower_cutoff: float, HIVDB: bool, LSR: bool, RF: bool):
    '''Generates all three outputs from the input data: ensemble resistance predictions, single mutation annotations and coverage data disclaimer.'''

    ensemble_preds = ensemble_predictions(mut_dict, hivdb_data, lsr_coeff, rf_models, higher_cutoff, HIVDB, LSR, RF)
    single_mut_annotations = HIVDB_single_mut(mut_dict, hivdb_data, lower_cutoff)
    coverage_disclaimer = get_coverage_dict(coverage_dict, robustness_data)
    
    return [ensemble_preds, single_mut_annotations, coverage_disclaimer]

# ###Main function
# def main(input_mut: Annotated[str ,typer.Argument(help="Path to a mutation_freq.tsv annotated file (output from annotate_vcf.py).")],
#          input_coverage: Annotated[str ,typer.Argument(help="Path to a coverage_annotated.tsv annotated file (output from annotate_vcf.py).")],
#          output: Annotated[Optional[str], typer.Option('--output', '-o', help='Output file name.')] = "report.pdf",
#             HIVDB: Annotated[Optional[bool], typer.Option('--HIVDB', '-H', help='Use HIVDB method for prediction')] = None,
#             LSR: Annotated[Optional[bool], typer.Option('--LSR', '-L', help='Use Linear Regression method for prediction')] = None,
#             RF: Annotated[Optional[bool], typer.Option('--RF', '-R', help='Use Random Forest method for prediction')] = None):
    
#     # Example usage
#     if HIVDB == None and LSR == None and RF == None:
#         write_report_md(input_mut, input_coverage, HIVDB = True, LSR = True, RF = True)
#     else:
#         if HIVDB == None:
#             HIVDB = False
#         if LSR == None:
#             LSR = False
#         if RF == None:
#             RF = False

#         write_report_md(input_mut, input_coverage, HIVDB = HIVDB, LSR = LSR, RF = RF)
#     report_to_pdf('example_files/report.md', output)

###Main function
def main(input_mut: Annotated[str ,typer.Argument(help="Path to a directory organised in per sample directories. Each sample directory must have a mutation_freq.tsv and coverage_annotated.tsv file.")],
         coverage_file: Annotated[str ,typer.Argument(help="Path to a directory organised in per sample directories. Each sample directory must have a mutation_freq.tsv and coverage_annotated.tsv file.")],
         hivdb_data: dict,
         lsr_coeff: dict,
         rf_models: dict,
         robustness_data: dict,
         output_ensemble_tsv: Annotated[Optional[str], typer.Option('--output', '-o', help='Path to directory where output files will be stored.')] = "HIV_prediction_output",
         output_single_tsv: Annotated[Optional[str], typer.Option('--output', '-o', help='Path to directory where output files will be stored.')] = "HIV_prediction_output",
         output_coverage_tsv: Annotated[Optional[str], typer.Option('--output', '-o', help='Path to directory where output files will be stored.')] = "HIV_prediction_output",
            HIVDB: Annotated[Optional[bool], typer.Option('--HIVDB', '-H', help='Use HIVDB method for prediction')] = None,
            LSR: Annotated[Optional[bool], typer.Option('--LSR', '-L', help='Use Linear Regression method for prediction')] = None,
            RF: Annotated[Optional[bool], typer.Option('--RF', '-R', help='Use Random Forest method for prediction')] = None):
    
    ##We convert the flags to boolean
    if HIVDB == None and LSR == None and RF == None:
        HIVDB = True
        LSR = True
        RF = True
    else:
        if HIVDB == None:
            HIVDB = False
        if LSR == None:
            LSR = False
        if RF == None:
            RF = False
    
    #We load all the samples data in a dictionary
    sample_data = load_sample_data(input_mut)
    coverage_data = load_coverage_data(coverage_file)
    # print(f"Loaded {len(sample_data)} sample(s).")

    # if HIVDB:        #We load the HIVDB data
    #     hivdb_data = load_hivdb_data()

    # if LSR:         #We load the linear regression coefficients in a dictionary
    #     lsr_coeff = load_lsr_coef()

    # if RF:          #We load the random forest models in a dictionary
    #     rf_models = load_random_forest()  
    
    # print("Loaded all models.")

    # robustness_data = load_robustness_data()    #We load the robustness data
    # print("Loaded robustness data.")

    ###8.546s it took to load everything for just one sample

    # if len(sample_data) > 1:
    #     print("Multiple samples detected. Running in parallel...")
    #     #We create a list of all the samples
    #     sample_keys = list(sample_data.keys())
    #     num_cores = 4

    #     pool = Pool(num_cores)
    #     print(f"Predicting {len(sample_keys)} samples with {num_cores} cores...")

    #     sample_outputs = pool.starmap(predict_annotate_and_coverage, tqdm([(sample_data[sample], coverage_data[sample], hivdb_data, lsr_coeff, rf_models, robustness_data, 0.15, 0.015, HIVDB, LSR, RF) for sample in sample_keys]))
    #     # ensemble_preds = pool.starmap(ensemble_predictions, tqdm([(sample_data[sample], hivdb_data, lsr_coeff, rf_models, 0.15, HIVDB, LSR, RF) for sample in sample_keys]))
    #     print(sample_outputs)
    #     pool.close()
    # else:
    print("Single sample detected. Running in serial...")
    sample_outputs = predict_annotate_and_coverage(sample_data, coverage_data, hivdb_data, lsr_coeff, rf_models, robustness_data, 0.15, 0.015, HIVDB, LSR, RF)
    ensemble_df = get_ensemble_df(sample_outputs[0])
    annotations_df = get_single_mut_annotations_df(sample_outputs[1])
    coverage_df = get_coverage_df(sample_outputs[2])
    #We save the dataframes to tsv files


    # ensemble_df.to_csv(snakemake.output.fname_vcf, sep="\t", index=False)
    # annotations_df.to_csv(snakemake.output.fname_single, sep="\t", index=False)
    # coverage_df.to_csv(snakemake.output.fname_cov, sep="\t", index=False)
    ensemble_df.to_csv(output_ensemble_tsv, sep="\t", index=False)
    annotations_df.to_csv(output_single_tsv, sep="\t", index=False)
    coverage_df.to_csv(output_coverage_tsv, sep="\t", index=False)

    # print(ensemble_preds)
    # print(f"Annotating {len(sample_keys)} samples with {num_cores} cores...")
    # single_mut_annotations = pool.starmap(HIVDB_single_mut, tqdm([(sample_data[sample], hivdb_data, 0.015) for sample in sample_keys]))
    # print(single_mut_annotations)
    
    # print(f"Getting coverage data disclaimer for {len(sample_keys)} samples with {num_cores} cores...")
    # coverage_disclaimer = pool.starmap(get_coverage_dict, tqdm([(coverage_data[sample], robustness_data) for sample in sample_keys]))
    # print(coverage_disclaimer)

    # coverage_disclaimer, single_mut_annotations = [], []
    # for sample_id in sample_keys:
    #     single_mut_annotations.append(HIVDB_single_mut(sample_data[sample_id], hivdb_data, 0.015))
    #     coverage_disclaimer.append(get_coverage_dict(coverage_data[sample_id], robustness_data))
    
    


    # pool.close()
    # ensemble_preds = [sample_prediction[0] for sample_prediction in sample_outputs]
    # single_mut_annotations = [sample_prediction[1] for sample_prediction in sample_outputs]
    # coverage_disclaimer = [sample_prediction[2] for sample_prediction in sample_outputs]
    
    # ensemble_df = get_ensemble_df(ensemble_preds, sample_keys)
    # print(ensemble_df)
    # annotations_df = get_single_mut_annotations_df(single_mut_annotations, sample_keys)
    # print(annotations_df)
    # coverage_df = get_coverage_df(coverage_disclaimer, sample_keys)
    # print(coverage_df)

    # #We create the output directory if it doesn't exist
    # if not os.path.exists(output):
    #     os.makedirs(output)
    #     print(f"Created directory {output}")
    
    # #We save the dataframes to tsv files
    # ensemble_df.to_csv(f"{output}/ensemble_predictions.tsv", sep="\t", index=False)
    # annotations_df.to_csv(f"{output}/single_mut_annotations.tsv", sep="\t", index=False)
    # coverage_df.to_csv(f"{output}/coverage_disclaimer.tsv", sep="\t", index=False)
    print(f"Saved HIV resistance files to {'/'.join(output_ensemble_tsv.split('/')[:-1])} directory.")



    

# # input_mut = "K20R, V35I, T39A, M41L, S68G, K103N, I135T, M184V, T200K, Q207E, T215Y" 
# input_mut = "input_samples"
# if __name__ == "__main__":
#     # main("input_samples")
#     main(
#         snakemake.input.fname_vcf,
#         snakemake.input.fname_cov
#     )
    # typer.run(main("input_samples"))
#     main(
#         snakemake.input.fname_vcf,
#         snakemake.input.fname_cov,
#         snakemake.output.fname_report
#     )

# table = pd.read_csv('example_files/mutation_freq.tsv', sep='\t')
# table = table[table['Freq'] > 0.05]
# print(extract_prot_mut_freq(table, 'RT'))
# input_list = [["K20R", "V35I", "T39A", "M41L", "S68G", "K103N", "I135T", "M184V", "T200K", "Q207E", "T215Y"], ["K20M", "V35I", "T39A", "M41L", "S68G", "K103N", "I135T", "M184V", "T200K", "Q207E", "T215Y"]]
# print(get_mutlist_comb(input_list))
# print(ensemble_table('example_files/mutation_freq.tsv', higher_cutoff = 0.05, HIVDB = True, LSR = True, RF = True))
# print(unify_mut_list(input_list))
# write_report_md('example_files/mutation_freq.tsv', 'example_files/CAP257/week_54/alignments/coverage.tsv', higher_cutoff = 0.15, lower_cutoff = 0.015, HIVDB = True, LSR = True, RF = True)
# robustness_step_plot("NNRTI", 78)
# main("input_samples")