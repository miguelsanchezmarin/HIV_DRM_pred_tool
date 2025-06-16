###Labels a list of mutations as Resistant or Susceptible for different anti-HIV drugs

import typer
from typing_extensions import Annotated, Optional


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

    print("Single sample detected. Running in serial...")
    sample_outputs = predict_annotate_and_coverage(sample_data, coverage_data, hivdb_data, lsr_coeff, rf_models, robustness_data, 0.15, 0.03, HIVDB, LSR, RF)
    ensemble_df = get_ensemble_df(sample_outputs[0])
    annotations_df = get_single_mut_annotations_df(sample_outputs[1])
    coverage_df = get_coverage_df(sample_outputs[2])
    #We save the dataframes to tsv files
    ensemble_df.to_csv(output_ensemble_tsv, sep="\t", index=False)
    annotations_df.to_csv(output_single_tsv, sep="\t", index=False)
    coverage_df.to_csv(output_coverage_tsv, sep="\t", index=False)
   
    print(f"Saved HIV resistance files to {'/'.join(output_ensemble_tsv.split('/')[:-1])} directory.")
