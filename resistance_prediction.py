###Labels a list of mutations as Resistant or Susceptible for different anti-HIV drugs
# import pandas as pd
# import re
import os
# import pickle
import typer
from typing_extensions import Annotated, Optional
from itertools import product
import matplotlib.pyplot as plt

from load_data import *
from drug_resistance_ensemble_prediction import *
from hivdb_single_mut_annotation import *
from coverage_analysis import *
from write_report_md import *


def report_to_pdf(path_to_md, path_to_pdf):
    ''' Converts the markdown report to pdf using pandoc.'''
    os.system(f"pandoc {path_to_md} -o {path_to_pdf} --pdf-engine=tectonic --template=resources/template.tex")
    print(f"Report saved at {path_to_pdf}")

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
         output: Annotated[Optional[str], typer.Option('--output', '-o', help='Output file name.')] = "report.pdf",
            HIVDB: Annotated[Optional[bool], typer.Option('--HIVDB', '-H', help='Use HIVDB method for prediction')] = None,
            LSR: Annotated[Optional[bool], typer.Option('--LSR', '-L', help='Use Linear Regression method for prediction')] = None,
            RF: Annotated[Optional[bool], typer.Option('--RF', '-R', help='Use Random Forest method for prediction')] = None):
    
    # Example usage
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
    print(f"Loaded {len(sample_data)} sample(s).")

    if HIVDB:        #We load the HIVDB data
        hivdb_data = load_hivdb_data()

    if LSR:         #We load the linear regression coefficients in a dictionary
        lsr_coeff = load_lsr_coef()

    if RF:          #We load the random forest models in a dictionary
        rf_models = load_random_forest()  
    
    print("Loaded all models.")
    ###6.827s it took to load everything for just one sample
    

# input_mut = "K20R, V35I, T39A, M41L, S68G, K103N, I135T, M184V, T200K, Q207E, T215Y" 

# if __name__ == "__main__":
#     # typer.run(main)
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
main("input_samples")