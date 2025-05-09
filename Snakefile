"""
Input to the workflow is a .vcf and coverage.tsv.
Output is a report on HIV drug resistance prediction of the sample.

1. Annotate .vcf and coverage.tsv bp files to Pol gene proteins of interest.
2. Predict HIV drug resistance and generate a .pdf report on these results.

"""

# rule definitions
rule all:
    input:
        "snakemake_trial/mutation_freq.tsv",
        "snakemake_trial/coverage_annotated.tsv",
        "snakemake_trial/report.pdf"


## 1. Annotated vcf and coverage.tsv files with aminoacids
rule annotate_vcf:
    input:
        fname_vcf="example_files/CAP257/week_80/variants/snvs.vcf",
        ref="ref/reference_gagpol_only.gb",
        fname_cov="example_files/CAP257/week_80/alignments/coverage.tsv"
    output:
        fname_vcf="snakemake_trial/mutation_freq.tsv",
        fname_cov="snakemake_trial/coverage_annotated.tsv"
    conda:
        "env/annotate_vcf_env.yml"
    script:
        "annotate_vcf.py"

## 2. Drig resistance prediction report
rule HIV_resistance:
    input:
        fname_vcf="snakemake_trial/mutation_freq.tsv",
        fname_cov="snakemake_trial/coverage_annotated.tsv"
    output:
        fname_report="snakemake_trial/report.pdf"
    conda:
        "env/resistance_pred_env.yml"
    script:
        "resistance_prediction.py"


