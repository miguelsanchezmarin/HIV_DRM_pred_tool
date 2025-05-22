"""
Input to the workflow is a .vcf and coverage.tsv.
Output is a report on HIV drug resistance prediction of the sample and these data on .tsv machine readable tables.

1. Annotate .vcf and coverage.tsv bp files to Gag-Pol gene proteins of interest.
2. Predict HIV drug resistance, estimate the balanced accuracy of the predictions based on the sample coverage and annotate single mutations with HIVDB comments.
2b. Aggregate all samples into a single machine readable file. One for each of the three tasks.
3. Generate a per sample .md report on the results of the analysis.
4. Convert the .md report to .pdf

"""
import os
sample_dir = "input_vcfs" #directory with a subdirectory for each sample with the sampleID as the name and the .vcf and coverage.tsv files inside.
sample_list = os.listdir(sample_dir)
sample_list = sample_list[:6]
# sample_list = ["CAP30", "CAP54", "CAP80", "CAP135", "CAP161", "CAP191"]

# rule definitions
rule all:
    input:
        "snakemake_trial/output_vcfs/ensemble_predictions.tsv",
        "snakemake_trial/output_vcfs/coverage_disclaimer.tsv",
        "snakemake_trial/output_vcfs/single_mut_annotations.tsv",
        expand("snakemake_trial/output_vcfs/{samples}/report_{samples}.pdf", samples=sample_list)


## 1. Annotated vcf and coverage.tsv files with aminoacids
rule annotate_vcf:
    input:
        fname_vcf= sample_dir + "/{samples}/mix_12_variants_chromchange.vcf",
        fname_cov= sample_dir + "/{samples}/coverage.tsv",
        ref="ref/reference_gagpol_only.gb",
    output:
        fname_vcf="snakemake_trial/output_vcfs/{samples}/mut_freq.tsv",
        fname_cov="snakemake_trial/output_vcfs/{samples}/coverage_annotated.tsv"
    conda:
        "env/annotate_vcf_env.yml"
    script:
        "annotate_vcf.py"

## 2. Drug resistance prediction to machine readable format
from resistance_prediction import load_hivdb_data, load_lsr_coef, load_random_forest, load_robustness_data #we first load the necessary data for the ensemble imputation and coverage disclaimer

hivdb_data_all = load_hivdb_data()
lsr_coeff_all = load_lsr_coef()
rf_models_all = load_random_forest()
robustness_data_all = load_robustness_data()

rule HIV_resistance:
    input:
        fname_vcf="snakemake_trial/output_vcfs/{samples}/mut_freq.tsv",
        fname_cov="snakemake_trial/output_vcfs/{samples}/coverage_annotated.tsv"
    output:
        fname_vcf="snakemake_trial/output_vcfs/{samples}/ensemble_predictions.tsv",
        fname_cov="snakemake_trial/output_vcfs/{samples}/coverage_disclaimer.tsv",
        fname_single="snakemake_trial/output_vcfs/{samples}/single_mut_annotations.tsv"
    run:
        from resistance_prediction import main

        main(
            input_mut=f"{input.fname_vcf}",
            coverage_file=f"{input.fname_cov}",
            hivdb_data=hivdb_data_all,
            lsr_coeff=lsr_coeff_all,
            rf_models=rf_models_all,
            robustness_data=robustness_data_all,
            output_ensemble_tsv=f"{output.fname_vcf}",
            output_single_tsv=f"{output.fname_single}",
            output_coverage_tsv=f"{output.fname_cov}",
            )

        # "resistance_prediction.py"

## 2B. Aggregate all samples into a single machine readable file.

rule aggregate_all_samples:
    input:
        fname_vcf=expand("snakemake_trial/output_vcfs/{samples}/ensemble_predictions.tsv", samples=sample_list),
        fname_cov=expand("snakemake_trial/output_vcfs/{samples}/coverage_disclaimer.tsv", samples=sample_list),
        fname_single=expand("snakemake_trial/output_vcfs/{samples}/single_mut_annotations.tsv", samples=sample_list)
    output:
        fname_vcf="snakemake_trial/output_vcfs/ensemble_predictions.tsv",
        fname_cov="snakemake_trial/output_vcfs/coverage_disclaimer.tsv",
        fname_single="snakemake_trial/output_vcfs/single_mut_annotations.tsv"
    run:
        import pandas as pd     #we append the tsv files of all samples into a single file adding a SeqID column

        for vcf, cov, single, sample in zip(input.fname_vcf, input.fname_cov, input.fname_single, sample_list):
            df_vcf, df_cov, df_single = pd.read_csv(vcf, sep="\t"), pd.read_csv(cov, sep="\t"), pd.read_csv(single, sep="\t")
            df_vcf["SeqID"] = sample
            df_cov["SeqID"] = sample
            df_single["SeqID"] = sample

            if sample == sample_list[0]:
                df_vcf_all = df_vcf
                df_cov_all = df_cov
                df_single_all = df_single
            else:
                df_vcf_all = pd.concat([df_vcf_all, df_vcf], ignore_index=True)
                df_cov_all = pd.concat([df_cov_all, df_cov], ignore_index=True)
                df_single_all = pd.concat([df_single_all, df_single], ignore_index=True)
        df_vcf_all.to_csv(output.fname_vcf, sep="\t", index=False)
        df_cov_all.to_csv(output.fname_cov, sep="\t", index=False)
        df_single_all.to_csv(output.fname_single, sep="\t", index=False)
            

## 3. Generate an .md report on the results
rule generate_md_report:
    input:
        fname_vcf="snakemake_trial/output_vcfs/{samples}/ensemble_predictions.tsv",
        fname_cov="snakemake_trial/output_vcfs/{samples}/coverage_disclaimer.tsv",
        fname_single="snakemake_trial/output_vcfs/{samples}/single_mut_annotations.tsv"
    output:
        report="snakemake_trial/output_vcfs/{samples}/report_{samples}.md"
    run:
        from write_report_md import main

        main(
            ensemble_tsv=f"{input.fname_vcf}",
            coverage_tsv=f"{input.fname_cov}",
            single_mut_annotations_tsv=f"{input.fname_single}",
            output_file =f"{output.report}",
            robustness_data=robustness_data_all,
            sample_id=f"{wildcards.samples}"
        )


## 4. Convert the .md report to .pdf
rule convert_md_to_pdf:
    input:
        report="snakemake_trial/output_vcfs/{samples}/report_{samples}.md"
    output:
        report="snakemake_trial/output_vcfs/{samples}/report_{samples}.pdf"
    conda:
        "env/resistance_pred_env.yml"
    shell:
        """
        pandoc {input.report} -o {output.report} --pdf-engine=tectonic --template=resources/template.tex
        rm {input.report}
        rm -rf tmpdir_hiv/{wildcards.samples}_PI_* tmpdir_hiv/{wildcards.samples}_INI_* tmpdir_hiv/{wildcards.samples}_NRTI_* tmpdir_hiv/{wildcards.samples}_NNRTI_* tmpdir_hiv/{wildcards.samples}_robustness_*
        """
        


# if os.path.exists("tmpdir_hiv") and os.path.isdir("tmpdir_hiv"):
#     if len(os.listdir("tmpdir_hiv")) == 0: #if the directoty is already empty we delete it
#         os.rmdir("tmpdir_hiv")
    

