"""
Input to the workflow is a .vcf and coverage.tsv.
Output is a report on HIV drug resistance prediction of the sample.

1. Annotate .vcf and coverage.tsv bp files to Pol gene proteins of interest.
2. Predict HIV drug resistance and generate a .pdf report on these results.

"""
import os
sample_dir = "input_vcfs" #directory with a subdirectory for each sample with the sampleID as the name and the .vcf and coverage.tsv files inside.
sample_list = os.listdir(sample_dir)
sample_list = [sample_list[0]]
# sample_list = ["CAP30", "CAP54", "CAP80", "CAP135", "CAP161", "CAP191"]

# rule definitions
rule all:
    input:
        # "HIV_resistance_reports/report_{samples}.pdf"
        # "HIV_resistance_reports/report_000001.pdf",
        # "HIV_resistance_reports/report_000002.pdf",
        # "HIV_resistance_reports/report_000003.pdf",
        # "HIV_resistance_reports/report_000004.pdf",
        # "HIV_resistance_reports/report_000005.pdf"
        # expand("snakemake_trial/output_vcfs/{samples}/mut_freq.tsv", samples=["CAP30", "CAP54", "CAP80", "CAP135", "CAP161", "CAP191"]),
        # expand("snakemake_trial/output_vcfs/{samples}/coverage_annotated.tsv", samples=["CAP30", "CAP54", "CAP80", "CAP135", "CAP161", "CAP191"])
        # expand("snakemake_trial/output_vcfs/{samples}/coverage_disclaimer.tsv", samples=sample_list),
        # expand("snakemake_trial/output_vcfs/{samples}/ensemble_predictions.tsv", samples=sample_list),
        # expand("snakemake_trial/output_vcfs/{samples}/single_mut_annotations.tsv", samples=sample_list),
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
        # dir_output = "output_vcfs"
        # 
    output:
        fname_vcf="snakemake_trial/output_vcfs/{samples}/mut_freq.tsv",
        fname_cov="snakemake_trial/output_vcfs/{samples}/coverage_annotated.tsv"
        # fname_cov="snakemake_trial/coverage_annotated.tsv"
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
        
        """
        #rm -r tmpdir_hiv/{wildcards.samples}*
# os.rmdir("tmpdir_hiv") #we remove the temporal directory we created
