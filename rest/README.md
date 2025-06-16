# HIV_DRM_pred_tool
Drug resistance prediction tool for HIV NGS .vcf files.

## Running instructions

0. (Variant calling should be done with ref/reference.fa gag-pol sequence, not done in this tool. Indexed on ref/reference.fa.fai ).
1. The current tool uses as input a directory with all samples data organised in directories with their sample ID as names. You should organise your input data the following way:
``` bash  
|--input_dir/  
|           |-- sample_ID_1/  
|           |              |-- coverage.tsv  
|           |              |--snvs.vcf  
|           |-- sample_ID_2/  
|           |              |-- coverage.tsv  
|           |              |-- snvs.vcf  
|           ...  
|   
```  

2. input_dir path should be specified on the Snakefile sample_dir variable (line 14).
3. If snvs.vcf or coverage.tsv files have a different naming, specify it on the Snakefile annotate_vcf rule input variables fname_vcf and fname_cov respectively (lines 31 and 32).
4. If **you have a proper Snakemake and mamba installation** then you can run the predictions and report generations with:
``` snakemake --cores 4 --use-conda ```
5. If **you don't want to use --use-conda option** (as I did). You can create a conda environment from the env/snakemake_all.yml file:
``` conda env create -f env/snakemake_all.yml ```.
Activate the environment:
``` conda activate snakemake_all ```.
And run:
``` snakemake --cores 4  ```.



## File description
* **env/**
   - env/annotate_vcf_env.yml : requirements file for annotate_vcf environment used in Snakefile for annotate_vcf rule.
   - env/resistance_pred_env.yml : requirements file for HIV_pred environment used in Snakefile for latter rules.
   - env/snakemake_all.yml : requirements file of both environments above. If active, Snakemake can be run without --use-conda option.
* **HIVDB_rules/**
    - HIVDB_rules/{drug_class}_combinations_score_Stanford_HIVDB : .tsv files with single mutation scores for each drug class from the HIVDB program (https://hivdb.stanford.edu/dr-summary/mut-scores/PI/). Files downloaded from StanfordHIV Data Base on 11th April 2025.
    - HIVDB_rules/{drug_class}_comments_Stanford_HIVDB : .tsv files with single mutation annotations and comments for each drug class from the HIVDB program (https://hivdb.stanford.edu/dr-summary/comments/PI/). Files downloaded from StanfordHIV Data Base on 11th April 2025.
    - HIVDB_rules/{drug_class}_muts_score_Stanford_HIVDB : .tsv files with single mutation scores for each drug class from the HIVDB program (https://hivdb.stanford.edu/dr-summary/mut-scores/PI/). Files downloaded from StanfordHIV Data Base on 11th April 2025.
* **linear_regression_coefficients/**
    - linear_regression_coefficients/OLS_{drug}_combinations_tsm_all_folds.txt : for each drug, files with coefficients for the different features (single and combined mutations) comming from Ordinary Least Squares Regression implementation based on Rhee et al. (https://doi.org/10.1073/pnas.0607274103).
* **random_forest_models/**
    - random_forest_models/random_forest_python_{drug}_RF_model_allfolds.pkl : for each drug, Random Forest models trained on all the project dataset. The Random Forest implementation for HIV drug resistance prediction is based on the implementation by Raposo et al. (https://doi.org/10.1007/978-3-030-38021-2_6).
* **ref/**
   - ref/GCF_000864765.1_ViralProj15476_genomic.gbff : GenBank file for the HIV genome reference.
   - ref/reference_gagpol_only.gb : GenBank file only with the GagPol gene from the HIV genome. This is the REFERENCE used for annotation by annotate_vcf.
* **resources/**
   - resources/template.tex : LaTex template used for final .pdf report generation. Copied from MinVar GitHub (https://github.com/medvir/MinVar/tree/master).
* **robustness_data/**
   - robustness_data/{drug_class}_step_data.tsv: for each drug class, balanced accuracy data in relation to DRM position coverage % (on 10% steps).
* **annotate_vcf.py** : annotates .vcf and coverage.tsv files from base pairs to aminoacid positions using a genbank reference (ref/reference_gagpol_only.gb). Adapted from Lara's script (https://github.com/LaraFuhrmann/Scan-for-mutations-of-interest-NGS-samples/blob/main/workflow/scripts/annotate_vcf.py)
* **coverage_analysis.py** : functions necessary for reading coverage robustness data. 
* **drug_resistance_ensemble_prediction.py** : functions necessary for drug resistance prediction with the ensemble method (HIVDB + LSR-I+ Random Forest) in resistance_prediction.py .
* **hivdb_single_mut_annotation.py** : functions necessary for single mutation annotation with HIVDB comments (from HIVDB_rules/{drug_class}_comments_Stanford_HIVDB). 
* **load_data.py** : functions necessary for loading input data (sample and coverage) and HIVDB, LSR-I, Random Forest and robustness data into dictionaries for the assessment generation.
* **md_to_pdf.py** : converts .md files to .pdf report using pandoc package with shell commands. Current Snakefile does not use it. 
*  **resistance_prediction.py** : predicts drug resistance for all drugs with the ensemble method.
* **Snakefile** : Snakemake pipeline that takes as input a directory with per sample .vcf and coverage.tsv files and ouputs:
   - snakemake_trial/output_vcfs/{sample_ID}/coverage_annotated.tsv : sample annotated aminoacid position coverage.tsv.
   - snakemake_trial/output_vcfs/{sample_ID}/coverage_disclaimer.tsv : sample per drug class estimated balanced accuracy based on sample coverage.
   - snakemake_trial/output_vcfs/{sample_ID}/ensemble_predictions.tsv : sample per drug drug resistance prediction from the ensemble and from the different methods in the ensemble for each mutation combination.
   - snakemake_trial/output_vcfs/{sample_ID}/report_{sample_ID}.pdf: final report assessing drug resistance on the sample.
   - snakemake_trial/output_vcfs/{sample_ID}/mut_freq.tsv : sample mutations per position on the different proteins (PR, RT and IN), coming from .vcf file.
   - snakemake_trial/output_vcfs/{sample_ID}/single_mut_annotations.tsv : sample single mutations annotated with HIVDB comments.
* **write_report_md.py** : writes the assessment report in .md format from .tsv per sample files. 

## Scaling analysis

<div align="center">

| Samples (n) | Time (4 cores) | Time (32 cores) |
|:-----------:|:--------------:|:---------------:|
|     1       |      22s       |      20s        |
|     20      |     1m 45s     |      29s        |
|     50      |     4m 21s     |      52s        |
|     500     |     54m 46s    |      7m 50s     |
|    1000     |    1h 36m 56s  |      12m 18s    |

</div>
Runs were done withouth --use-conda option, with the snakemake_all environment already active.
(1GB RAM / core)

## *Running annotation for all subtypes
Standard and tested annotation is only done for subtype B samples using the whole gag-pol sequence (see instructions above). Annotation for all subtypes was developed but not tested, available to use at your own risk. In this case, variant calling should be done with just the pol sequence from the defined subtype (available in ref/pol_refs). Then, the Snakefile rule 1 should be silenced and the rule 1B should be uncommented and the reference file original used set on the 'ref' variable (line 46). The subtype is automatically retrieved from the FASTA reference. Annotation is done with annotate_pol.py script.