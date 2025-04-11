# HIV_DRM_pred_tool
Drug resistance prediction tool for HIV NGS .vcf files.

## File description
* HIVDB_rules/
  - HIVDB_rules/{drug_class}_muts_score_Stanford_HIVDB : .tsv files with single mutation scores for each drug class from the HIVDB program (https://hivdb.stanford.edu/dr-summary/mut-scores/PI/). Files downloaded from StanfordHIV Data Base on 11th April 2025.
  - HIVDB_rules/{drug_class}_combinations_score_Stanford_HIVDB : .tsv files with single mutation scores for each drug class from the HIVDB program (https://hivdb.stanford.edu/dr-summary/mut-scores/PI/). Files downloaded from StanfordHIV Data Base on 11th April 2025.
* linear_regression_coefficients/
  - linear_regression_coefficients/OLS_{drug}_combinations_tsm_all_folds.txt : for each drug, files with coefficients for the different features (single and combined mutations) comming from Ordinary Least Squares Regression implementation based on Rhee et al. (https://doi.org/10.1073/pnas.0607274103).
* random_forest_models/
  - random_forest_models/random_forest_python_{drug}_RF_model_allfolds.pkl : for each drug, Random Forest models trained on all the project dataset. The Random Forest implementation for HIV drug resistance prediction is based on the implementation by Raposo et al. (https://doi.org/10.1007/978-3-030-38021-2_6)
