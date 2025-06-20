# Comparative analysis and integration of HIV drug resistance mutation testing and prediction tools
This repository contains the data and code needed to reproduce the figures and results from the Master project report held in the Computational Biology Group (D-BSSE, ETH Zürich, Basel) by Miguel Sánchez Marín (MSc Bioinformatics and Biocomplexity, Utrecht University).

![img1](figures/method_performance_f1_score_barplot.png)

## File description
* **HIV_prediction_models/**: contains the scripts used for training and prediction of the different HIV drug resistance prediction methods compared.
    - **cnn_models/**: contains the trained CNN models for the different drugs and folds. 
        - ```{drug}_model_fold_{fold}.keras```: CNN model for drug 'drug' and trained on all folds but number 'fold'.
    - ```HIVDB_prediction.py```: script used to predict the benchmark dataset with HIVDB program (via [sierrapy client](https://github.com/hivdb/sierra-client.git)). Outputs HIVDB predictions in ```method_predictions/HIVDB```.
    - ```cnn_prediction.R```: script used to predict the benchmark dataset with the CNN implementation from [Steiner et al.](https://doi.org/10.3390/v12050560) . Outputs CNN predictions in ```method_predictions/cnn```.
    - ```cnn_training.R```: script used to train the CNN models following the implementation from [Steiner et al.](https://doi.org/10.3390/v12050560) . Trained models stored in ```HIV_prediction_models/cnn_models```.
    - ```ensemble_prediction.py```: script used to retrieve the benchmark dataset predictions from the ensemble approaches explained in the report. Ensemble predictions are stored in ```method_predictions/small_ensemble_predictions.tsv``` and ```method_predictions/large_ensemble_predictions.tsv```.
    - ```linear_regression_prediction.R```: script used to predict the benchmark dataset with the linear regression approximations explained in the report (LSR, LARS, LSR-I and LARS-I). The script is adapted from the original implementation by [Rhee et al.](https://doi.org/10.1073/pnas.0607274103). Linear regression predictions are stored in ```method_predictions/linear_regression```.
    - ```random_forest_prediction.py```: script used to predict the benchmark dataset following the random forest implementation by [Raposo et al.](https://doi.org/10.1007/978-3-030-38021-2_6). Random forest predictions are stored in ```method_predictions/random_forest```.
* **datasets/**: contains the datasets used for the method comparison and robustness analysis. Drug susceptibility and mutation data were all obtained from High-quality filtered genotype-phenotype datasets from the public [Stanford HIV drug resistance database](https://hivdb.stanford.edu/pages/genopheno.dataset.html). 
    - **{drug_class}/**: contains independent drug-specific datasets including the fold distribution.
        - ```{drug_class}_{drug}_5folds.tsv```
    - **fasta/**: contains the drug-specific dataset in FASTA format, as required for the CNN implementation.
        - ```{drug}.fasta```
    - **sumsampled_datasets/**: randomly subsampled datasets from the original datasets. Data for increasingly larger missing sequence positions is used for sequence coverage robustness analysis. Obtained with ```utilities/subsample_dataset.py```.
        - ```{drug_class}_subsampled_dataset.tsv```
    - ```{drug_class}_dataset.tsv```: original Stanford HIV drug resistance database datasets after prefiltering done for the report.
    - ```dataset_fold_distribution.tsv```: summary with the sizes of the different folds for every drug dataset and the folds distribution in Resistant and Susceptible labels. Obtained with ```dataset_distribution_plot.py```. (Supplementary Table 1)
* **figures/**: contains the figure images from the report.
    - **roc_pr_curves/**: contains the ROC and PR AUC graphs.
        - ```{drug_class}_roc_pr_curves_5folds_combined.png```: Supplementary Figures 5-8.
    - ```{performance_metric}_robustness_percent_plot.png```: Figure 3 and Supplementary Figures 9-11.
    - ```{performance_metric}_robustness_total_number_plot.png```: Supplementary Figures 12-15.
    - ```dataset_distribution_plot.png```: Figure 1.
    - ```method_performance_{performance_metric}_barplot.png```: Figure 2 and extra figures for accuracy, balanced accuracy and MCC-
    - ```method_performance_{performance_metric}_heatmap.png```: Supplementary Figures 1-4.
    - ```scaling_analysis_plot.png```: Figure 4. 
* **method_predictions/**: contains the benchmark single predictions for the different methods and their performance metrics.
    - **HIVDB/**: contains HIVDB predictions, generated with ```HIV_prediction_models/HIVDB_prediction.py```.
        - ```HIVDB_{drug_class}_predictions.tsv```
    - **cnn/**: contains CNN predictions, generated with ```HIV_prediction_models/cnn_prediction.R```.
        - ```cnn_{drug}_fold_{fold}_predictions.tsv```
    - **ensemble_subsampled_predictions/**: contains the *Small ensemble* predictions for the subsampled dataset.
        - **{drug_class}/**
            - ```{drug}_ensemble_predictions.tsv```
    - **linear_regression/**: contains the linear regression predictions for LSR, LARS, LSR-I and LARS-I, obtained with ```HIV_prediction_models/linear_regression_prediction.R```.
        - **{drug_class}/**
            - **{drug}/**
                - ```LARS_{drug}_fold_{fold}_predictions.tsv```: LARS predictions for drug 'drug' and fold number 'fold'.
                - ```LARS_I_{drug}_fold_{fold}_predictions.tsv```: LARS-I predictions for drug 'drug' and fold number 'fold'.
                - ```LSR_{drug}_fold_{fold}_predictions.tsv```: LSR predictions for drug 'drug' and fold number 'fold'.
                - ```LSR_I_{drug}_fold_{fold}_predictions.tsv```: LSR-I predictions for drug 'drug' and fold number 'fold'.
    - **random_forest/**: contains the random forest predictions, obtained with ```HIV_prediction_models/random_forest_prediction.py```.
        - **{drug_class}/**
            - **{drug}/**
                - ```random_forest_{drug}_fold_{fold}_predictions.tsv```: random forest predictions for drug 'drug' and fold number 'fold'.
    - ```HIVDB_performance_5folds.tsv```: performance metrics for HIVDB predictions on the benchmark.
    - ```LARS_I_performance_5folds.tsv```: performance metrics for LARS-I predictions on the benchmark.
    - ```LARS_performance_5folds.tsv```: performance metrics for LARS predictions on the benchmark.
    - ```LSR_I_performance_5folds.tsv```: performance metrics for LSR-I predictions on the benchmark.
    - ```LSR_performance_5folds.tsv```: performance metrics for LSR predictions on the benchmark.
    - ```cnn_performance_5folds.tsv```: performance metrics for CNN predictions on the benchmark.
    - ```large_ensemble_performance_5folds.tsv```: performance metrics for *Large ensemble* predictions on the benchmark.
    - ```large_ensemble_predictions.tsv```: *Large ensemble* predictions on the benchmark, generated with ```HIV_prediction_models/ensemble_prediction.py```.
    - ```random_forest_performance_5folds.tsv```: performance metrics for random forest predictions on the benchmark.
    - ```roc_pr_auc_data.tsv```: probability data for LSR, LSR-I, LARS, LARS-I and CNN for ROC and PR AUC calculation and plotting.
    - ```small_ensemble_performance_5folds.tsv```: performance metrics for *Small ensemble* predictions on the benchmark.
    - ```small_ensemble_predictions.tsv```: *Small ensemble* predictions on the benchmark, generated with ```HIV_prediction_models/ensemble_prediction.py```.
* **utilities/**:
    - ```dataset_fold_creation.py```: scipt used for the 5 fold creation from the original dataset. Takes as input ```datasets/{drug_class}_dataset.tsv``` and generates ```datasets/{drug_class}/{drug_class}_{drug}_5folds.tsv```.
    - ```subsample_dataset.py```: generates a subsampled dataset from the original dataset, with different % of sequence coverage, used for robustness to sequence coverage analysis. Takes as input ```datasets/{drug_class}_dataset.tsv``` and generates ```datasets/subsampled_datasets/{drug_class}_subsampled_dataset.tsv```.
- ```dataset_distribution_plot.py```: creates dataset distribution plots ```figures/dataset_distribution_plot.png``` and fold distribution table ```datasets/dataset_fold_distribution.tsv```.
- ```method_performance_plots.py```: creates comparison plots for the different method performances: ```figures/method_performance_{performance_metric}_barplot.png``` and ```figures/method_performance_{performance_metric}_heatmap.png```.
- ```roc_pr_auc_plot.py```: creates ROC and PR AUC plots in ```figures/roc_pr_curves/*```.
- ```scaling_analysis_data.tsv```: data for the scaling analysis of the software tool. 
- ```scaling_analysis_plot.py```: creates a plot with the scaling analysis data: ```figures/scaling_analysis_plot.png```.
- ```seq_coverage_robustness_plot.py```: generates robustness to sequence coverage performance plots: ```figures/{performance_metric}_robustness_percent_plot.png``` and ```figures/{performance_metric}_robustness_total_number_plot.png```.
    



