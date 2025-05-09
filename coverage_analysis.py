import pandas as pd

import matplotlib.pyplot as plt

def robustness_step_plot(dataset: str , coverage_value: float, handle):
    '''Writes a step plot showing the balanced accuracy on the Y axis and the % of covered DRM positions on the X axis.
    We also plot with a red line the position of the analyzed sample coverage.
    INPUT:
    dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.
    coverage_value: coverage value of the analyzed sample as a %. i.e. 80% coverage = 80 .
    '''
    if coverage_value > 100 or coverage_value < 0:
        raise ValueError(f'Coverage value should be between 0 and 100. {coverage_value} provided.')
    
    if dataset == "INI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    elif dataset == "NNRTI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    elif dataset == "NRTI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    elif dataset == "PI":
        step_data = pd.read_csv(f"robustness_data/{dataset}_step_data.tsv", sep="\t")
    else:
        raise ValueError(f'Unknown dataset: {dataset}')

    fig, axs = plt.subplots(1, 1, figsize=(10, 5))
    x = [0.5*m+1 for m in range(11)]
    y = step_data.drop('Drug', axis=1).groupby(["Miss_muts"]).mean()
    max_accuracy = y["Balanced_Accuracy"].max()
    coverage_step = round(1-((int(coverage_value/10)+1)/10),5)
    if coverage_value == 100:
        coverage_step = 0
    if coverage_value == 0:
        coverage_step = 1
    #we get the accuracy for the coverage step
    y_coverage = y.loc[coverage_step]["Balanced_Accuracy"]
    # print(f"Coverage step: {coverage_step}, y_coverage: {y_coverage}, max_accuracy: {max_accuracy}, coverage_value: {coverage_value}")
    handle.write(f"\nThe estimated balanced accuracy for the {dataset} dataset drug resistance prediction is **{round(y_coverage, 2)}** at {round(coverage_value, 2)}% DRM positions coverage. The reported balanced accuracy for a 100% coverage is {round(max_accuracy,2)}.\n\n")
    
    axs.axvline(x=((100-coverage_value)/10)*0.5 +1, color='red', linewidth = 2, label=f'Your sample ({round(coverage_value, 2)}% coverage)')#we plot a red line at the coverage value
    axs.step(x, y["Balanced_Accuracy"], where = 'post', alpha = 0.7, linewidth = 2, color='blue')
    axs.axhline(y=y_coverage, color='blue', linestyle='--', alpha = 0.35)
    axs.set_xticks([0.5*i+1 for i in range(11)])
    axs.set_xticklabels([100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0], fontsize=12)
    axs.set_ylim(0, 1)
    axs.set_xlabel('% of covered DRM positions', fontsize=15)
    axs.set_ylabel('Balanced Accuracy', fontsize=15)
    axs.set_title(f'Balanced accuracy for {dataset} class', fontsize=15)
    axs.legend(loc='upper right', fontsize=10)    #and we plot a small legend 

    path_to_plot = f'example_files/robustness_{dataset}.pdf'
    plt.savefig(path_to_plot)

    handle.write(r"\begin{center}")
    handle.write("\n")
    handle.write(r"\includegraphics[width=\textwidth]")
    handle.write("{%s}\n" % path_to_plot)
    handle.write("\end{center}\n")

    return f'example_files/coverage_robustness_{dataset}.pdf'

