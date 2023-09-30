from run_simulation import run
from generate_params import generate_params

import utilities

import pickle
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# Named analyses
    
def monte_carlo_analysis_24aug2023():
    #Make a single monte-carlo heatmap for (b, samplenum) given
    
    generate_params(seed=None,
             num_lab_covariates = 2,
             avg_samples_per_lab = 40,
             num_metabolites = 1000)
    init_dict = utilities.get_init_dict()
    
    out_path = run(num_iter=1)
    with open(f"{out_path}/outputs.pkl","rb") as file:
        est_corr_matrices = pickle.load(file)
        
    output_folder_name = out_path.split('/')[-1]
    results_path = f"./analysis_results/monte_carlo_analysis_24aug2023___{output_folder_name}"
    os.mkdir(results_path)

    #plot the average of the biases
    
    true_corr_matrix = utilities.get_init_corr()
    bias_matrices = [est_corr_matrix-true_corr_matrix for est_corr_matrix in est_corr_matrices]
    avg_biases = sum(bias_matrices)/len(bias_matrices)
    
    sns.set()
    plt.figure(figsize=(8,6))

    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    
    sns.heatmap(avg_biases,annot=False,cmap=cmap,center=0,fmt=".2f",linewidth=0.5)
    plt.title(f'Average Intercept = {np.mean(init_dict["true_b"])}, Average samples per lab = {np.mean(init_dict["samples_per_lab"])}')
    
    plt.savefig(results_path + "/bias_heatmap.png")

    #also plot the average of just the estimates
    
    avg_ests = sum(est_corr_matrices)/len(est_corr_matrices)
    sns.set()
    plt.figure(figsize=(8,6))

    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    
    sns.heatmap(avg_ests,annot=False,cmap=cmap,center=0,fmt=".2f",linewidth=0.5)
    plt.title(f'Average Intercept = {np.mean(init_dict["true_b"])}, Average samples per lab = {np.mean(init_dict["samples_per_lab"])}')
    
    plt.savefig(results_path + "/average_estimated_correlation_heatmap.png")

def monte_carlo_analysis_7sept2023():
    #Make a single monte-carlo heatmap for (b, samplenum) given
    #use imshow instead of seaborn
    
    generate_params(seed=None,
             num_lab_covariates = 2,
             avg_samples_per_lab = 40,
             num_metabolites = 1000)
    init_dict = utilities.get_init_dict()
    
    out_path = run(num_iter=5)
    with open(f"{out_path}/outputs.pkl","rb") as file:
        est_corr_matrices = pickle.load(file)
        
    output_folder_name = out_path.split('/')[-1]
    results_path = f"./analysis_results/monte_carlo_analysis_7sept2023___{output_folder_name}"
    os.mkdir(results_path)

    #plot the average of the biases
    
    true_corr_matrix = utilities.get_init_corr()
    bias_matrices = [est_corr_matrix-true_corr_matrix for est_corr_matrix in est_corr_matrices]
    avg_biases = sum(bias_matrices)/len(bias_matrices)
    
    plt.figure(figsize=(8,6))
    
    plt.imshow(avg_biases)
    plt.title(f'Average Intercept = {np.mean(init_dict["true_b"])}, Average samples per lab = {np.mean(init_dict["samples_per_lab"])}')
    
    plt.savefig(results_path + "/bias_heatmap.png")

    #also plot the average of just the estimates
    
    avg_ests = sum(est_corr_matrices)/len(est_corr_matrices)

    plt.figure(figsize=(8,6))
    
    plt.imshow(avg_ests)
    plt.title(f'Average Intercept = {np.mean(init_dict["true_b"])}, Average samples per lab = {np.mean(init_dict["samples_per_lab"])}')
    
    plt.savefig(results_path + "/average_estimated_correlation_heatmap.png")

def compare_real_to_generated_data_29aug2023():
    #generate data based on estimated correlation matrices, and compare with "stacked_met_data" i.e. "true data"
    num_iter = 2
        
    generate_params(seed=None,
             num_lab_covariates = 2,
             avg_samples_per_lab = 40,
             num_metabolites=1000)
    init_dict = utilities.get_init_dict()
    
    out_path = run(num_iter=num_iter)
    with open(f"{out_path}/outputs.pkl","rb") as file:
        est_corr_matrices = pickle.load(file)

    stacked_met_data = []
    samples_per_lab = init_dict['samples_per_lab']

    #retrieve data for each iteration
    for i in range(num_iter):
        iteration_data = []
        stacked_data = np.loadtxt(f"{out_path}/aux/stacked_met_data{i}.csv",delimiter=",")
        for lab in range(len(samples_per_lab)):
            iteration_data.append(stacked_data[samples_per_lab[lab]:(samples_per_lab[lab+1]-1)])
        stacked_met_data.append(iteration_data)
            
    
        
    output_folder_name = out_path.split('/')[-1]
    results_path = f"./analysis_results/compare_real_to_generated_data_29aug2023___{output_folder_name}"
    os.mkdir(results_path)

# Main entry point
    
if __name__ == "__main__":
    monte_carlo_analysis_7sept2023()
    #compare_real_to_generated_data_29aug2023()
    pass


