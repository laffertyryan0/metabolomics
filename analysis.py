from run_simulation import run
from generate_params import generate

import utilities

import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# Named analyses
    
def monte_carlo_analysis_24aug2023():
    #Make a single monte-carlo heatmap for (b, samplenum) given
    
    generate(seed=None)
    out_path = run(num_iter=10)
    with open(f"{out_path}/outputs.pkl","rb") as file:
        est_corr_matrices = pickle.load(file)
        
    true_corr_matrix = utilities.get_init_corr()
    bias_matrices = [est_corr_matrix-true_corr_matrix for est_corr_matrix in est_corr_matrices]
    avg_biases = sum(bias_matrices)/len(bias_matrices)

    sns.set()
    plt.figure(figsize=(8,6))
    sns.heatmap(avg_biases,annot=False,cmap="Blues",fmt=".2f",linewidth=0.5)
    plt.title('Title here')
    plt.show()


# Main entry point
    
if __name__ == "__main__":
    monte_carlo_analysis_24aug2023()
    pass


