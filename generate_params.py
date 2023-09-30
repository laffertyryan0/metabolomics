import json
import numpy as np
import shutil
import os
import sys

def generate_params(num_labs=50, #num labs
             num_lab_covariates=2, #num lab covariates
             num_metabolites=50, #num metabolites
             avg_samples_per_lab = 40,
             samples_per_lab_deviation = 2.5,
             seed = None
             ):
    if seed is not None:
        np.random.seed(seed)

    N = num_labs 
    P = num_lab_covariates  
    K = num_metabolites 



    init_fn = "./data/init.json"

    init_dict = dict()
    init_dict["covariate_files"] = '~/projects/metabolomics/data/covariates'

    #delete and regenerate the covariate files
    try:
        shutil.rmtree("./data/covariates")
    except:
        pass
    os.mkdir("./data/covariates")
    for i in range(K):
        np.savetxt(fname=f"./data/covariates/covariates_metabolite#{i}.csv",
                   X = np.random.randn(N,P),
                   delimiter=",")

    init_dict["true_beta"] = np.random.randn(P).tolist()
    init_dict["true_b"] = (2*np.random.randn(N)/3 + 2).tolist() 

    init_dict["samples_per_lab"] = [max(0,int(x)) for x in
                                    (avg_samples_per_lab +
                                     samples_per_lab_deviation*(np.random.randn(N)))]
    init_dict["mean"] = np.zeros(K).tolist()
    #init_dict["cov"] = ((np.ones((K,K))+np.diag(np.ones(K)))/2).tolist() interclass correlation
    mat = np.zeros((K,K))
    mat[:int(K/3),:int(K/3)] += 1
    mat[int(K/3):int(2*K/3+1),int(K/3):int(2*K/3+1)] += 1
    mat[int(2*K/3):,int(2*K/3):] += 1
    mat[mat>1]=1
    mat /= 2
    mat = mat + np.diag(np.ones(K))/2
    
    init_dict["corr"] = mat.tolist()
    


    with open(init_fn,"w") as file:
        json.dump(init_dict,file)

if __name__ == "__main__":

    command_args = sys.argv

    N = 17
    K = 12
    P = 2

    for i, arg in enumerate(command_args):
        if arg == "-N" or arg == "-n":
            N = int(command_args[i+1])
        if arg == "-P" or arg == "-p":
            P = int(command_args[i+1])
        if arg == "-K" or arg == "-k":
            K = int(command_args[i+1])
            
    generate_params(N,P,K,seed=42)
