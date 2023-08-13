import json
import numpy as np
import shutil
import os
import sys

def generate(command_args = None):
    np.random.seed(5)

    N = 50
    P = 2
    K = 50

    for i, arg in enumerate(command_args):
        if arg == "-N" or arg == "-n":
            N = int(command_args[i+1])
        if arg == "-P" or arg == "-p":
            P = int(command_args[i+1])
        if arg == "-K" or arg == "-k":
            K = int(command_args[i+1])

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
    init_dict["true_b"] = (2*np.random.randn(N)/3).tolist()

    init_dict["samples_per_lab"] = [int(x) for x in (40 + 2.5*np.random.rand(N))]
    init_dict["mean"] = np.zeros(K).tolist()
    #init_dict["cov"] = ((np.ones((K,K))+np.diag(np.ones(K)))/2).tolist() interclass correlation
    mat = np.zeros((K,K))
    mat[:int(K/3),:int(K/3)] += 1
    mat[int(K/3):int(2*K/3+1),int(K/3):int(2*K/3+1)] += 1
    mat[int(2*K/3):,int(2*K/3):] += 1
    mat[mat>1]=1
    mat /= 2
    mat = mat + np.diag(np.ones(K))/2
    print(mat)
    print(np.real(np.linalg.eig(mat)[0]))
    
    init_dict["cov"] = mat.tolist()
    


    with open(init_fn,"w") as file:
        json.dump(init_dict,file)

if __name__ == "__main__":
    generate(sys.argv)
