import subprocess
import numpy as np
import os
from scipy.io import loadmat
import pickle
import sys
import shutil
import datetime


"""
This assumes that certain parameters are unknown but fixed:
These include the true correlation matrix, covariates, the parameters of the
of the missingness probit random effects model.
To alter the fixed parameters, use generate_params.py, or change them manually.
Non-fixed parameters include the metabolite data itself, which is generated
from a multivariate distribution with the given correlation matrix, as well as
the missingness distribution. 
Every time the run function is called, the correlation matrix will be estimated
num_iter times, on different generated datasets with the same true correlation and
other fixed parameters.
If the estimator is consistent, the averages of these matrices should converge
to the true correlation.
The resulting matrices will not be averaged automatically but will be stored separately
in the cells of a cell array, so that there are num_iter cells. Also, a python pkl version
of this will be generated. 
"""

def run(num_iter = None, save_folder = None):
    
    if not os.path.exists("outputs"):
        os.mkdir("./outputs")

    if num_iter is None:
        num_iter = 10
    if save_folder is None:
        save_folder =  "./outputs/"+datetime.datetime.now().isoformat().replace(':',"_").replace(".","_").replace("-","_")
    
    if os.path.exists("pr_vecs.csv"):
        os.remove("pr_vecs.csv")
    if os.path.exists("missing_vecs.csv"):
        os.remove("missing_vecs.csv")
    if os.path.exists("min_eigs.csv"):
        os.remove("min_eigs.csv")

    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
        
    if not os.path.exists(save_folder+"/aux"):
        os.mkdir(save_folder+"/aux")

    matlab_code = "outputs = {};\ndrops = {};"


    for i in range(1,num_iter+1):
        #call r script, saves x_matrix.csv and h_matrix.csv
        subprocess.run("/usr/bin/Rscript get_x_and_h.R",
                        shell=True)
        X = np.loadtxt("x_matrix.csv",delimiter=",")
        H = np.loadtxt("h_matrix.csv",delimiter=",")
        drops = np.loadtxt("./mets_in_some_lab_index.csv")
        
        
        shutil.move("./stacked_met_data.csv","./"+save_folder+f"/aux/stacked_met_data{i}.csv")
        shutil.move("./x_matrix.csv","./"+save_folder+f"/aux/x_matrix{i}.csv")
        x_path = "./"+save_folder+f"/aux/x_matrix{i}.csv"
        shutil.move("./h_matrix.csv","./"+save_folder+f"/aux/h_matrix{i}.csv")
        h_path = "./"+save_folder+f"/aux/h_matrix{i}.csv"
        shutil.move("./missing_vecs.csv","./"+save_folder+f"/aux/missing_vecs{i}.csv")
        shutil.move("./pr_vecs.csv","./"+save_folder+f"/aux/pr_vecs{i}.csv")
        shutil.move("./min_eigs.csv","./"+save_folder+f"/aux/min_eigs{i}.csv")
        shutil.move("./mets_in_some_lab_index.csv","./"+save_folder+f"/aux/mets_in_some_lab_index{i}.csv")
        
        x_str = np.array2string(X,precision=9,separator=",",threshold=np.inf)\
                .replace("],\n","];")\
                .replace("\n","...\n")
        h_str = np.array2string(H,precision=9,separator=",",threshold=np.inf)\
                .replace("],\n","];")\
                .replace("\n","...\n")
        drops_str = np.array2string(drops,precision=9,separator=",",threshold=np.inf)\
                .replace("],\n","];")\
                .replace("\n","...\n")

        matlab_code = matlab_code + "\n"
        matlab_code = matlab_code + f"X = readmatrix('{x_path}');\n"
        matlab_code = matlab_code + "for i=1:length(X);X(i,i) = 1;end\n"
        matlab_code = matlab_code + f"H = readmatrix('{h_path}');\n"
        matlab_code = matlab_code + "drops{"+str(i)+"} = "+drops_str+";\n"
        matlab_code = matlab_code + "disp('start nc')\n"
        matlab_code = matlab_code + "out = nearcorr(X,'Weights',H,'MaxIterations',10000);\n"
        matlab_code = matlab_code + "disp('stop nc')\n"
        matlab_code = matlab_code + "outputs{"+str(i)+"} = out;\n"

        matlab_code = matlab_code + "disp('hi')\n"
        
    matlab_code = matlab_code + "\nfinal_drops = ones(1,length(drops{1}));for i=1:length(drops);final_drops=final_drops.*drops{i};end;"
    matlab_code = matlab_code + "\nfor(i = 1:length(outputs));"
    matlab_code = matlab_code + "\nnew_out = repmat(NaN,length(final_drops),length(final_drops));"
    matlab_code = matlab_code + "\nnew_out(drops{i}'*drops{i}==1)=outputs{i};"
    matlab_code = matlab_code + "\noutputs_final{i} = new_out(final_drops==1,final_drops==1);end;"
    matlab_code = matlab_code + "\noutputs = outputs_final;"

    matlab_code = matlab_code + "\nsave('outputs.mat','outputs');"
    
    if os.path.exists("simulation.m"):
        os.remove("simulation.m")
                  
    with open("simulation.m","w") as file:
        file.write(matlab_code)

    subprocess.run("bash matlab_call.sh",shell=True)

    shutil.move("./simulation.m","./"+save_folder+"/aux/simulation.m")

    matfile = loadmat('outputs.mat')
    outputs = list(matfile['outputs'].flatten('A'))

    with open("outputs.pkl","wb") as file:
        pickle.dump(outputs, file)

    
    shutil.move("./outputs.mat","./"+save_folder+"/outputs.mat")
    shutil.move("./outputs.pkl","./"+save_folder+"/outputs.pkl")
    shutil.copytree("./data","./"+save_folder+f"/aux/data")

    return save_folder

if __name__ == "__main__":
    num_iter = None
    save_folder = None

    for i,arg in enumerate(sys.argv):
        if arg == "-i":
            num_iter = int(sys.argv[i+1])

    run(num_iter,save_folder)
