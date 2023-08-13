import subprocess
import numpy as np
import os
from scipy.io import loadmat
import pickle
import sys
import shutil
import datetime

def run(num_iter, save_folder = "."):
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

    matlab_code = "outputs = {};"


    for i in range(1,num_iter+1):
        #call r script, saves x_matrix.csv and h_matrix.csv
        subprocess.run("/usr/bin/Rscript get_x_and_h.R",
                        shell=True)
        X = np.loadtxt("x_matrix.csv",delimiter=",")
        H = np.loadtxt("h_matrix.csv",delimiter=",")
        
        shutil.move("./x_matrix.csv","./"+save_folder+"/aux/x_matrix.csv")
        shutil.move("./h_matrix.csv","./"+save_folder+"/aux/h_matrix.csv")
        
        x_str = np.array2string(X,precision=9,separator=",",threshold=np.inf)\
                .replace("],\n","];")\
                .replace("\n","...\n")
        h_str = np.array2string(H,precision=9,separator=",",threshold=np.inf)\
                .replace("],\n","];")\
                .replace("\n","...\n")

        matlab_code = matlab_code + "\n"
        matlab_code = matlab_code + f"X = {x_str};\n"
        matlab_code = matlab_code + "for i=1:length(X);X(i,i) = 1;end\n"
        matlab_code = matlab_code + f"H = {h_str};\n"
        matlab_code = matlab_code + "out = nearcorr(X,'Weights',H,'MaxIterations',10000);\n"
        matlab_code = matlab_code + "outputs{"+str(i)+"} = out;\n"
        
    matlab_code = matlab_code + "\nsave('outputs.mat','outputs');"
    if os.path.exists("simulation.m"):
        os.remove("simulation.m")
                  
    with open("simulation.m","w") as file:
        file.write(matlab_code)

    print("A")
    subprocess.run("bash matlab_call.sh",shell=True)
    print("B")

    shutil.move("./simulation.m","./"+save_folder+"/aux/simulation.m")

    matfile = loadmat('outputs.mat')
    outputs = list(matfile['outputs'].flatten('A'))

    with open("outputs.pkl","wb") as file:
        pickle.dump(outputs, file)

    
    shutil.move("./outputs.mat","./"+save_folder+"/outputs.mat")
    shutil.move("./outputs.pkl","./"+save_folder+"/outputs.pkl")
    shutil.move("./missing_vecs.csv","./"+save_folder+"/missing_vecs.csv")
    shutil.move("./pr_vecs.csv","./"+save_folder+"/pr_vecs.csv")
    shutil.move("./min_eigs.csv","./"+save_folder+"/min_eigs.csv")

if __name__ == "__main__":
    num_iter = 10
    if not os.path.exists("outputs"):
        os.mkdir("./outputs")
    save_folder =  "./outputs/"+datetime.datetime.now().isoformat().replace(':',"_").replace(".","_").replace("-","_")

    for i,arg in enumerate(sys.argv):
        if arg == "-i":
            num_iter = int(sys.argv[i+1])

    run(num_iter,save_folder)
