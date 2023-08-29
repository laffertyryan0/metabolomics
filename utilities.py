import json
import numpy as np

# Basic utilities
def get_date_str():
    return datetime.datetime.now().isoformat().replace(':',"_").replace(".","_").replace("-","_")


# Initial parameter manipulation
def shift_b(amt=-.1):
    init = json.load(open("./data/init.json","r"))
    init["true_b"] = [p + amt for p in init["true_b"]]
    json.dump(init,open("./data/init.json","w"))

def shift_samplenum(amt=50):
    init = json.load(open("./data/init.json","r"))
    init["samples_per_lab"] = [p + amt for p in init["samples_per_lab"]]
    json.dump(init,open("./data/init.json","w"))

def get_init_dict():
    init = json.load(open("./data/init.json","r"))
    return init

def get_init_corr():
    return np.array(get_init_dict()['corr'])
