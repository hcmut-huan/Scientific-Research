import sys
import json
from enum import Enum

class Error(Enum):
    INVALID_ALGO=1
    MISSING_PARAM=2
    INVALID_PARAM=3
    

def validate(algo):
    if algo["name"] == "pmc":
        if "error" not in algo or "mode" not in algo:
            exit(Error.MISSING_PARAM.value)
        elif algo["mode"] != "midrange" and algo["mode"] != "mean":
            exit(Error.INVALID_PARAM.value)
            
    elif algo["name"] == "hybrid-pca":
        if "error" not in algo or "w_size" not in algo or "n_window" not in algo:
            exit(Error.MISSING_PARAM.value)
            
    elif algo["name"] == "normal-equation":
        if "error" not in algo or "mode" not in algo or "degree" not in algo:
            exit(Error.MISSING_PARAM.value)
        elif algo["mode"] != "individual" and algo["mode"] != "accumulate":
            exit(Error.INVALID_PARAM.value)
    
    elif algo["name"] == "optimal-pla":
        if "error" not in algo:
            exit(Error.MISSING_PARAM.value)

    elif algo["name"] == "semi-optimal-pla":
        if "error" not in algo:
            exit(Error.MISSING_PARAM.value)
    
    elif algo["name"] == "semi-mixed-pla":
        if "error" not in algo:
            exit(Error.MISSING_PARAM.value)
        
    elif algo["name"] == "swing-filter":
        if "error" not in algo:
            exit(Error.MISSING_PARAM.value)
        
    elif algo["name"] == "slide-filter":
        if "error" not in algo:
            exit(Error.MISSING_PARAM.value)

    elif algo["name"] == "mix-piece":
        if "error" not in algo or "n_segment" not in algo: 
            exit(Error.MISSING_PARAM.value)
    
    elif algo["name"] == "smart-grid-compression":
        if "error" not in algo or "degree" not in algo:
            exit(Error.MISSING_PARAM.value)

    elif algo["name"] == "unbounded":
        if "error" not in algo:
            exit(Error.MISSING_PARAM.value)

    elif algo["name"] == "bounded":
        if "error" not in algo or "degree" not in algo:
            exit(Error.MISSING_PARAM.value)

    else:
        exit(Error.INVALID_ALGO.value)
    
    
def parse(conf):
    print(conf["data"])
    print(conf["compress"])
    print(conf["decompress"])
    print(conf["interval"])
    
    algo = conf["algorithm"]
    if algo["name"] == "pmc":
        print("{} {} {}".format(algo["name"], algo["error"], algo["mode"]))
    elif algo["name"] == "hybrid-pca":
        print("{} {} {} {}".format(algo["name"], algo["error"], algo["w_size"], algo["n_window"]))
    elif algo["name"] == "normal-equation":
        print("{} {} {} {}".format(algo["name"], algo["error"], algo["mode"], algo["degree"]))
    elif algo["name"] == "optimal-pla":
        print("{} {}".format(algo["name"], algo["error"]))
    elif algo["name"] == "semi-optimal-pla":
        print("{} {}".format(algo["name"], algo["error"]))
    elif algo["name"] == "semi-mixed-pla":
        print("{} {}".format(algo["name"], algo["error"]))
    elif algo["name"] == "swing-filter":
        print("{} {}".format(algo["name"], algo["error"]))
    elif algo["name"] == "slide-filter":
        print("{} {}".format(algo["name"], algo["error"]))
    elif algo["name"] == "mix-piece":
        print("{} {} {}".format(algo["name"], algo["error"], algo["n_segment"]))
    elif algo["name"] == "smart-grid-compression":
        print("{} {} {}".format(algo["name"], algo["error"], algo["degree"]))
    elif algo["name"] == "unbounded":
        print("{} {}".format(algo["name"], algo["error"]))
    elif algo["name"] == "bounded":
        print("{} {} {}".format(algo["name"], algo["error"], algo["degree"]))
        

    
if __name__ == "__main__":
    with open(sys.argv[1], 'r') as file:
        conf = json.load(file)
    
    validate(conf['algorithm'])
    parse(conf)
    
    exit(0)
        