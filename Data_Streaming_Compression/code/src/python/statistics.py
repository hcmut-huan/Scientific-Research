import os
import csv
import sys
import numpy as np


def load(origin, approx):
    origin_data = []
    approx_data = []

    with open(origin, 'r') as file:
        csvFile = csv.reader(file)
        for line in csvFile:
            origin_data.append([int(line[0]), float(line[1])])
        
        origin_data = sorted(origin_data, key=lambda a_entry: a_entry[0]) 
                
    with open(approx, 'r') as file:
        csvFile = csv.reader(file)
        for line in csvFile:
            approx_data.append([int(line[0]), float(line[1])])
        
        approx_data = sorted(approx_data, key=lambda a_entry: a_entry[0]) 
    
    return np.array(origin_data)[:, 1], np.array(approx_data)[:, 1]


def load_monitor(file):
    max_vsz, max_rss = float("-inf"), float("-inf")
    
    with open(file, 'r') as file:
        csvFile = csv.reader(file)
        next(csvFile, None)
        for line in csvFile:
            c_vsz = int(line[2])
            c_rss = int(line[3])

            max_vsz = max_vsz if max_vsz > c_vsz else c_vsz
            max_rss = max_rss if max_rss > c_rss else c_rss
        
    return max_vsz, max_rss


def rmse(origin_data, approx_data):
    return np.sqrt(np.square(origin_data - approx_data).mean())


def mse(origin_data, approx_data):
    return np.square(origin_data - approx_data).mean()


def mae(origin_data, approx_data):
    return np.abs(origin_data - approx_data).mean()


def maxdiff(origin_data, approx_data):
    index = np.argmax(np.abs(origin_data - approx_data))
    return np.max(np.abs(origin_data - approx_data)) #, index, origin_data[index], approx_data[index]


def mindiff(origin_data, approx_data):
    return np.min(np.abs(origin_data - approx_data))
   
   
def compressratio(origin, approx):
    return os.path.getsize(origin) / os.path.getsize(approx)
    

if __name__ == "__main__":
    DATA = sys.argv[1]
    DECOMPRESS = sys.argv[2]
    COMPRESS = sys.argv[3]
    MONITOR = ".mon"

    origin_data, approx_data = load(DATA, DECOMPRESS)
    max_vsz, max_rss = load_monitor(MONITOR)
    count_min = min(origin_data.shape[0], approx_data.shape[0])
    
    origin_data = origin_data[:count_min]
    approx_data = approx_data[:count_min]
        
    print("Compress Ratio:", compressratio(DATA, COMPRESS))
    print("MSE:", mse(origin_data, approx_data))
    print("RMSE:", rmse(origin_data, approx_data))
    print("MAE:", mae(origin_data, approx_data))
    print("Max_E:", maxdiff(origin_data, approx_data))
    print("Min_E:", mindiff(origin_data, approx_data))
    print("max_vsz:", max_vsz/1024/1024)
    print("max_rss:", max_rss/1024/1024)