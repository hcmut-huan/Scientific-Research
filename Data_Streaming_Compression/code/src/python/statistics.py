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


def load_monitor(compress, decompress):
    c_max_vsz, c_max_rss = float("-inf"), float("-inf")
    d_max_vsz, d_max_rss = float("-inf"), float("-inf")
    
    with open(compress, 'r') as file:
        csvFile = csv.reader(file)
        next(csvFile, None)
        for line in csvFile:
            c_vsz = int(line[2])
            c_rss = int(line[3])

            c_max_vsz = c_max_vsz if c_max_vsz > c_vsz else c_vsz
            c_max_rss = c_max_rss if c_max_rss > c_rss else c_rss
                    
    with open(decompress, 'r') as file:
        csvFile = csv.reader(file)
        next(csvFile, None)
        for line in csvFile:
            d_vsz = int(line[2])
            d_rss = int(line[3])

            d_max_vsz = d_max_vsz if d_max_vsz > d_vsz else d_vsz
            d_max_rss = d_max_rss if d_max_rss > d_rss else d_rss
        
    return c_max_vsz, c_max_rss, d_max_vsz, d_max_rss


def rmse(origin_data, approx_data):
    return np.sqrt(np.square(origin_data - approx_data).mean())


def mse(origin_data, approx_data):
    return np.square(origin_data - approx_data).mean()


def mae(origin_data, approx_data):
    return np.abs(origin_data - approx_data).mean()


def maxdiff(origin_data, approx_data):
    index = np.argmax(np.abs(origin_data - approx_data))
    return np.max(np.abs(origin_data - approx_data)) , index, origin_data[index], approx_data[index]


def mindiff(origin_data, approx_data):
    return np.min(np.abs(origin_data - approx_data))
   
   
def compressratio(origin, approx):
    return os.path.getsize(origin) / os.path.getsize(approx)
    

if __name__ == "__main__":
    DATA = sys.argv[1]
    DECOMPRESS = sys.argv[2]
    COMPRESS = sys.argv[3]
    C_MONITOR = COMPRESS + ".mon" 
    D_MONITOR = DECOMPRESS + ".mon"
    
    origin_data, approx_data = load(DATA, DECOMPRESS)
    c_max_vsz, c_max_rss, d_max_vsz, d_max_rss = load_monitor(C_MONITOR, D_MONITOR)
    count_min = min(origin_data.shape[0], approx_data.shape[0])
    
    origin_data = origin_data[:count_min]
    approx_data = approx_data[:count_min]
    
    
    # print(compressratio(DATA, COMPRESS), rmse(origin_data, approx_data), maxdiff(origin_data, approx_data)(0), c_max_vsz/1024/1024, c_max_rss/1024/1024, d_max_vsz/1024/1024,)
    
    print("Compress Ratio:", compressratio(DATA, COMPRESS))
    print("MSE:", mse(origin_data, approx_data))
    print("RMSE:", rmse(origin_data, approx_data))
    print("MAE:", mae(origin_data, approx_data))
    print("Max_E:", maxdiff(origin_data, approx_data))
    print("Min_E:", mindiff(origin_data, approx_data))
    print("c_max_vsz:", c_max_vsz/1024/1024)
    print("c_max_rss:", c_max_rss/1024/1024)
    print("d_max_vsz:", d_max_vsz/1024/1024)
    print("d_max_rss:", d_max_rss/1024/1024)