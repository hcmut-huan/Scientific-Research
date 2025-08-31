import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

ERROR = 4.355

def load(file, index):
    data = []

    with open(file, 'r') as file:
        csvFile = csv.reader(file)
        for line in csvFile:
            data.append([int(line[0]), float(line[index])])
            
        data = sorted(data, key=lambda a_entry: a_entry[0]) 
    
    return np.array(data)[:, 0], np.array(data)[:, 1]


def plot(s, e, label, time, data):
    plt.figure(figsize=(30, 8))

    for i in range(len(label)):
        plt.plot(time[i][s:e], data[i][s:e], label=label[i])

    plt.legend()
    plt.show()

if __name__ == "__main__":
    label, x, y = [], [], []
    for i in range(1, len(sys.argv) - 1, 2):
        time, value = load(sys.argv[i+1], 1)
        
        label.append(sys.argv[i])
        x.append(time)
        y.append(value)

    SIZE = len(x[0])
    INTERVAL = 10000

    plot(0, SIZE, label, x, y)
    start = 0
    while start < SIZE:
        plot(start, start+INTERVAL, label, x, y)
        start += INTERVAL