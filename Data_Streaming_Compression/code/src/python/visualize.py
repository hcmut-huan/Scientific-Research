import os
import sys
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt


def readCSV(inFile):
    data = []
    with open(inFile, 'r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for row in reader:
            data.append([int(row[0]), float(row[1])])
    
    return np.array(data)

NUM = 10000
EXTENSION = "pdf"

csv_files = glob.glob("data/*/*.csv")
for file in csv_files:
    if "Synthesis" in file:
        continue
        
    name = os.path.splitext(os.path.basename(file))[0]
    data = readCSV(file)
    # Extract x and y
    x = data[:NUM, 0]  # Timestamps
    y = data[:NUM, 1]  # SOC

    # Plot
    plt.figure(figsize=(12, 8))
    plt.plot(x, y, color='blue', linewidth=0.1)
    plt.title(name, fontsize=25, pad=15)
    # plt.xlabel("Timestamp")      
    # plt.ylabel("Value")
    plt.savefig("{}.{}".format(name, EXTENSION), format=EXTENSION, bbox_inches="tight")
