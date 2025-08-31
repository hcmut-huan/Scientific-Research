import sys
import numpy as np

file = sys.argv[1]
data = np.genfromtxt(file, delimiter=',')

value_range = np.max(data[:,1]) - np.min(data[:,1])


print("1% -> epsilon =", value_range/100)
print("2.5% -> epsilon =", value_range*2.5/100)
print("5% -> epsilon =", value_range*5/100)
print("7.5% -> epsilon =", value_range*7.5/100)
print("10% -> epsilon =", value_range*10/100)
