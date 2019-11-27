import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

infile = open("../data3D.nosync" + "/test_occupations.txt", "r")
data = get_data("../data3D.nosync" + "/test_final.txt", ["c", "L", "b", "E", "EE"])
maxb = np.max(data["b"])
i = 0
for line in infile:
    if data["b"][i] == maxb:
        plt.figure()
        f = np.array(line.split(), float)
        phi = np.linspace(0, 2*np.pi, len(f), endpoint = False)
        plt.plot(phi, f, label = "cLL = %.2f" % (data["c"][i]*data["L"][i]*data["L"][i]))
        plt.legend()
    i += 1
plt.show()
