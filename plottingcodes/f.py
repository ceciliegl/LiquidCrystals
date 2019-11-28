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
data = get_data("../data3D.nosync" + "/test_final.txt", ["c", "L", "D", "b", "E", "EE"])
maxb = np.max(data["b"])
i = 0
j = 0
for line in infile:
    if data["b"][i] == maxb:
        if j%1 == 0:
            plt.figure()
            f = np.array(line.split(), float)
            x = np.linspace(0, np.pi, len(f), endpoint = False)
            plt.plot(x, f, label = "cLL = %.2f" % (data["c"][i]*data["L"][i]*data["L"][i]))
            plt.ylim(0, 2*np.max(f))
            plt.legend()
        j += 1
    i += 1
plt.show()
