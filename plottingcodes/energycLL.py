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

data = get_data("../data3D.nosync" + "/test_final.txt", ["c", "L", "b", "E", "EE"])

f = []
cLL = []

maxb = np.max(data["b"])

for i in range(len(data["c"])):
    if data["b"][i] == maxb:
        f.append(data["E"][i])
        cLL.append(data["c"][i]*data["L"][i]**2)


cLL = np.array(cLL)

#plt.plot(1/cLL, cLL/np.pi)
plt.plot(cLL, -cLL*np.pi/32)
plt.plot(cLL, f, 'o')
plt.xlabel(r"$cL^2$")
plt.ylabel(r"$f$")
plt.show()
