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

data = get_data("../data3D.nosync" + "/test_final.txt", ["c", "L", "D", "b", "E", "EE"])

f = []
cLLD = []

maxb = np.max(data["b"])

for i in range(len(data["c"])):
    if data["b"][i] == maxb:
        f.append(data["E"][i])
        cLLD.append(data["c"][i]*data["L"][i]**2)


cLLD = np.array(cLLD)
f = np.array(f)

indices = np.argsort(cLLD)
cLLD = cLLD[indices]
f = f[indices]

#plt.plot(1/cLL, cLL/np.pi)
#plt.plot(cLLD, cLLD*3*np.pi/16) #+ np.log(cLLD)
plt.plot(cLLD, f, 'o-')
plt.xlabel(r"$cL^2D$")
plt.ylabel(r"$f$")
plt.show()
