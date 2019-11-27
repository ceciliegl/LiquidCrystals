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

data = get_data("../data.nosync" + "/test_expect.txt", ["c", "L", "b", "E", "EE"])


plt.semilogx(range(0, 100*len(data["E"]), 100), data["E"])
plt.show()
