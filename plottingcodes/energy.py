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

data = get_data("../data.nosync" + "/test_expect.txt", ["b", "E", "EE"])


plt.semilogx(range(len(data["E"])), data["E"])
plt.show()
