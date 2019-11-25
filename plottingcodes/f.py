import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

infile = open("../data.nosync" + "/test_occupations.txt", "r")
N = 0
for line in infile:
    plt.figure()
    f = np.array(line.split(), float)
    phi = np.linspace(0, 2*np.pi, len(f), endpoint = False)
    plt.plot(phi, f)
    plt.show()
