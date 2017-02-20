import numpy as np


def next2pow(x):
    return int(np.ceil(np.log(float(x))/np.log(2.0)))