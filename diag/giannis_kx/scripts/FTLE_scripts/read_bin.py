from scipy.io import FortranFile
import numpy as np

filename = "ic_test.bin"

def read():
    A = np.fromfile(filename,dtype='d')
    return A

