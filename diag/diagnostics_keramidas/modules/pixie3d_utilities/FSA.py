import numpy as np
from scipy.integrate import romb

"""Module that contains functions that do flux surface averaging."""

def Dromb(Arr, r_ind):
    udim = Arr.shape[1]
    phidim = Arr.shape[2]
    I = np.zeros(udim) # phi integral list over all us
    for i in range(udim):
        I[i] = romb(Arr[r_ind,i,:],dx=2*np.pi/(phidim-1))
    Int = romb(I,dx=2*np.pi/(udim-1)) # u integral
    return Int

def FSA_denominator(jac):
    denominator = np.zeros(jac.shape[0])
    for i in range(jac.shape[0]):
        denominator[i] = Dromb(jac,i)
    return denominator

def FSA_numerator(A,jac):
    integrand = jac*A
    numerator = np.zeros(integrand.shape[0])
    for i in range(integrand.shape[0]):
        numerator[i] = Dromb(integrand,i)
    return numerator

def FSA_total(A,jac):
    numerator = FSA_numerator(A,jac)
    denominator = FSA_denominator(jac)
    fsa = numerator/denominator
    return fsa
