#!/bin/python
import poinc

"""Script that provides a list of poinc_t=******.bin files into the nemato submit script."""
#Inputs
directory = os.getcwd()

def make_plots(directory):
    poinc.list_bin(directory)
    
#EXECUTION
make_plots(directory)
