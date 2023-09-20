#!/bin/python
import poinc

"""Script that makes individual Poincare plots. Gives you the option to either make a plot with horizontal lines on different 
   flux surfaces or make a plot with the X-axis labeled by the initial poloidal flux. Works when the ssrz_poinc_t=****.in 
   files have already been converted to text files."""
#Inputs
directory = os.getcwd()
eqdsk_file = "/users/giannis_kx/eqdsks/eqdsk_9MA_SS.Gpolevoa"
position = 0 # choice of which of the sorted text files to process

def make_intersection_plots(directory,eqdsk_file,position):
    A = poinc.srtd_fl(directory)
    L=[A[position]]
    poinc.sfig(L,option=3,directory,eqdsk_file)

def make_plot_w_flux_label(directory,eqdsk_file,position):
    SA = poinc.srtd_fl(directory)
    L = SA[position]
    poinc.sfig(L,option=2,directory,eqdsk_file,position)

    
#EXECUTION
make_intersection_plots(directory,eqdsk_file,position)
make_plot_w_flux_label(directory,eqdsk_file,position)