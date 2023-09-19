#!/bin/python
import os
import poinc

"""Script that creates Poincare plot movie from ssrz_poinc* nemato files. Images are resized before made into animation."""

#Inputs
directory = os.getcwd()
eqdsk_file = "/users/giannis_kx/eqdsks/eqdsk_9MA_SS.Gpolevoa"
title = "double_tearing_no_diffusion.mp4"

def movie(title):
    os.system('ffmpeg -r 5 -f image2 -s 1080x1080 -i img%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p {}'.format(title))


def make_movie(directory,eqdsk_file,title):
    poinc.textmaker(directory)
    A = poinc.find_png2cnv(directory)
    poinc.sfig(A,option=1,directory,eqdsk_file)
    A = poinc.srtd_im(directory)
    poinc.resize_im(A)
    movie(title)

# EXECUTION    
make_movie(directory,eqddsk_file,title)