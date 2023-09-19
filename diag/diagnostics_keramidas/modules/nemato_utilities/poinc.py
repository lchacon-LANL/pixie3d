#!/bin/python
import os
import sys
import glob
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import pylab as py
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import animation, rc
import matplotlib.animation as animation
import PIL
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 
import imageio
import master_read as m
import efit_plots as e

Image.MAX_IMAGE_PIXELS = 10000000000

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.labelweight'] = 'heavy'
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['text.usetex']=True

"""Collection of functions for processing poinc_t=***.bin and ssrz_poinc_t=****.in files. 

   Plotting customization is done at the top and in function sfig(). 
   In function resize_im() we specify the resolution of downsized images.
   
   CAUTION: Below functions are very sensitive to the file names and the format of the textfiles produced by xdraw2text.
   Common sourse of problems is when text files produced by xdraw2text change column label from 'R' to 'r'. """


"""Function that finds all the bin files."""
def find_bin(directory):
    bin_list = [f for f in os.listdir(directory) if (f.startswith("poinc_t=") and f.endswith(".bin"))]
    return bin_list

"""Function that finds which bin files have already been converted in .in files."""
def find_existing_infile(directory):
    bin_list = find_bin(directory)
    in_list = find_in(directory)
    new_bin_list = find_bin(directory)
    for i,filebin in enumerate(bin_list):
        tstamp = filebin[8:-4]
        for filein in in_list:
            if "drawssrz_poinc_t="+tstamp+".in" == str(filein):
                new_bin_list.remove(filebin)
    return new_bin_list

"""Lists all bin files that need to be converted to in files."""
def list_bin(directory):
	bin_list = find_existing_infile(directory)
	file_string = ""
	for counter,file in enumerate(bin_list):
		if counter%1==0:
			file_string = file_string + str(file) +" " + "*"
		else:
			file_string = file_string + str(file) + " "
	print(file_string)

"""Number of bin files that need to be converted to in files."""
def bin_num(directory):
    bin_list = find_existing_infile(directory)
    return len(bin_list)

def find_in(directory):
	in_list = []
	for file in os.listdir(directory):
		if file.startswith("drawssrz"):
			in_list.append(file)
	return in_list

def srtd(list_bin):
	fileDict = {}
	for file in list_bin:
		key = file.split('.')[0][8:]
		fileDict[int(key)] = {file}
	file_list = []
	for key in sorted(fileDict.keys()):
		file_list.append(fileDict[key])
	new_list = []
	for file in file_list:
		new_list.append(str(file)[2:-2])
	return list(new_list)

"""Return all txt files."""
def textlist(directory):
    txt_list = []
    for file in os.listdir(directory):
        if file.endswith(".txt") and file.startswith("ssrz_poinc"):
            txt_list.append(file)
    return txt_list

"""Turns all .in files into text files. Ignores the ones that have already been done."""
def textmaker(directory):
    txt_list = textlist(directory)
    in_list = find_in(directory)
    new_in_list = find_in(directory)
    for i,file1 in enumerate(in_list):
        for file2 in txt_list:
            if file1[4:-3]+".txt" == str(file2):
                new_in_list.remove(file1)
    for file in new_in_list:
        os.system("xdraw2ascii "+file)

"""Function that returns the txt files sorted by time stamp."""
def srtd_fl(directory):
    txtFiles = []
    for file in glob.glob(directory + '/ssrz_poinc*.txt'):
        txtFiles.append(file[len(directory)+1:])
        
    fileDict = {}
    for file in txtFiles:
        key = file.split('=')[1].split('.')[0]
        fileDict[int(key)] = {file}
    file_list = []
    for key in sorted(fileDict.keys()):
        file_list.append(fileDict[key])
    return list(file_list)

"""Function that finds which txt files have not already been made into .png files."""
def find_png2cnv():
    txt_files = srtd_fl(directory)
    img_files = srtd_im(directory)
    txt2cnvrt = srtd_fl(directory)
    for i, tfile in enumerate(txt_files):
        tstamp_txt = str(tfile)[15:-6]
        for imfile in img_files:
            tstamp_im = str(imfile)[4:-6]
            if tstamp_txt == tstamp_im:
                txt2cnvrt.remove(tfile)
    return txt2cnvrt


"""Function that reads a Poincare text file."""
def poinc_reader(file,directory):    
    lstcnt = 0
    with open(directory+'/'+file,"r") as p_file:
        for line in p_file:
            if line[0] == 'R':
                lstcnt = lstcnt+1
            else:
                pass
        R = [[] for i in range(1,lstcnt+4)]
        Z = [[] for i in range(1,lstcnt+4)]
        j = 1
    
    p_file = open(directory+'/'+file,"r")
    for columns in (raw.strip().split() for raw in p_file):
        try:
            R[j].append(float(columns[0])*e.a)
            Z[j].append(float(columns[1])*e.a)
        except ValueError:
            j = j + 1
    p_file.close()
    return R,Z,lstcnt        

def poinc_reader_normed(file,directory):    
    lstcnt = 0
    with open(directory+'/'+file,"r") as p_file:
        for line in p_file:
            if line[0] == 'R':
                lstcnt = lstcnt+1
            else:
                pass
        R = [[] for i in range(1,lstcnt+4)]
        Z = [[] for i in range(1,lstcnt+4)]
        j = 1
    
    p_file = open(directory+'/'+file,"r")
    for columns in (raw.strip().split() for raw in p_file):
        try:
            R[j].append(float(columns[0]))
            Z[j].append(float(columns[1]))
        except ValueError:
            j = j + 1
    p_file.close()
    return R,Z,lstcnt        


"""Function that produces the scatter plots and saves them."""
def sfig(file_list,option,directory,eqdsk_file):
    
    eqdsk_info(eqdsk_file)
    Rs, Zs = shaping_plot()
    x_r = max(Rs)+0.1
    x_l = min(Rs)-0.1
    z_u = max(Zs)+0.1
    z_d = min(Zs)-0.1
    r_bnd = [x for x in m.DS.rlim]
    z_bnd = [y for y in m.DS.zlim]
    r_sep = [x for x in m.DS.rbbbs]
    z_sep = [y for y in m.DS.zbbbs]
   
    fig,ax = plt.subplots(figsize=(4,4), dpi=600)
    
    
    for file in file_list:
        R,Z,lstcnt = poinc_reader(str(list(file)[0]),directory)
        for i in range(lstcnt):
            poinc = ax.scatter(R[i][:],Z[i][:],s=0.1,linewidths=0.2)
            ax.set_xlabel('R(m)')
            ax.set_ylabel('Z(m)')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylim(z_d,z_u)
            ax.set_xlim(x_l,x_r)
            ax.set_aspect(aspect=1)
        ax.plot(Rs,Zs, color='blue', linewidth=2)
        ax.plot(r_bnd,z_bnd, color='black', linewidth=3)
        ax.plot(r_sep,z_sep, color='red', linewidth=1.5, linestyle='--')

        text = str(file).split("_")[2].split(".")[0]+"."+str(file).split("_")[2].split(".")[1]
        text_anot = str(file).split("_")[2].split(".")[0] + r"$\;\tau_A$"
        plt.annotate(text_anot, xy=(0.6, 0.94), xycoords='axes fraction',fontsize=14)
        if option == 1:
            plt.savefig(directory+'/'+'%s.png'%(text))
        if option == 2:
            newlab = [0.0,0.025,0.1,0.3,0.6,0.8,1.0]
            newpos = e.psin2R(newlab)
            newlabstr = [str(x) for x in newlab]
            ax.set_xticks(newpos)
            ax.set_xticklabels(newlabstr)
            plt.savefig(directory+'/'+'%s.png') 
        if option == 3:
            #ax.get_xaxis().set_visible(False)
            trans=ax.get_xaxis_transform()
            for idxi in e.idx32:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx32:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'3/2',fontsize=12,transform=trans)
            for idxi in e.idx43:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx43:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.06,'4/3',fontsize=12,transform=trans)
            for idxi in e.idx2:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx2:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'2/1',fontsize=12,transform=trans)
            for idxi in e.idx54:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx54:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'5/4',fontsize=12, transform=trans)
            for idxi in e.idx65: #start of small values
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]), ymax=0.75, color='green', linestyle='--',linewidth=2)
            for idxi in e.idx65:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'6/5',fontsize=12,transform=trans)
            for idxi in e.idx76:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]), ymax=0.75, color='green', linestyle='--',linewidth=2)
            for idxi in e.idx76:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'7/6',fontsize=12,transform=trans)
            for idxi in e.idx87:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]), ymax=0.75, color='green', linestyle='--',linewidth=2)
            for idxi in e.idx87:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'8/7',fontsize=12,transform=trans) #end
            for idxi in e.idx4:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx4:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'4/1',fontsize=12,transform=trans)
            for idxi in e.idx3:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx3:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'3/1',fontsize=12,transform=trans)
            for idxi in e.idx83:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx83:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'8/3',fontsize=12,transform=trans)
            for idxi in e.idx52:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx52:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.06,'5/2',fontsize=12,transform=trans)
            for idxi in e.idx73:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx73:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.09,'7/3',fontsize=12,transform=trans)
            for idxi in e.idx72:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx72:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.09,'7/2',fontsize=12,transform=trans)
            for idxi in e.idx5:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=2)
            for idxi in e.idx5:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'5/1',fontsize=12,transform=trans)
            plt.savefig(directory+'/'+'inter_''%s.png'%(text))
        plt.cla()


"""Function that returns the image files sorted."""
def srtd_im(directory):
    imFiles = []
    for file in glob.glob(directory + '/t=*.png'):
       # print(file)
        imFiles.append(file[len(directory)+1:])
    
    fileDict = {}
    for file in imFiles:
        key = file.split('=')[1].split('.')[0]
        fileDict[int(key)] = {file}
    file_list = []
    for key in sorted(fileDict.keys()):
        file_list.append(fileDict[key])
    return list(file_list)

"""Function that returns how many large image files have already been resized."""
def num_resized(directory):
    resized = []
    for file in glob.glob(directory + '/img*.png'):
        resized.append(file)
    return len(resized)

"""Function that resizes the images to lower resolution."""
def resize_im(file_list):
    n = num_resized(directory) # Number of already resized images
    for fn in file_list[n:]:
        fnf = str(list(fn)[0])
        im = Image.open(directory+'/' + fnf)
        im = im.resize((1080,1080), Image.ANTIALIAS)
        a = str(list(fn)[0])
        b = str(file_list.index(fn))
        c = b.rjust(4,'0')
        im.save("%s/img%s.png"%(directory,c))
      
"""Functions that read the equilibrium file and produce the boundary of the plot."""
def eqdsk_info(eqdsk_file):
    global a
    sys.stdout = open(os.devnull, 'w') # supressing print statements from master_read
    m.read_geqdsk(eqdsk_file)
    m.struct_hor_ax_det()
    e.gridsNunits()
    e.intersections()
    sys.stdout = sys.__stdout__ # reattach stdout

def shaping_plot():
    Ro = 6.219577546
    Zo = 0.5143555944
    alpha = 2.24
    kappa = 1.9
    delta = 0.6
    zeta = 0.06

    theta = np.linspace(0,2*np.pi,100)
    Rs = Ro + alpha*np.cos(theta + np.arcsin(delta*np.sin(theta)))
    Zs = Zo + kappa*alpha*np.sin(theta + zeta*np.sin(2*theta))

    return Rs, Zs

def norm_shape(directory):
    Rs, Zs = shaping_plot()
    plt.plot(Rs/e.a,Zs/e.a)
    plt.savefig(directory+"/"+"normplot.png")

