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
#import imageio
import master_read as m
import efit_plots as e

Image.MAX_IMAGE_PIXELS = 10000000000

plt.rcParams['axes.labelsize'] = 120
plt.rcParams['axes.labelweight'] = 'heavy'
plt.rcParams['xtick.labelsize'] = 90
plt.rcParams['ytick.labelsize'] = 90
#plt.rcParams['text.usetex']=True

directory = os.getcwd()

if len(sys.argv) == 1:		# Check if user did not specify an EQDSK file for plotting boundaries 
    eqdsk_file = "../SN_fr_11_sh.geqdsk"		# Default EQDSK file

else:
    eqdsk_file = str(sys.argv[1])	# Use user-specified file, with *relative path* to current-working directory

eqdsk_read = os.path.isfile(directory+"/"+eqdsk_file)  # If no file is found, do not plot boundaries

if not eqdsk_read:
    for fname in os.listdir(directory):
        if fname.endswith('.geqdsk'):	# One final check if some other .geqdsk file is found
            eqdsk_file = fname
            eqdsk_read = os.path.isfile(directory+"/"+eqdsk_file)
            break

if not eqdsk_read:
    for fname in os.listdir(os.path.dirname(directory)):
        if fname.endswith('.geqdsk'):   # Also check in the parent directory
            eqdsk_file = '../'+fname
            eqdsk_read = os.path.isfile(directory+"/"+eqdsk_file) 
            break

"""Function that finds all the bin files."""
def find_bin():
    bin_list = [f for f in os.listdir(directory) if (f.startswith("poinc_t=") and f.endswith(".bin"))]
    return bin_list

"""Function that finds which bin files have already been converted in .in files."""
def find_existing_infile():
    bin_list = find_bin()
    in_list = find_in()
    new_bin_list = find_bin()
    for i,filebin in enumerate(bin_list):
        tstamp = filebin[8:-4]
        for filein in in_list:
            if "drawssrz_poinc_t="+tstamp+".in" == str(filein):
                new_bin_list.remove(filebin)
    return new_bin_list

"""Lists all bin files that need to be converted to in files."""
def list_bin():
	bin_list = find_existing_infile()
	file_string = ""
	for counter,file in enumerate(bin_list):
		if counter%1==0:
			file_string = file_string + str(file) +" " + "*"
		else:
			file_string = file_string + str(file) + " "
	print(file_string)

"""Number of bin files that need to be converted to in files."""
def bin_num():
    bin_list = find_existing_infile()
    return len(bin_list)

def find_in():
	in_list = []
	for file in os.listdir(directory):
		if file.startswith("drawssrz"):
			in_list.append(file)
	return in_list

"""DEPRECATED. Draws poincare plots with xdraw."""
def poinc_plots():
	draw_files = find_in()
	for file in draw_files:
		string = file[17:-3]+".ps"
		os.system("xdraw"+" "+file[4:-3])
		os.rename("xdrawa.ps", string)	

def srtd(list_bin):
	fileDict = {}
	for file in list_bin:
		key = file.split('.')[0][8:]
#		fileDict[int(key)] = {file} #  Do not use curly brackets or each element will be a python 'set'
		fileDict[int(key)] = file
	file_list = []
	for key in sorted(fileDict.keys()):
		file_list.append(fileDict[key])
#	new_list = []
#	for file in file_list:
#		new_list.append(str(file)[2:-2]) # No need to clip brackets from file name if treated properly above
#	return list(new_list)
	return file_list

"""Return all txt files."""
def textlist():
    txt_list = []
    for file in os.listdir(directory):
        if file.endswith(".txt") and file.startswith("ssrz_poinc"):
            txt_list.append(file)
    return txt_list

"""Turns all .in files into text files. Ignores the ones that have already been done."""
def textmaker():
    txt_list = textlist()
    in_list = find_in()
    new_in_list = find_in()
    for i,file1 in enumerate(in_list):
        for file2 in txt_list:
            if file1[4:-3]+".txt" == str(file2) or "ssrz_poinc_t="+file1[8:-3]+".txt" == str(file2):
                new_in_list.remove(file1)
    for file in new_in_list:
        os.system("mv "+file+" drawssrz"+file[17:]) # Renames xdraw file to shorter name
        file = str("drawssrz"+file[17:])
        os.system("xdraw2ascii "+file)

"""Function that returns the txt files sorted by time stamp."""
def srtd_fl():
    txtFiles = []
    for file in glob.glob(directory + '/ssrz_poinc*.txt'):
        txtFiles.append(file[len(directory)+1:])
        
    fileDict = {}
    for file in txtFiles:
        key = file.split('=')[1].split('.')[0]
#        fileDict[int(key)] = {file}
        fileDict[int(key)] = file
    file_list = []
    for key in sorted(fileDict.keys()):
        file_list.append(fileDict[key])
    return list(file_list)

"""Function that finds which txt files have not already been made into .png files."""
def find_png2cnv():
    txt_files = srtd_fl()
    img_files = srtd_im()
    txt2cnvrt = srtd_fl()
    for i, tfile in enumerate(txt_files):
#        tstamp_txt = str(tfile)[15:-6] # Since the curly brackets were removed from the sorting fns, char count is updated
        tstamp_txt = str(tfile)[13:-4]
        for imfile in img_files:
#            tstamp_im = str(imfile)[4:-6]
            tstamp_im = str(imfile)[2:-4]
            if tstamp_txt == tstamp_im:
                txt2cnvrt.remove(tfile)
    return txt2cnvrt


"""Function that reads a Poincare text file."""
def poinc_reader(file):    
    lstcnt = 0
    with open(directory+'/'+file,"r") as p_file:
        for line in p_file:
            if line[0] == 'R':
                lstcnt = lstcnt+1
            else:
                pass
#        R = [[] for i in range(1,lstcnt+4)] # Not sure why these coordinate lists were increased in dimension by 4
#        Z = [[] for i in range(1,lstcnt+4)]
        R = [[] for i in range(lstcnt)]
        Z = [[] for i in range(lstcnt)]
        j = -1 # Start at -1 so that the header at the top of the file pushes the first data point to address 0
    
    p_file = open(directory+'/'+file,"r")
    for columns in (raw.strip().split() for raw in p_file):
        try:
            R[j].append(float(columns[0])*e.a)
            Z[j].append(float(columns[1])*e.a)
        except ValueError:
            j = j + 1
    p_file.close()
    return R,Z,lstcnt

def poinc_reader_normed(file):    
    lstcnt = 0
    with open(directory+'/'+file,"r") as p_file:
        for line in p_file:
            if line[0] == 'R':
                lstcnt = lstcnt+1
            else:
                pass
        R = [[] for i in range(lstcnt)]
        Z = [[] for i in range(lstcnt)]
        j = -1
    
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
def sfig(file_list,option):

    if eqdsk_read: 
        eqdsk_info()
        r_bnd = [x for x in m.DS.rlim]
        z_bnd = [y for y in m.DS.zlim] 
        r_sep = [x for x in m.DS.rbbbs]
        z_sep = [y for y in m.DS.zbbbs]
    Rs, Zs = shaping_plot()
    x_r = max(Rs)+0.1
    x_l = min(Rs)-0.1
    z_u = max(Zs)+0.1
    z_d = min(Zs)-0.1
   
    fig,ax = plt.subplots(figsize=(150,150))
    
    
    for file in file_list:
#        R,Z,lstcnt = poinc_reader(str(list(file)[0])) # Fixed typo
        if eqdsk_read: 
            R,Z,lstcnt = poinc_reader(str(file))
        else:
            R,Z,lstcnt = poinc_reader_normed(str(file))
        for i in range(lstcnt):
            poinc = ax.scatter(R[i][:],Z[i][:],s=240)
            ax.set_xlabel('R', fontsize=300)
            ax.set_ylabel('Z', fontsize=300)
            ax.tick_params(axis='both', which='major', labelsize=200)
            ax.tick_params(axis='both', which='minor', labelsize=200)
            if eqdsk_read:
                ax.set_ylim(z_d,z_u)
                ax.set_xlim(x_l,x_r)
            else:						# Without an EQDSK file, use one plot's max/mins for all future plots
                if file_list.index(file) == 0:
                    x_r = np.max(sum(R,[])) + 0.1
                    x_l = np.min(sum(R,[])) - 0.1
                    z_u = np.max(sum(Z,[])) + 0.1
                    z_d = np.min(sum(Z,[])) - 0.1
                ax.set_ylim(z_d,z_u)
                ax.set_xlim(x_l,x_r)
            ax.set_aspect(aspect=1)
        if eqdsk_read: ax.plot(Rs,Zs, color='blue', linewidth=20)
        if eqdsk_read: ax.plot(r_bnd,z_bnd, color='black', linewidth=30)
        if eqdsk_read: ax.plot(r_sep,z_sep, color='red', linewidth=30, linestyle='--')

        text = str(file).split("_")[2].split(".")[0]+"."+str(file).split("_")[2].split(".")[1]
        plt.annotate(text+' '+r'$\tau_{A}$', xy=(0.6, 0.94), xycoords='axes fraction',fontsize=200)
        if option == 1 or not eqdsk_read:
            plt.savefig(directory+'/'+'%s.png'%(text))
        if option == 2 and eqdsk_read:
            newlab = [0.0,0.025,0.1,0.3,0.6,0.8,1.0]
            newpos = e.psin2R(newlab)
            newlabstr = [str(x) for x in newlab]
            ax.set_xticks(newpos)
            ax.set_xticklabels(newlabstr)
            plt.savefig(directory+'/'+'%s.png'%(text)) 
        if option == 3 and eqdsk_read:
            #ax.get_xaxis().set_visible(False)
            trans=ax.get_xaxis_transform()
            for idxi in e.idx32:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx32:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'3/2',fontsize=120,transform=trans)
            for idxi in e.idx43:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx43:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.06,'4/3',fontsize=120,transform=trans)
            for idxi in e.idx2:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx2:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'2/1',fontsize=120,transform=trans)
            for idxi in e.idx54:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx54:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'5/4',fontsize=120, transform=trans)
            for idxi in e.idx65: #start of small values
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]), ymax=0.75, color='green', linestyle='--',linewidth=20)
            for idxi in e.idx65:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'6/5',fontsize=120,transform=trans)
            for idxi in e.idx76:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]), ymax=0.75, color='green', linestyle='--',linewidth=20)
            for idxi in e.idx76:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'7/6',fontsize=120,transform=trans)
            for idxi in e.idx87:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]), ymax=0.75, color='green', linestyle='--',linewidth=20)
            for idxi in e.idx87:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'8/7',fontsize=120,transform=trans) #end
            for idxi in e.idx4:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx4:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'4/1',fontsize=120,transform=trans)
            for idxi in e.idx3:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx3:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'3/1',fontsize=120,transform=trans)
            for idxi in e.idx83:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx83:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.01,'8/3',fontsize=120,transform=trans)
            for idxi in e.idx52:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx52:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.06,'5/2',fontsize=120,transform=trans)
            for idxi in e.idx73:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx73:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.09,'7/3',fontsize=120,transform=trans)
            for idxi in e.idx72:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx72:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.09,'7/2',fontsize=120,transform=trans)
            for idxi in e.idx5:
                ax.axvline(x=e.psin2R(e.psi_norm[idxi]),ymax=0.75, color='black', linestyle='--',linewidth=20)
            for idxi in e.idx5:
                ax.text(e.psin2R(e.psi_norm[idxi]),-0.03,'5/1',fontsize=120,transform=trans)
            plt.savefig(directory+'/'+'inter_''%s.png'%(text))
        plt.cla()



        plt.cla()

"""Function that returns the image files sorted."""
def srtd_im():
    imFiles = []
    for file in glob.glob(directory + '/t=*.png'):
       # print(file)
        imFiles.append(file[len(directory)+1:])
    
    fileDict = {}
    for file in imFiles:
        key = file.split('=')[1].split('.')[0]
#        fileDict[int(key)] = {file}
        fileDict[int(key)] = file
    file_list = []
    for key in sorted(fileDict.keys()):
        file_list.append(fileDict[key])
    mylist = list(file_list)
    for file in mylist:
        if str(file) == 't=0.000.png': # Skip first poincare.bin file as it contains no perturbation information
            mylist.remove(file)
    return mylist

"""Function that returns how many large image files have already been resized."""
def num_resized():
    resized = []
    for file in glob.glob(directory + '/img*.png'):
        resized.append(file)
    return len(resized)

def resize_im(file_list):
    n = num_resized() # Number of already resized images
    for fn in file_list[n:]:
#        fnf = str(list(fn)[0]) # Fixed typo
        fnf = str(fn)
        im = Image.open(directory+'/' + fnf)
        im = im.resize((1080,1080), Image.ANTIALIAS)
        a = str(list(fn)[0])
        b = str(file_list.index(fn))
        c = b.rjust(4,'0')
#        if c != '0000':	# Skip first poincare.bin file as it contains no perturbation information
#            im.save("%s/img%s.png"%(directory,c))
#        else:	# Save initial poincare.bin separately for diagnostic purposes
#            im.save("%s/img_0.png"%(directory))
        im.save("%s/img%s.png"%(directory,c))

"""Function that creates an mp4 from a file list of images."""
def mp4_write(file_list):
    images = []
    for fn in file_list:
        fnf = str(fn)
        images.append(imageio.imread(directory + '/' + fnf))
    imageio.mimsave(directory2+'/'+'sawtooth.mp4', images, fps = 1)


"""Functions that read the equilibrium file and produce the boundary of the plot."""
def eqdsk_info():
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

def norm_shape():
    Rs, Zs = shaping_plot()
    plt.plot(Rs/e.a,Zs/e.a)
    plt.savefig(directory+"/"+"normplot.png")

#After this point you need to run the ffmpeg command:
def movie():
    os.system('ffmpeg -y -r 5 -f image2 -s 1080x1080 -i img%04d.png -vcodec mjpeg -qscale:v 5 -pix_fmt yuv420p poincare.mp4')

"""Functions that perform sets of instructions."""
def make_movie():
    textmaker()
    A = find_png2cnv()
    sfig(A,option=1)
    A = srtd_im()
    resize_im(A)
    movie()

def make_plots():
    list_bin()

def make_intersection_plots():
    A = srtd_fl()
    L=[A[0]]
    sfig(L,option=3)

def make_plot_w_flux_label():
    A = find_bin()
    SA = srtd(find_bin)
    L = SA[109]
    sfig(L,option=2)


#make_intersection_plots()
#make_plots()
make_movie()
