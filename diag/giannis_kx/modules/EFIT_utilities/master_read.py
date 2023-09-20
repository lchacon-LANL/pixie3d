''''Master Reader for EFIT files'''

import numpy as np
import math
import numpy.ma as ma
import os

#filename = '/Users/Giannis/freegs/example.geqdsk'

class ds():
    def __init__(self,nw, nh, rdim, zdim, rleft, zmid, rmaxis, zmaxis, simag, sibry, rcentr, bcentr, 
        current, fpol, pres, ffprim, pprime, psirz, qpsi, nbbbs, limitr, rbbbs, zbbbs, rlim, zlim, kvtor, 
        rvtor, nmass, pressw, pwprim, dmion, rhovn):
        self.nw = nw
        self.nh = nh
        self.rdim = rdim
        self.zdim = zdim
        self.rleft = rleft
        self.zmid = zmid
        self.rmaxis = rmaxis
        self.zmaxis = zmaxis
        self.simag = simag
        self.sibry = sibry
        self.rcentr = rcentr
        self.bcentr = bcentr
        self.current = current
        self.fpol = fpol
        self.pres = pres
        self.ffprim = ffprim
        self.pprime = pprime
        self.psirz = psirz
        self.qpsi = qpsi
        self.nbbbs = nbbbs
        self.limitr = limitr
        self.rbbbs = rbbbs
        self.zbbbs = zbbbs
        self.rlim = rlim
        self.zlim = zlim
        self.kvtor = kvtor
        self.rvtor = rvtor
        self.nmass = nmass
        self.pressw = pressw
        self.pwprim = pwprim
        self.dmion = dmion
        self.rhovn = rhovn

def num_chunk(s):
    n = 16
    for ind in range(0,len(s),n):
        yield s[ind:ind+n]
        
def l2f(l):
    nl = []
    for elem in l:
        nl.append(float(elem))
    return nl
        
def fpol(nw_l):
    s = ''
    for line in all_lines[5:5+nw_l]:
        s = s + line.rstrip() 
    fpol = list(num_chunk(s))[0:DS.nw]
    return l2f(fpol)

def pres(nw_l):
    s = ''
    for line in all_lines[5+nw_l:5+2*nw_l]:
        s = s + line.rstrip() 
    pres = list(num_chunk(s))[0:DS.nw]
    return l2f(pres)

def ffprim(nw_l):
    s = ''
    for line in all_lines[5+2*nw_l:5+3*nw_l]:
        s = s + line.rstrip() 
    ffprim = list(num_chunk(s))
    return l2f(ffprim)[0:DS.nw]

def pprime(nw_l):
    s = ''
    for line in all_lines[5+3*nw_l:5+4*nw_l]:
        s = s + line.rstrip() 
    pprime = list(num_chunk(s))[0:DS.nw]
    return l2f(pprime)

def psirz(arr_l):
    s = ''
    for line in all_lines[5+4*nw_l:5+4*nw_l+arr_l]: # Jason: fixed error where nw, arr_l were rounded down by 1
        s = s + line.rstrip() 
    psi = list(num_chunk(s))[0:DS.nw*DS.nh]
    psirz = []
    for elem in psi:
        psirz.append(float(elem))
    psirz = np.asarray(psirz)
    psirz = np.reshape(psirz,(int(DS.nh),int(DS.nw)))
    return psirz

def qpsi(nw_l):
    s = ''
    for line in all_lines[5+4*nw_l+arr_l:5+5*nw_l+arr_l]:
        s = s + line.rstrip() 
    qpsi = list(num_chunk(s))[0:DS.nw]
    return l2f(qpsi)

def boundary_and_limiter():
    line = all_lines[5+5*nw_l+arr_l].rstrip()
    spl = line.split() # Jason: default split is any white space, do not specify a specific white space character
    non_empty = [s for s in spl if s]
    nbbbs = non_empty[0]    
    limitr = non_empty[1]
    return int(nbbbs), int(limitr)

def rbbbs_and_zbbbs(bbbs_l):
    s = ''
    for line in all_lines[6+5*nw_l+arr_l:6+5*nw_l+arr_l+bbbs_l]:
        s = s + line.rstrip() 
    rbbbs = list(num_chunk(s))[0:2*DS.nbbbs:2]
    zbbbs = list(num_chunk(s))[1:2*DS.nbbbs:2]
    return l2f(rbbbs), l2f(zbbbs)

def rlim_and_zlim(lim_l):
    s = ''
    for line in all_lines[6+5*nw_l+arr_l+bbbs_l:6+5*nw_l+arr_l+bbbs_l+lim_l]:
        s = s + line.rstrip() 
    rlim = list(num_chunk(s))[0:2*DS.limitr:2]
    zlim = list(num_chunk(s))[1:2*DS.limitr:2]
    return l2f(rlim), l2f(zlim)

def kvtor():
    kvtor = 0
    try:
        line = all_lines[6+5*nw_l+arr_l+bbbs_l+lim_l].rstrip()
        spl = line.split()
        non_empty = [s for s in spl if s]
        kvtor = non_empty[0]
    except IndexError:
        print('End of File')
    except:
        print('Something went wrong')
    return float(kvtor)

def rvtor_nmass():
    rvtor = 0
    nmass = 0
    try:
        line = all_lines[7+5*nw_l+arr_l+bbbs_l+lim_l].rstrip()
        rvtor = list(num_chunk(line))[0]
        nmass = list(num_chunk(line))[1]
    except IndexError:
        print('End of File')
    except:
        print('Something went wrong')
    return float(rvtor), float(nmass)

def rest_of_file():
    s = ''
    for line in all_lines[7+5*nw_l+arr_l+bbbs_l+lim_l:]:
        s = s + line.rstrip()
    rest = list(num_chunk(s))[2:]
    return rest

def pressw():
    pressw = rest[0:DS.nw]
    return l2f(pressw)

def pwprim():
    pwprim = rest[int(DS.nw):2*DS.nw]
    return l2f(pwprim)

def dmion(select):
    if select == 1: #kvtor > 0
        if len(rest) >= 3*DS.nw:
            dmion = rest[2*DS.nw:3*DS.nw]
        else:
            dmion = rest[2*DS.nw:]
    if select == 2: #kvtor = 0
        if len(rest) >= DS.nw:
            dmion = rest[0:DS.nw]
        else:
            dmion = rest[0:]
    return l2f(dmion)
        
def rhovn(select):
    if select == 1: #kvtor = 0, nmass = 0
        if len(rest) >= DS.nw:
            rhovn = rest[0:DS.nw]
        else:
            rhovn = rest[0:]
    if select == 2: #kvtor = 0, nmass > 0
        if len(rest) >= 2*DS.nw:
            rhovn = rest[DS.nw:2*DS.nw]
        else:
            rhovn = rest[DS.nw:]
    if select == 3: #kvtor > 0, nmass = 0
        if len(rest) >= 3*DS.nw:
            rhovn = rest[2*DS.nw:3*DS.nw]
        else:
            rhovn = rest[2*DS.nw:]
    if select == 4: #kvtor > 0, nmass > 0
        if len(rest) >= 4*DS.nw:
            rhovn = rest[3*DS.nw:4*DS.nw]
        else:
            rhovn = rest[3*DS.nw:]
    return l2f(rhovn)


def read_geqdsk(geqdsk):
    
    global nw, nh, rdim, zdim, rleft, zmid, rmaxis, zmaxis, simag, sibry, rcentr, bcentr 
    global current, fpol, pres, ffprim, pprime, psirz, qpsi, nbbbs, limitr, rbbbs, zbbbs, rlim, zlim, kvtor 
    global rvtor, nmass, pressw, pwprim, dmion, rhovn
    global all_lines, DS, nw_l, arr_l, bbbs_l, lim_l, rest
    
    DS = ds(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    
    fid = open(geqdsk,'r')

    all_lines = fid.readlines()
    
    first_line = all_lines[0].rstrip()
    second_line = all_lines[1].rstrip()
    third_line = all_lines[2].rstrip()
    fourth_line = all_lines[3].rstrip()
    fifth_line = all_lines[4].rstrip()

    flne = [s for s in first_line.split() if s]
    info = flne[0:-2]
    print('Case Info:', info)
    DS.nw = int(flne[-2])
    print('nw:', DS.nw)
    DS.nh = int(flne[-1])
    print('nh:', DS.nh)
    
    line2 = list(num_chunk(second_line))
    DS.rdim = float(line2[0])
    print('rdim:', DS.rdim)
    DS.zdim = float(line2[1])
    print('zdim:', DS.zdim)
    DS.rcentr = float(line2[2])
    DS.rleft = float(line2[3])
    print('rleft:', DS.rleft)
    DS.zmid = float(line2[4])
    print('zmid:', DS.zmid)
    
    line3 = list(num_chunk(third_line))
    DS.rmaxis = float(line3[0])
    print('rmaxis:', DS.rmaxis)
    DS.zmaxis = float(line3[1])
    print('zmaxis:', DS.zmaxis)
    DS.simag = float(line3[2])
    print('simag:', DS.simag)
    DS.sibry = float(line3[3])
    print('sibry:', DS.sibry)
    DS.bcentr = float(line3[4])
    print('rcentr:', DS.rcentr)
    print('bcentr:', DS.bcentr)
    
    line4 = list(num_chunk(fourth_line))
    DS.current = float(line4[0])
    print('current:', DS.current)
    #DS.simag = float(line4[1])
    #DS.rmaxis = float(line4[3])
    
    line5 = list(num_chunk(fifth_line))
    #DS.zmaxis = float(line5[0])
    #DS.sibry = float(line5[2])
    
    '''Calculation of lines with nw characters'''
    nw_l = int(math.ceil(float(DS.nw)/5)) # Jason: changed to force integer type
    
    DS.fpol = fpol(nw_l)
    DS.pres = pres(nw_l)
    DS.ffprim = ffprim(nw_l)
    DS.pprime = pprime(nw_l)
    print('fpol:', len(DS.fpol))
    print('pres:', len(DS.pres))
    print('ffprim:', len(DS.ffprim))
    print('pprime:', len(DS.pprime))
    
    """Calculation of lines for psi matrix"""
    arr_l = int(math.ceil(float(DS.nw*DS.nh)/5)) # Jason: changed to force integer type
    
    DS.psirz = psirz(arr_l)
    print('psirz:', DS.psirz.shape)
    DS.qpsi = qpsi(nw_l)
    print('qpsi:', len(DS.qpsi))
    
    DS.nbbbs, DS.limitr = boundary_and_limiter()
    print('nbbbs:', DS.nbbbs)
    print('limitr:', DS.limitr)
    
    '''Calculation of lines for bbbs lists'''
    bbbs_l = int(math.ceil(float(2*DS.nbbbs)/5)) # Jason: changed to force integer type
    
    '''Calculation for lim lists'''
    lim_l = int(math.ceil(float(2*DS.limitr)/5)) # Jason: changed to force intger type
    
    DS.rbbbs, DS.zbbbs = rbbbs_and_zbbbs(bbbs_l)
    print('rbbbs:', len(DS.rbbbs))
    print('zbbbs:', len(DS.zbbbs))
    DS.rlim, DS.zlim = rlim_and_zlim(lim_l)
    print('rlim:', len(DS.rlim))
    print('zlim:', len(DS.zlim))
    DS.kvtor = kvtor()
    print('kvtor:', DS.kvtor)
    DS.rvtor, DS.nmass = rvtor_nmass()
    print('rvtor:', DS.rvtor)
    print('nmass:', DS.nmass)
    
    rest = rest_of_file()
    
    if float(DS.kvtor) > 0:
        rest = rest_of_file()
        DS.pressw = pressw()
        print('pressw:', len(DS.pressw))
        DS.pwprim = pwprim()
        print('pwprim:', len(DS.pwprim))
    else:
        pass
    
    if float(DS.kvtor) > 0 and float(DS.nmass) > 0:
        DS.dmion = dmion(1)
        DS.rhovn = rhovn(4)
    if float(DS.kvtor) == 0 and float(DS.nmass) > 0:
        DS.dmion = dmion(2)
        DS.rhovn = rhovn(2)
    if float(DS.kvtor) == 0 and float(DS.nmass) == 0:
        DS.rhovn = rhovn(1)
    if float(DS.kvtor) > 0 and float(DS.nmass) == 0:
        DS.rhovn = rhovn(1)
    try:    
        print('dmion:', len(DS.dmion))
        print('rhovn:', len(DS.rhovn))
    except:
        pass
    
    fid.close()

def struct_hor_ax_det():
    
    global rr, zz
    global grid_mag_center_r_index, grid_mag_center_z_index, lim_left_int, lim_right_int
    global sep_left_int, sep_right_int, r_lim_left_intersection, r_lim_right_intersection
    global r_lim_index_left_intersection, r_lim_index_right_intersection
    global r_sep_left_intersection, r_sep_right_intersection
    global r_sep_index_left_intersection, r_sep_index_right_intersection
    
    rr = np.linspace(DS.rleft, DS.rleft + DS.rdim, DS.nw)
    zz = np.linspace(-DS.zdim/2 + DS.zmid, DS.zdim/2 +DS.zmid, DS.nh)
    
    '''structural horizontal axis detector'''
    if DS.limitr > 0:
        ## Index Ranges in which the intersection is searched
        right_int_range = [x for x in range(0,int(math.floor(DS.limitr/2.))+1)]
        left_int_range = [x for x in range(10,DS.limitr)]

        ## Index Ranges in which the intersection with the separatrix is searched
        right_int_range_sep = [x for x in range(0,10)]
        left_int_range_sep = [x for x in range(10, len(DS.rbbbs)-10)]

        ## grid indexes of magnetic center
        a1 = list(np.abs(rr - DS.rmaxis))
        grid_mag_center_r_index = a1.index(min(a1))
        a2 = list(np.abs(zz - DS.zmaxis))
        grid_mag_center_z_index = a2.index(min(a2))

        ## r indexes of intersection of z = zmaxis with device
        ## distance between the z device boundary coords and zmaxis
        dist_vect_from_z_axis  = np.abs(DS.zmaxis - np.array(DS.zlim))
        ## r_left and r_right index of device boundary on mag axis line
        a3 = [dist_vect_from_z_axis[x] for x in left_int_range]
        lim_left_int = a3.index(min(a3))
        a4 = [dist_vect_from_z_axis[x] for x in right_int_range]
        lim_right_int = a4.index(min(a4))
        lim_left_int = lim_left_int + left_int_range[0]
        lim_right_int = lim_right_int + right_int_range[0]
    
        ## r indexes of of intersection of z =zmaxis with separatrix
        ## distance between the z separatrix coords and zmaxis
        dist_vect_from_z_axis2 = np.abs(DS.zmaxis - np.array(DS.zbbbs))
        ## r_left and r_right index of separatrix boundary on mag axis line
        a5 = [dist_vect_from_z_axis2[x] for x in left_int_range_sep]
        sep_left_int = a5.index(min(a5))
        a6 = [dist_vect_from_z_axis2[x] for x in right_int_range_sep]
        sep_right_int = a6.index(min(a6))
        sep_left_int = sep_left_int + left_int_range_sep[0]
        sep_right_int = sep_right_int + right_int_range_sep[0]
    
        ## external_zero_index = external_zero_index_a(end)
    
        ## intersection grid indexes
    
        r_lim_left_intersection = DS.rlim[lim_left_int] 
        a7 = list(np.abs(np.array(rr)-r_lim_left_intersection))
        r_lim_index_left_intersection = a7.index(min(a7))
        r_lim_right_intersection = DS.rlim[lim_right_int] 
        a8 = list(np.abs(np.array(rr)-r_lim_right_intersection))
        r_lim_index_right_intersection = a8.index(min(a8))
    
        r_sep_left_intersection = DS.rbbbs[sep_left_int] 
        a9 = list(np.abs(np.array(rr)-r_sep_left_intersection))
        r_sep_index_left_intersection = a9.index(min(a9))
        r_sep_right_intersection = DS.rbbbs[sep_right_int] 
        a10 = list(np.abs(np.array(rr)-r_sep_right_intersection))
        r_sep_index_right_intersection = a10.index(min(a10))
    else:
            pass

