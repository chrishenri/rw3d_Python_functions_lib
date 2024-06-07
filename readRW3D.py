# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd


#--------------read BTCs


#--------------read CBTCs
def read_cbtc(filename):
    fid = open(filename,'r')
    Lines = fid.readlines()
    
    #get zone info
    zoneline = [i for i in Lines if i.startswith('ZONE')]
    x = np.array([zoneline[i].split() for i in range(len(zoneline))])
    if x[0,3]=='Spe': #We're dealing with a sink term
        control = np.array(x[:,2]).astype(str)
        species = np.array(x[:,4]).astype(int)
    else:                       #We're dealing with a well ora control plane
        control = np.array(x[:,3]).astype(float)
        species = np.array(x[:,5]).astype(int)
    
    ncbtc = len(x)

    #get data
    df_cbtc = {}
    zoneind = [Lines.index(l) for l in Lines if l.startswith('ZONE')]
    zoneind.append(len(Lines)-1)
    for i in range(ncbtc):
        data = Lines[zoneind[i]+1:zoneind[i+1]-1]
        data_val = np.array([data[k].split() for k in range(len(data))])
        df_cbtc[i] = pd.DataFrame(data_val, columns = ['time','mass','partID','xp','yp','zp','x0','y0','z0'])
        df_cbtc[i] = df_cbtc[i].astype(float)
        
    class output:
        pass
    cbtc = output()
    cbtc.ncbtc = ncbtc
    cbtc.control = control
    cbtc.species = species
    cbtc.data = df_cbtc
    return cbtc


#--------------Particle evolution
def read_partevol(filename):
    fid = open(filename,'r')
    Lines = fid.readlines()
    
    #get zone info
    zoneline = [i for i in Lines if i.startswith('zone')]
    x = [zoneline[i].split() for i in range(len(zoneline))]
    time = np.array([float(x[i][2]) for i in range(len(zoneline))])
    
    #create plume object
    class output:
        pass
    plume = output()
    plume.times = time
    
    df_t = {}
    spe = {}
    
    #get data
    zoneind = [Lines.index(l) for l in Lines if l.startswith('zone')]
    zoneind.append(len(Lines)-1)
    for i in range(len(time)):
        data = Lines[zoneind[i]+1:zoneind[i+1]]
        data_val = np.array([data[i].split() for i in range(len(data))])
        spe[i] = np.unique(data_val[:,6])
        spe[i] = spe[i].astype(int)
        df_t[i] = pd.DataFrame(data_val, columns = ['xp','yp','zp','mp','rp','zone','specie','id'])
        df_t[i] = df_t[i].astype(float)
        
    plume.data = df_t
    plume.spe = spe
    return plume

#--------------Particle evolution
def read_partevol_noID(filename):
    fid = open(filename,'r')
    Lines = fid.readlines()
    
    #get zone info
    zoneline = [i for i in Lines if i.startswith('zone')]
    x = [zoneline[i].split() for i in range(len(zoneline))]
    time = np.array([float(x[i][2]) for i in range(len(zoneline))])
    
    #create plume object
    class output:
        pass
    plume = output()
    plume.times = time
    
    df_t = {}
    spe = {}
    
    #get data
    zoneind = [Lines.index(l) for l in Lines if l.startswith('zone')]
    zoneind.append(len(Lines)-1)
    for i in range(len(time)):
        data = Lines[zoneind[i]+1:zoneind[i+1]]
        data_val = np.array([data[i].split() for i in range(len(data))])
        spe[i] = np.unique(data_val[:,6])
        spe[i] = spe[i].astype(int)
        df_t[i] = pd.DataFrame(data_val, columns = ['xp','yp','zp','mp','rp','zone','specie'])
        df_t[i] = df_t[i].astype(float)
        
    plume.data = df_t
    plume.spe = spe
    return plume

#--------------Particle evolution from DBG file
def read_partDBG(filename):
    fid = open(filename,'r')
    data = fid.readlines()
    zoneind = [data.index(l) for l in data if l.startswith('      TIME')]
    zoneind = np.array(zoneind).astype(int)
    x = range(zoneind[0]+2, len(data))
    test = np.array([data[i].split() for i in x])
    dim = len(test[0])
    
    data_it = []
    for it in range(len(test)):
        if len(test[it])==dim:
            data_it.append(np.array(test[it]).astype(float))
    data_it = np.array(data_it)
    
    return data_it


#--------------Exit particles
def read_exit_particles(filename):
    fid = open(filename,'r')
    Lines = fid.readlines()
    data = Lines[2:]
    data_val = np.array([data[i].split() for i in range(len(data))])
    df = pd.DataFrame(data_val, columns = ["Xp","Yp","Zp","SPECIES","ZONE","PARTID","TIME","LOCATION"])
    df['Xp'] = df['Xp'].astype(float)
    df['Yp'] = df['Yp'].astype(float)
    df['Zp'] = df['Zp'].astype(float)
    df['SPECIES'] = df['SPECIES'].astype(int)
    df['ZONE'] = df['ZONE'].astype(int)
    df['PARTID'] = df['PARTID'].astype(int)
    df['TIME'] = df['TIME'].astype(float)
    return df

#--------------Spatial moments


#get in which cell are particles located
#       mesh: from mike dfs input file
#       x0: dataframe from particle plume reader (single time step)
# lower left cell at index (0,0)
def locate_cell_part(partloc,disc,mesh):
    
    partloc["ix"] = " "
    partloc["iy"] = " "
    partloc["iz"] = " "
    partloc["cellID"] = " "
    
    # xcoord = mesh.mesh_x[0,0,:]
    # ycoord = mesh.mesh_y[0,:,0]
    partloc.ix = np.floor(partloc["xp"]/disc.dx).astype(int)
    partloc.iy = np.floor(partloc["yp"]/disc.dy).astype(int)
    for ip in range(len(partloc)):
        # diff_x = xcoord-partloc.xp[ip]
        # partloc.ix[ip] = np.where(diff_x == max([n for n in diff_x if n<0]))[0][0]
        
        # diff_y = ycoord-partloc.yp[ip]
        # partloc.iy[ip] = np.where(diff_y == max([n for n in diff_y if n<0]))[0][0]
        
        if len(mesh.mesh_z)>2:
            zcoord = mesh.mesh_z[:,partloc.iy[ip],partloc.ix[ip]]
            diff_z = zcoord-partloc.zp[ip]
            partloc.iz[ip] = np.where(diff_z == max([n for n in diff_z if n<0]))[0][0]
        else:
            partloc.iz[ip] = 0
        
        partloc.cellID[ip] = partloc.ix[ip]*100000+partloc.iy[ip]*100+partloc.iz[ip]
        
    return partloc