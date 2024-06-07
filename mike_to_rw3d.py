# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:59:59 2023

@author: cvh
"""

import hydroinform.DFS as DFS
import numpy as np
import mikeio
import pandas as pd
import warnings
import struct
from scipy.io import FortranFile

import sys
sys.path.append(r'\\netapp1p\Dkmodel-hydro\Diverse\Python\dfstools')
from LTR_pfs import wel2df

# from mikeio_spatial_grid_geometry import _index_to_Grid3D

#------------------------------------------------------------------------------
#--------------get domain discretization
#-----OBS------(works for variable layer thickness)
def get_discretization(fp,file,pr):
    data = DFS.DFS3.from_file(fp)
    
    #variable layer thickness
    items = [item.name for item in data.items]
    to_ana = 'Thickness of computational layers in the saturated zone'
    item = items.index(to_ana) + 1 #get the item index for a specific item (remember +1!)
    dz = np.empty((data.number_of_layers, data.number_of_rows, data.number_of_columns))
    dz[:,:,:] = data.get_data(0, item)
    dz[dz<0] = 9999
    
    # dz2 = np.empty((nlay, nrow, ncol))
    # (tl,ts,meta)=RS_dfs.extract_dfs_data_with_meta(fp,item_name='Thickness of computational layers in the saturated zone', static=False)
    # dz2[:,:,:] = tl[0,:,:,:]
    # dz2[np.isnan(dz2)] = 9999
    
    class output:
        pass
    disc = output()
    disc.nrow = data.number_of_rows
    disc.ncol = data.number_of_columns
    disc.nlay = data.number_of_layers
    disc.dx = data.dx
    disc.dy = data.dy
    disc.dz = dz
    
    #print
    if pr=='T':
        print_gslib_3D(dz,file,'thickness of each layer','dz')
    
    return disc


#------------------------------------------------------------------------------
#--------------get mesh (lower left coordinates of each cell)
def get_mesh(disc,floor):
    mesh_x = np.zeros([disc.nlay, disc.nrow, disc.ncol+1], dtype=float)
    mesh_y = np.zeros([disc.nlay, disc.nrow+1, disc.ncol], dtype=float)
    mesh_z = np.zeros([disc.nlay+1, disc.nrow, disc.ncol], dtype=float)
    for ix in range(disc.ncol):
        for iy in range(disc.nrow):
            mesh_z[0,iy,ix] = floor.elevation[iy,ix]
            for iz in range(disc.nlay):
                mesh_x[iz,iy,ix+1] = mesh_x[iz,iy,ix]+disc.dx
                mesh_y[iz,iy+1,ix] = mesh_y[iz,iy,ix]+disc.dy
                mesh_z[iz+1,iy,ix] = mesh_z[iz,iy,ix]+disc.dz[iz,iy,ix]
    
    class output:
        pass
    mesh = output()
    mesh.mesh_x = mesh_x
    mesh.mesh_y = mesh_y
    mesh.mesh_z = mesh_z
    
    return mesh
    
#------------------------------------------------------------------------------
#--------------get floor elevation
def get_floor_elevation(fp,file,pr):
    data = DFS.DFS3.from_file(fp)
    
    items = [item.name for item in data.items]
    to_ana = 'Lower level of computational layers in the saturated zone'
    item = items.index(to_ana) + 1
    lowlevel = np.empty((data.number_of_layers, data.number_of_rows, data.number_of_columns))
    lowlevel[:,:,:] = data.get_data(0, item)
    # lowlevel[lowlevel<0] = 0

    class output:
        pass
    floor = output()
    floor.nrow = data.number_of_rows
    floor.ncol = data.number_of_columns
    floor.elevation = lowlevel[0,:,:]
    
    #print
    if pr=='T':
        print_gslib_2D(floor.elevation,file,'floor elevation','z0')
    
    return floor


#------------------------------------------------------------------------------
#--------------get inactive cells
def get_inactive_cells(fp,file,pr):
    dfs = mikeio.read(fp)
    bc = bc2 = dfs.Boundary_conditions_for_the_saturated_zone
    bc_val = bc2.values
    if len(bc.shape)==4:
        bc_val = bc_val[0,:,:,:]
    if len(bc.shape)==3:
        bc_val = bc_val[:,:,:]
    InactCell = bc_val
    InactCell[~np.isnan(InactCell)] = 1 #all non-NaN values are active
    InactCell[np.isnan(InactCell)] = 0
    
    #print
    if pr=='T':
        print_gslib_3D(InactCell,file,'Inactive cells','InactCell')
    
    return InactCell


#------------------------------------------------------------------------------
#--------------get flow from/to rivers
def get_Qsink(fp,sink_type,file,it,pr):
    dfs = mikeio.read(fp)
    
    if sink_type == 'river':
        data = dfs.SZ_exchange_flow_with_river
        to_ana = 'SZ exchange flow with river'
        header = 'flux from/to rivers'
        varname = 'Qriver'
    if sink_type == 'drain':
        data = dfs.SZ_drainage_flow_from_point
        to_ana = 'SZ drainage flow from point'
        header = 'flux from/to drains'
        varname = 'Qdrain'
            
    #process data for each time step
    if len(data.shape)==4:
        Qsink = data[it,:,:,:].values*86400*365 #m3/s to m3/year
        #print
        if pr=='T':   
            print_gslib_3D(Qsink,file,header,varname)
    if len(data.shape)==3: 
        Qsink = data[data.shape[0],:,:].values*86400*365 #m3/s to m3/year
        #print
        if pr=='T':   
            print_gslib_2D(Qsink,file,header,varname)
    
    # data = DFS.DFS3.from_file(fp)
    # items = [item.name for item in data.items]
    # item = items.index(to_ana) + 1
    
    # nt = data.number_of_time_steps
    # dt = data.time_step_size
    # dts = dt.total_seconds()/86400 #seconds to days
    # times = np.linspace(0,nt*dts,nt+1)
    # itp = 0
    
    # for it in range(nt-1):
    #     Qsink_lay = data.get_data(it+1, item)
    #     Qsink = np.empty((data.number_of_layers, data.number_of_rows, data.number_of_columns))
    #     for iz in range(data.number_of_layers):
    #         Qsink_layer = Qsink_lay[iz]
    #         Qsink_layer.mask = np.ma.nomask
    #         Qsink[iz,:,:] = Qsink_layer*86400*365 #m3/s to m3/year
            
    #     print_gslib_3D_transient(Qsink,file,header,varname,it+1,times,itp)
    #     itp = itp+1 
        
    # print_gslib_3D(Qsink,file,header,varname)

    return Qsink


#------------------------------------------------------------------------------
#--------------get groundwater (Darcy) fluxes
def get_flux(fp,fp_uz,direction,file,binary,it,pr):
    dfs = mikeio.read(fp)
    # dfs_uz = mikeio.read(fp_uz)
    if direction == 'x':
        data = dfs.groundwater_flux_in_x_direction
        header = 'flux in x-direction'
        varname = 'qx'
    if direction == 'y':
        data = dfs.groundwater_flux_in_y_direction
        header = 'flux in y-direction'
        varname = 'qy'
    if direction == 'z':
        data = dfs.groundwater_flux_in_z_direction
        header = 'flux in z-direction'
        varname = 'qz'
    
    if len(data.shape)==4:
        Qp = data[it-1,:,:,:].values/1000*365 #mm/d to m/year
        #readjust the dimension of the flux matrix depending on the direction
        if direction == 'x':
            # Q_add = Qp[:, :, data.geometry.nx-1]
            # Qprint = np.concatenate([Qp, Q_add[:,:,None]],axis=2)
            Q_add = Qp[:,:,0]
            Qprint = np.concatenate([Q_add[:,:,None], Qp],axis=2)
            Qprint = np.nan_to_num(Qprint) #replace nan values by 0

        if direction == 'y':
            # Q_add = Qp[:, data.geometry.ny-1, :]
            # Qprint = np.concatenate([Qp, Q_add[:,None,:]],axis=1)
            Q_add = Qp[:,0,:]
            Qprint = np.concatenate([Q_add[:,None,:], Qp],axis=1)
            Qprint = np.nan_to_num(Qprint) #replace nan values by 0

        if direction == 'z':
            # Q_add = -dfs_uz.Total_recharge_to_SZ_pos_down_.values[it-1,:,:]/1000*365
            # Qprint = np.concatenate([Qp, Q_add[None,:,:]],axis=0)
            Q_add = Qp[0,:,:]
            Q_add[~np.isnan(Q_add)]=0
            Qprint = np.concatenate([Q_add[None,:,:], Qp],axis=0)
        
        #print
        if pr=='T':   
            print_gslib_3D(Qprint,file,header,varname)
            
    if len(data.shape)==3:
        Qp = data[data.shape[0]-1,:,:].values/1000*365 #mm/d to m/year
        if direction == 'x':
            # Q_add = Qp[:, data.geometry.nx-1]
            # Qprint = np.concatenate([Qp, Q_add[:,None]],axis=1)
            Q_add = Qp[:,0]
            Qprint = np.concatenate([Q_add[:,None], Qp],axis=1)
            Qprint = np.nan_to_num(Qprint)
            
            if pr=='T':   
                print_gslib_2D(Qprint,file,header,varname)    
                
        if direction == 'y':
            # Q_add = Qp[data.geometry.ny-1, :]
            # Qprint = np.concatenate([Qp, Q_add[None,:]],axis=0)
            Q_add = Qp[0,:]
            Qprint = np.concatenate([Q_add[None,:], Qp],axis=0)
            Qprint = np.nan_to_num(Qprint)
            
            if pr=='T':   
                print_gslib_2D(Qprint,file,header,varname)
                
        if direction == 'z':
            Q_bot = np.copy(Qp)
            Q_bot[~np.isnan(Q_bot)]=0
            Q_bot = Q_bot[np.newaxis,:,:]
            Qprint = np.concatenate([Q_bot, Qp[None,:,:]],axis=0)
            # Q_add = -dfs_uz.Total_recharge_to_SZ_pos_down_.values[it-1,:,:]/1000*365
            # Qp = Qp[np.newaxis,:,:]
            # Qprint = np.concatenate([Qp, Q_add[None,:,:]],axis=0)
            if pr=='T':   
                print_gslib_3D(Qprint,file,header,varname)
        

    
    return Qprint
    
    #print
    # if binary == 0:
    #     print_gslib_3D_transient(Qprint,file,header,varname,it+1,times,itp)
    # if binary == 1:
    #     print_binary_3D_transient(Qprint,file,header,varname,it+1,times,itp)
    # itp = itp+1

    
    #FOR TRANSIENT
    
    # #get times
    # nt = dfs.n_timesteps
    # dt = dfs.timestep/86400/365 #time step from s in years
    # times = np.linspace(0,nt*dt,nt+1)
    
    # #process data for each time step
    # itp = 0
    # for it in range(nt-1):
    #     Qp = data[it+1,:,:,:].values/1000*365 #mm/d to m/year
    #     #readjust the dimension of the flux matrix depending on the direction
    #     if direction == 'x':
    #         Q_add = Qp[:, :, data.geometry.nx-1]
    #         Qprint = np.concatenate([Qp, Q_add[:,:,None]],axis=2)
    #     if direction == 'y':
    #         Q_add = Qp[:, data.geometry.ny-1, :]
    #         Qprint = np.concatenate([Qp, Q_add[:,None,:]],axis=1)
    #     if direction == 'z':
    #         Q_add = Qp[data.geometry.nz-1, :, :]
    #         Qprint = np.concatenate([Qp, Q_add[None,:,:]],axis=0)
    #     #print    
    #     if binary == 0:
    #         print_gslib_3D_transient(Qprint,file,header,varname,it+1,times,itp)
    #     if binary == 1:
    #         print_binary_3D_transient(Qprint,file,header,varname,it+1,times,itp)
    #     itp = itp+1


    # data = DFS.DFS3.from_file(fp) 
    # items = [item.name for item in data.items]
    # if direction == 'x':
    #     to_ana = 'groundwater flux in x-direction'
    #     header = 'flux in x-direction'
    #     varname = 'qx'
    # if direction == 'y':
    #     to_ana = 'groundwater flux in y-direction'
    #     header = 'flux in y-direction'
    #     varname = 'qy'
    # if direction == 'z':
    #     to_ana = 'groundwater flux in z-direction'
    #     header = 'flux in z-direction'
    #     varname = 'qz'
    # item = items.index(to_ana) + 1
    # nt = data.number_of_time_steps
    # dt = data.time_step_size
    # dts = dt.total_seconds()/86400/365 #seconds to year
    # times = np.linspace(0,nt*dts,nt+1)
    # itp = 0
    # for it in range(nt-1):
    #     Qit = data.get_data(it+1, item)
    #     Qp = np.empty((data.number_of_layers, data.number_of_rows, data.number_of_columns))
    #     for iz in range(data.number_of_layers):
    #         Q_layer = Qit[iz]
    #         Q_layer.mask = np.ma.nomask
    #         Qp[iz,:,:] = Q_layer/1000*365 #mm/d to m/year
        
    #     #readjust the dimension of the flux matrix depending on the direction
    #     if direction == 'x':
    #         Q_add = Qp[:, :, data.number_of_columns-1]
    #         Qprint = np.concatenate([Qp, Q_add[:,:,None]],axis=2)
    #     if direction == 'y':
    #         Q_add = Qp[:, data.number_of_rows-1, :]
    #         Qprint = np.concatenate([Qp, Q_add[:,None,:]],axis=1)
    #     if direction == 'z':
    #         Q_add = Qp[data.number_of_layers-1, :, :]
    #         Qprint = np.concatenate([Qp, Q_add[None,:,:]],axis=0)
            
    #     if binary == 0:
    #         print_gslib_3D_transient(Qprint,file,header,varname,it+1,times,itp)
    #     if binary == 1:
    #         print_binary_3D_transient(Qprint,file,header,varname,it+1,times,itp)
    #     itp = itp+1 

#------------------------------------------------------------------------------
#--------------get groundwater (volumetric) flow
def get_flow(fp,direction,file,it,pr):
    dfs = mikeio.read(fp)
    if direction == 'x':
        data = dfs.groundwater_flux_in_x_direction
        header = 'flow in x-direction'
        varname = 'qx'
    if direction == 'y':
        data = dfs.groundwater_flux_in_y_direction
        header = 'flow in y-direction'
        varname = 'qy'
    if direction == 'z':
        print("error: no flows in z-direction, only fluxes")
        sys.exit()
    
    Qp = data[it-1,:,:,:].values*60*60*24*365 #m3/s to m3/year
    #readjust the dimension of the flux matrix depending on the direction
    if direction == 'x':
        Q_add = Qp[:, :, data.geometry.nx-1]
        Qprint = np.concatenate([Qp, Q_add[:,:,None]],axis=2)
    if direction == 'y':
        Q_add = Qp[:, data.geometry.ny-1, :]
        Qprint = np.concatenate([Qp, Q_add[:,None,:]],axis=1)
    
    #print
    if pr=='T':   
        print_gslib_3D(Qprint,file,header,varname)
    
    return Qprint

#------------------------------------------------------------------------------
#--------------get porosity
def get_porosity(fp,file,pr):
    dfs = mikeio.read(fp)
    poro = dfs.Porosity_in_the_saturated_zone[0].values
    if pr:
        header = 'porosity in the saturated zone'
        varname = 'poro'
        if len(poro.shape)==3:
            print_gslib_3D(poro,file,header,varname)
        if len(poro.shape)==2:
            print_gslib_2D(poro,file,header,varname)        
    
    return poro


#------------------------------------------------------------------------------
#--------------get head elevation in saturated zone
def get_heads(fp,file,it,pr):
    dfs = mikeio.read(fp)
    header = 'head elevation in saturated zone'
    varname = 'hds'
    data = dfs.head_elevation_in_saturated_zone
    if len(data.shape)==4:
        heads = data.values[it-1,:,:,:]    
        #print
        if pr=='T':   
            print_gslib_3D(heads,file,header,varname)
            
    if len(data.shape)==3:
        heads = data.values[it-1,:,:]    
        #print
        if pr=='T':   
            print_gslib_2D(heads,file,header,varname)        
            
    return heads


#------------------------------------------------------------------------------
#--------------get total recharge
def get_total_recharge(fp,file,it,pr):
    varname = 'Total recharge to SZ (pos.down)'
    dfs = mikeio.read(fp,items=varname)
    mult = -1/1000*365*100*100 
    dfs[varname].values = dfs[varname].values*mult #conversion of units, if necessary
    
    class output:
        pass
    rch = output()
    rch.values = dfs[varname].values[it,:,:]
    
    #print
    if pr=='T':
        print_gslib_2D(rch.values,file,'total_recharge','rech')
    
    return rch


#------------------------------------------------------------------------------
#--------------get total recharge
def get_seepage_flow(fp,file,it,pr):
    varname = 'seepage flow SZ -overland'
    dfs = mikeio.read(fp,items=varname)
    mult = 1/1000*365*100*100 
    dfs[varname].values = dfs[varname].values*mult #conversion of units, if necessary
    
    class output:
        pass
    rch = output()
    rch.values = dfs[varname].values[it,:,:]
    
    #print
    if pr=='T':
        print_gslib_2D(rch.values,file,'total_recharge','rech')
    
    return rch

#------------------------------------------------------------------------------
#--------------get total recharge
def get_qUZtoSZ(fp,file,it,pr):
    varname = 'exchange between UZ and SZ (pos.up)'
    dfs = mikeio.read(fp,items=varname)
    mult = 1/1000*365*100*100 
    dfs[varname].values = dfs[varname].values*mult #conversion of units, if necessary
    
    class output:
        pass
    rch = output()
    rch.values = dfs[varname].values[it,:,:]
    
    #print
    if pr=='T':
        print_gslib_2D(rch.values,file,'total_recharge','rech')
    
    return rch

#------------------------------------------------------------------------------
#--------------get all information for all wells
def get_wellinfo(well_file,dfs0_well,dfs3_flow,logfile,mesh,file_out,pr):  
    
    class output:
        pass
    well = output()
    
    dfs = mikeio.read(dfs3_flow) 
    times = dfs.time-dfs.time[0] #set first time flag as t=0
    times = times.days/365 #days to years
    
    #get information from log file (IDs of wells in the studied area)
    df_welllog = process_WMlogfile(logfile)
    
    #get prescribed extraction rate time series from dfs0 file
    df_Qpresc_timeseries = process_well_dfs0(dfs0_well,df_welllog)
    
    #get location and screen top/bottom (for specific set of wells) &PRINT
    df_well_info = process_wellfile(well_file,df_welllog,dfs)
    if pr=='T':
        print_well_location(df_well_info,file_out)
    
    #locate in each cell are the wells
    df_well_info, df_well_cell = locate_cell_well(df_well_info,mesh)
    
    #for each timestep (in dfs3_flow), get the effective extraction rate per cell (accounting for the relative contribution of each well) &PRINT
    df_well_cell = get_Qeff_timestep(df_well_cell,dfs,df_Qpresc_timeseries,times,file_out)        
    
    return well


#--------------series of required functions
#Process the WM_Print log file: get the well IDs
def process_WMlogfile(logfile):
    file = open(logfile)
    lines = file.readlines()
    start = np.asarray([lines.index(l) for l in lines if l.startswith('Well')])+2
    end   = np.asarray([lines.index(l) for l in lines if l.startswith('    TOTAL:')])-2
    data = np.asarray(lines[start[0]:end[0]])
    wellinfo = np.empty([len(data),5])
    wellID_dict = {}
    for iwll in range(len(data)):
        datawll = data[iwll].split()
        wellinfo[iwll,:] = datawll[0:5]
        wellID_dict[iwll] = datawll[6]
    wellID = list(wellID_dict.values())    
    
    #reorganize in a dataframe
    df_welllog = pd.DataFrame(wellinfo, columns = ['Well','Filter','Specified_Q','Applied_Q','Efficiency'])
    df_welllog['WellID'] = wellID
    df_welllog['Filter'] = df_welllog['Filter'].astype(int)
    df_welllog['WellID_dfs0'] = df_welllog['WellID'] + "_" + df_welllog['Filter'].astype(str)
    df_welllog.loc[df_welllog['Filter'] == 1, 'WellID_dfs0'] = df_welllog['WellID']
    
    return df_welllog


#Process the dfs0 file: get prescribed extraction rates time series
def process_well_dfs0(dfs0_file,df_welllog):
    ds = mikeio.read(dfs0_file)
    df = ds.to_dataframe()
    df_Qpresc_timeseries = pd.DataFrame(index=df.index,columns=range(0))
    wellIDunique = list(set(df_welllog.WellID))
    for iwll in range(len(wellIDunique)): 
        df_singlewell = df[df.columns[pd.Series(df.columns).str.startswith(wellIDunique[iwll])]]
        df_Qpresc_timeseries = df_Qpresc_timeseries.join(df_singlewell)
    df_Qpresc_timeseries = df_Qpresc_timeseries.reset_index()
    df_Qpresc_timeseries = df_Qpresc_timeseries.rename(columns={"index": "time"})
    
    return df_Qpresc_timeseries


#Process the well file: get well location (x,y)
def process_wellfile(well_file,df_welllog,dfs):
    df_well_DK = wel2df(well_file) #slow to debug!
    df_well_info = df_well_DK[df_well_DK['dkm_id'].isin(df_welllog['WellID'])]
    df_well_info['xcor'] = df_well_info['xcor']-dfs.geometry.origin[0] #project location to a grid with lower left corner coordinates [0,0]
    df_well_info['ycor'] = df_well_info['ycor']-dfs.geometry.origin[1] #project location to a grid with lower left corner coordinates [0,0]
    df_well_info = df_well_info.reset_index()
    df_well_info['rw'] = 1.0
    
    return df_well_info


#get in which cell are wells located
def locate_cell_well(df_well_info,mesh):
    df_well_info["ix"] = " "
    df_well_info["iy"] = " "
    df_well_info["iz_top"] = " "
    df_well_info["iz_bottom"] = " "
    
    xcoord = mesh.mesh_x[0,0,:]
    ycoord = mesh.mesh_y[0,:,0]
    
    df_well_cell = pd.DataFrame(columns=['ID','ix','iy','iz','cellind','Qtotcell','Qprescribed','Qeff'])
    
    iwellcell = 0
    for iwll in range(len(df_well_info)):
        diff_x = xcoord-df_well_info.xcor[iwll]
        df_well_info.ix[iwll] = np.where(diff_x == max([n for n in diff_x if n<0]))[0][0]
        
        diff_y = ycoord-df_well_info.ycor[iwll]
        df_well_info.iy[iwll] = np.where(diff_y == max([n for n in diff_y if n<0]))[0][0]
        
        zcoord = mesh.mesh_z[:,df_well_info.iy[iwll],df_well_info.ix[iwll]]
        diff_z = zcoord-df_well_info.top[iwll]
        df_well_info.iz_top[iwll] = np.where(diff_z == max([n for n in diff_z if n<0]))[0][0]
        diff_z = zcoord-df_well_info.bottom[iwll]
        df_well_info.iz_bottom[iwll] = np.where(diff_z == max([n for n in diff_z if n<0]))[0][0]
        
        #create unique indexes of cells where screen is located
        for ilay in range(df_well_info.iz_bottom[iwll],df_well_info.iz_top[iwll]+1):
            df_well_cell.at[iwellcell,'ix'] = df_well_info.ix[iwll]
            df_well_cell.at[iwellcell,'iy'] = df_well_info.iy[iwll]
            df_well_cell.at[iwellcell,'iz'] = ilay
            df_well_cell.at[iwellcell,'ID'] = df_well_info.dkm_id[iwll]
            df_well_cell.at[iwellcell,'cellind'] = df_well_info.ix[iwll]*10000000+df_well_info.iy[iwll]*10000+df_well_info.iz_top[iwll]*100+ilay
            iwellcell = iwellcell+1
            
    return df_well_info, df_well_cell


#for each timestep (in dfs3_flow), get the effective extraction rate per cell (accounting for the relative contribution of each well) &PRINT
def get_Qeff_timestep(df_well_cell,dfs,df_Qpresc_timeseries,times,file_out):
    for it in range(dfs.n_timesteps-1):
        t = dfs.time[it+1]
        tprev = dfs.time[it]
        Qwell_eff = dfs.groundwater_extraction.values[it+1,:,:,:] #extraction fluxes averaged over the last timestep
        
        df_Qpresc_t = df_Qpresc_timeseries.iloc[(df_Qpresc_timeseries['time']-t).abs().argsort()[:1]] #get nearest time from the prescribed extraction time series
        if df_Qpresc_t.time[df_Qpresc_t.index[0]]>t: #get the lower nearest time step
            df_Qpresc_t = df_Qpresc_timeseries.loc[df_Qpresc_t.index-1]
            
        df_Qpresc_tprev = df_Qpresc_timeseries.iloc[(df_Qpresc_timeseries['time']-tprev).abs().argsort()[:1]] #get nearest time from the prescribed extraction time series
        if df_Qpresc_tprev.time[df_Qpresc_tprev.index[0]]>tprev: #get the lower nearest time step
            df_Qpresc_tprev = df_Qpresc_timeseries.loc[df_Qpresc_tprev.index-1]
        
        df_Qpresc_dt = df_Qpresc_timeseries[df_Qpresc_tprev.index[0]:df_Qpresc_t.index[0]+1]
        # df_Qpresc_dt_transposed = df_Qpresc_dt.T
        
        for iwll in range(len(df_well_cell)):
            df_well_cell.Qtotcell[iwll] = Qwell_eff[df_well_cell.iz[iwll],df_well_cell.iy[iwll],df_well_cell.ix[iwll]]*60*60*24*365 #convert m3/s to m3/yr
            df_well_cell.Qprescribed[iwll] = df_Qpresc_dt[df_well_cell.ID[iwll]].mean() #df_Qpresc_t[df_well_cell.ID[iwll]].values[0]
            
        df_well_cell_shared = df_well_cell.loc[df_well_cell.duplicated(subset='cellind', keep=False)]
        df_well_cell_notshared = df_well_cell.drop_duplicates(subset='cellind', keep=False)
        df_well_cell_notshared["Qeff"] = df_well_cell_notshared["Qtotcell"]
        df_well_cell["Qeff"] = df_well_cell_notshared["Qeff"]
        
        wellcellunique_shared = list(set(df_well_cell_shared.cellind))
        for icell in range(len(wellcellunique_shared)):
            df_well_icell = df_well_cell_shared.loc[df_well_cell_shared['cellind']==wellcellunique_shared[icell]]
            df_well_icell['Qeff'] = df_well_icell['Qtotcell']*df_well_icell['Qprescribed']/df_well_icell['Qprescribed'].sum()
            df_well_cell.loc[df_well_icell.index.values] = df_well_icell.loc[df_well_icell.index.values]
        
        df_well_cell['Qeff'] = df_well_cell['Qeff'].fillna(0)
        
        print_well_flux_time(df_well_cell,file_out,it+1,times)

        return df_well_cell
    
#------------------------------------------------------------------------------
#--------------print file
#------------------------------------------------------------------------------
def print_gslib_3D(var,file,header,varname):
    f = open(file, 'w')
    f.write(f"{header}\n")
    f.write("1\n")
    f.write(f"{varname}\n")
    dim = var.shape
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                print(var[i,j,k], file = f)
    f.close()
    
    
def print_gslib_2D(var,file,header,varname):
    f = open(file, 'w')
    f.write(f"{header}\n")
    f.write("1\n")
    f.write(f"{varname}\n")
    dim = var.shape
    for i in range(dim[0]):
        for j in range(dim[1]):
                print(var[i,j], file = f)
    f.close()
    
    
def print_gslib_3D_transient(var,file,header,varname,it,times,itp):
    if itp==0:
        f = open(file, 'w')
        f.write(f"{header}\n")
        f.write("1\n")
    else:
        f = open(file, 'a')
    f.write(f"{varname}; time = {times[it]}\n")
    dim = var.shape
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                print(var[i,j,k], file = f)
    if it==len(times)-1:
        f.close()

def print_binary_3D_transient(var,file,header,varname,it,times,itp):
    f = open(file, 'wb')
    dim = var.shape
    headpr = [it,times[it],dim[2],dim[1],dim[0]]
    headpa = pack_all(headpr)
    f.write(headpa)
    f.write('\n'.encode('utf-8'))
    # varbin = bytearray(np.float32(var))
    # f.write(varbin)
    
    for i in range(dim[0]):
        varlay = bytearray(np.float32(var[i,:,:]))
        f.write(varlay)
        f.write('\n'.encode('utf-8'))
    if it==len(times)-1:
        f.close()
        
def pack_all(lst):
    fmt = ''.join('i' if isinstance(x, int) else 'f' for x in lst)
    return struct.pack(fmt, *lst)

def print_well_location(var,file):
    var = var.sort_values(by=['dkm_id'])
    var = var.reset_index()
    
    f = open(file, 'w')
    f.write("Well information (ID, x, y, rw, zbot, ztop and Qw for each time step)\n")
    nwll = len(var)
    print(nwll, file = f)
    for iwll in range(nwll):
        print(var['dkm_id'].loc[iwll], var['xcor'].loc[iwll], var['ycor'].loc[iwll], var['rw'].loc[iwll], var['bottom'].loc[iwll], var['top'].loc[iwll], "10", "F", file = f)
    f.write("\n")

    
def print_well_flux_time(var,file,it,times):
    f = open(file, 'a')
    f.write(f"Qw; time = {times[it]}\n")
    wellIDunique = list(set(var.ID))
    wellIDunique = sorted(wellIDunique)
    for iwll in range(len(wellIDunique)):
        df_wll = var.loc[var['ID']==wellIDunique[iwll]]
        df_wll = df_wll.reset_index()
        for ilay in range(len(df_wll)):
            print(df_wll['ID'].loc[ilay], df_wll['ix'].loc[ilay], df_wll['iy'].loc[ilay], df_wll['iz'].loc[ilay], round(df_wll['Qeff'].loc[ilay],4), file = f)
        
    if it==len(times)-1:
        f.close()