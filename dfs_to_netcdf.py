# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 14:27:18 2023

@author: cvh
"""
import mikeio
import numpy as np
import xarray as xr
import gc

#:-----------------------------------------------------------------------------
def _add_column_to_DataArray_(arr):
    val = arr.to_numpy()
    val_add = val[:,:,:,0] #add column with same values than the first column
    val_mod = np.concatenate([val_add[:,:,:,None], val],axis=3)
    np.nan_to_num(val_mod,copy=False)
    dx = arr.x[1] - arr.x[0]
    time = arr.time.to_numpy()
    x = arr.x.to_numpy()
    y = arr.y.to_numpy()
    z = arr.z.to_numpy()
    x = np.append(x, x[-1]+dx)
    new_arr = xr.DataArray(data=val_mod, 
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr.name,
                           attrs=dict(name=arr.attrs['name']))
    del arr
    gc.collect()
    return new_arr


#:-------------
def _add_row_to_DataArray_(arr):
    val = arr.to_numpy()
    val_add = val[:,:,0,:]
    val_mod = np.concatenate([val_add[:,:,None,:], val],axis=2)
    np.nan_to_num(val_mod,copy=False)
    dy = arr.y[1] - arr.y[0]
    time = arr.time.to_numpy()
    x = arr.x.to_numpy()
    y = arr.y.to_numpy()
    z = arr.z.to_numpy()
    y = np.append(y, y[-1]+dy)
    new_arr = xr.DataArray(data=val_mod, 
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr.name,
                           attrs=dict(name=arr.attrs['name']))
    del arr
    gc.collect()
    return new_arr


#:-------------
def _add_layer_to_DataArray_(arr):
    val = arr.to_numpy()
    val_add = val[:,0,:,:]
    val_mod = np.concatenate([val_add[:,None,:,:], val],axis=1)
    del val,val_add
    np.nan_to_num(val_mod,copy=False)
    dz = arr.z[1] - arr.z[0]
    time = arr.time.to_numpy()
    x = arr.x.to_numpy()
    y = arr.y.to_numpy()
    z = arr.z.to_numpy()
    z = np.append(z, z[-1]+dz)
    new_arr = xr.DataArray(data=val_mod, 
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr.name,
                           attrs=dict(name=arr.attrs['name']))
    del arr,val_mod
    gc.collect()
    return new_arr


#:-------------
def _DataArray_it_(arr,it_start,it_end):
    val = np.copy(arr.values) #arr.to_numpy()
    time = np.copy(arr.time) #.to_numpy()
    
    time_it = time[it_start:it_end+1]
    val_it = val[it_start:it_end+1,:,:,:]
    
    x = arr.x.to_numpy()
    y = arr.y.to_numpy()
    z = arr.z.to_numpy()
    new_arr = xr.DataArray(data=val_it, 
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time_it),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr.name,
                           attrs=dict(name=arr.attrs['name']))
    del arr
    gc.collect()
    return new_arr
    

#:-------------
def _top2D_to_3D_DataArray_(arr,nlay):
    x = arr.x.to_numpy()
    y = arr.y.to_numpy()
    z = np.linspace(0, nlay-1, num=nlay)
    nx, ny, nz = len(x), len(y), len(z)
    
    time = arr.time.to_numpy()
    nt = len(arr[:,0,0])
    
    val3D = np.empty([nt,nlay,ny,nx])
    val2D = arr.to_numpy()
    val3D[:,nlay-1,:,:] = val2D
    val3D[val3D<0] = 0
    
    new_arr = xr.DataArray(data=val3D,
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr.name,
                           attrs=dict(name=arr.attrs['name']))
    del arr
    gc.collect()
    return new_arr

#:-------------
def _topvalues_only_3D_DataArray_(arr,nlay):
    x = arr.x.to_numpy()
    y = arr.y.to_numpy()
    z = np.linspace(0, nlay-1, num=nlay)
    nx, ny, nz = len(x), len(y), len(z)
    
    time = arr.time.to_numpy()
    nt = len(arr[:,0,0])
    
    val3D = np.empty([nt,nlay,ny,nx])
    valin = arr.to_numpy()
    val3D[:,nlay-1,:,:] = valin[:,nlay-1,:,:]
    val3D[val3D<0] = 0
    
    new_arr = xr.DataArray(data=val3D,
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr.name,
                           attrs=dict(name=arr.attrs['name']))
    del arr
    gc.collect()
    return new_arr

def _add_two_3D_DataArray_(arr1,arr2):
    x = arr1.x.to_numpy()
    y = arr1.y.to_numpy()
    z = arr1.z.to_numpy()
    time = arr1.time.to_numpy()
    
    val = arr1.to_numpy() + arr2.to_numpy()
    
    new_arr = xr.DataArray(data=val,
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr1.name,
                           attrs=dict(name=arr1.attrs['name']))
    del arr1,arr2
    gc.collect()
    return new_arr


#:-------------
def get_time_discretization(fp,file,it_start=-1,it_end=-1):
    dfs = mikeio.read(fp,time=slice(it_start,it_end))
    time = dfs.time.to_numpy()
    #time_it = time[it_start:it_end+1]
    nt = len(time)
    
    dt = dfs.timestep/60/60/24/365
    dt_disc = np.full((nt,), dt)
    
    _print_time_discretization_(dt_disc,file,nt,'time_disc')
    del dfs
    gc.collect()
    
    
def _print_time_discretization_(var,file,nt,header):
    f = open(file, 'w')
    f.write(f"{header}\n")
    f.write(f"{nt}\n")
    dim = var.shape
    for i in range(dim[0]):
        print(var[i], file = f)
    f.close()


#:-----------------------------------------------------------------------------
#:-------------get netCDF file
def get_netcdf(fp,varname,mult,file,it_start=-1,it_end=-1,nlay=-1):
    dfs = mikeio.read(fp,items=varname,time=slice(it_start,it_end))
    xr_da = dfs[varname].to_xarray()
    del dfs
    
    #convert 2D dataarray to 3D dataarray
    if varname == 'Total recharge to SZ (pos.down)':
        xr_da = _top2D_to_3D_DataArray_(xr_da,nlay)
    
    if mult!=1: 
        xr_da.values = xr_da.values*float(mult)
        # for it in range(len(xr_da[:,0,0,0])):
        #     xr_da[it,:,:,:] = xr_da[it,:,:,:]*mult #conversion of units, if necessary

    #extract top layer values    
    if varname == 'groundwater flux in z-direction' and nlay>0:
        xr_da = _topvalues_only_3D_DataArray_(xr_da,nlay)
        
    #select specific time if steady state problem
    #if it_start>0:
    #    xr_da = _DataArray_it_(xr_da,it_start,it_end)
    
    #replace all NaN with 0
    np.nan_to_num(xr_da.values,copy=False)
    # for it in range(len(xr_da[:,0,0,0])):
    #    xr_da[it,:,:,:] = xr_da[it,:,:,:].fillna(0)
    
    #add values, if necessary (groundwater fluxes only)
    if varname == 'groundwater flux in x-direction':
        xr_da = _add_column_to_DataArray_(xr_da)
        xr_da.to_netcdf(file, format='NETCDF3_64BIT')
        
    elif varname == 'groundwater flux in y-direction':
        xr_da = _add_row_to_DataArray_(xr_da)
        xr_da.to_netcdf(file, format='NETCDF3_64BIT')
        
    elif varname == 'groundwater flux in z-direction':
        xr_da = _add_layer_to_DataArray_(xr_da)
        xr_da.to_netcdf(file, format='NETCDF3_64BIT')
        
    else:
        xr_da.to_netcdf(file, format='NETCDF3_64BIT')
        
    del xr_da #,xr_da_mod
    gc.collect()
        
        
def get_netcdf_add(data1,data2,file,it_start=-1,it_end=-1,nlay=-1):
    dfs1 = mikeio.read(data1[0],items=data1[1])
    dfs1[data1[1]].values = dfs1[data1[1]].values*data1[2] #conversion of units, if necessary
    xr_da1 = dfs1[data1[1]].to_xarray()
    xr_da1 = _top2D_to_3D_DataArray_(xr_da1,nlay)
    
    dfs2 = mikeio.read(data2[0],items=data2[1])
    dfs2[data2[1]].values = dfs2[data2[1]].values*data2[2] #conversion of units, if necessary
    xr_da2 = dfs2[data2[1]].to_xarray()
    xr_da2 = _top2D_to_3D_DataArray_(xr_da2,nlay)
    
    xr_da = _add_two_3D_DataArray_(xr_da1,xr_da2)
    
    #select specific time if steady state problem
    if it_start>0:
        xr_da = _DataArray_it_(xr_da,it_start,it_end)
    
    #replace all NaN with 0
    xr_da = xr_da.fillna(0)
    
    xr_da.to_netcdf(file, format='NETCDF3_64BIT')
    
    
#:-----------------------------------------------------------------------------
#:-------------get netCDF file for horizontal fluxes from flows
def get_netcdf_q_from_Q(fp,fph,fpre,varname,mult,file,it_start=-1,it_end=-1,nlay=-1):
    #read flows
    dfs = mikeio.read(fp,items=varname,time=slice(it_start,it_end))
    if varname == 'groundwater flow in x-direction':
        cross_length = dfs.geometry.dy
        newname = 'groundwater flux in x-direction'
    if varname == 'groundwater flow in y-direction':
        cross_length = dfs.geometry.dx
        newname = 'groundwater flux in y-direction'
    xr_da = dfs[varname].to_xarray()
    if mult!=1: 
        xr_da.values = xr_da.values*float(mult)
    del dfs
    
    dz_sat = get_dz_sat(fph,fpre,it_start,it_end)
    
    xr_da = from_Q_to_q(xr_da,dz_sat,cross_length,newname)
        
    #replace all NaN with 0
    np.nan_to_num(xr_da.values,copy=False)
    
    #add values, if necessary (groundwater fluxes only)
    if varname == 'groundwater flux in x-direction' or varname == 'groundwater flow in x-direction':
        xr_da = _add_column_to_DataArray_(xr_da)
        xr_da.to_netcdf(file, format='NETCDF3_64BIT')
        
    elif varname == 'groundwater flux in y-direction' or varname == 'groundwater flow in y-direction':
        xr_da = _add_row_to_DataArray_(xr_da)
        xr_da.to_netcdf(file, format='NETCDF3_64BIT')
                
    else:
        xr_da.to_netcdf(file, format='NETCDF3_64BIT')
        
    del xr_da #,xr_da_mod

#:-----------------------------------------------------------------------------
#:-------------get saturated thickness of all cells
def get_dz_sat(fph,fpre,it_start,it_end):
    #read heads
    varname_hds = 'head elevation in saturated zone'
    dfs_hds = mikeio.read(fph,items=varname_hds,time=slice(it_start,it_end))
    heads = dfs_hds.head_elevation_in_saturated_zone.values
    nx = dfs_hds.geometry.nx
    ny = dfs_hds.geometry.ny
    nz = dfs_hds.geometry.nz
    nt = dfs_hds.n_timesteps
    del dfs_hds  
    
    #read lower elevation of cells
    dfs_pre = mikeio.read(fpre,items=['Lower level of computational layers in the saturated zone','Thickness of computational layers in the saturated zone'])
    z_bot = dfs_pre['Lower level of computational layers in the saturated zone'].values
    dz = dfs_pre['Thickness of computational layers in the saturated zone'].values
    del dfs_pre
    
    head_above_bottom = heads-z_bot
    dz_sat = head_above_bottom
    for it in range(nt):
        for iz in range(nz):
            for ix in range(nx):
                for iy in range(ny):
                    if dz_sat[it,iz,iy,ix]>dz[0,iz,iy,ix]:
                         dz_sat[it,iz,iy,ix]=dz[0,iz,iy,ix]
                             
    return dz_sat

#:-----------------------------------------------------------------------------
#:-------------get saturated thickness of all cells
def from_Q_to_q(arr1,dz_sat,cross_length,varname):
    flow = arr1.to_numpy() #meter pow 3 per year
    A = cross_length*dz_sat #cross section for flow: meter pow 2
    flux = flow/A #meter per year
    
    x = arr1.x.to_numpy()
    y = arr1.y.to_numpy()
    z = arr1.z.to_numpy()
    time = arr1.time.to_numpy()
    
    arr1.attrs['name'] = varname
    arr1.name = varname
    
    new_arr = xr.DataArray(data=flux,
                           dims=['time','z','y','x'], 
                           coords=dict(time=('time', time),
                                       z=('z', z),
                                       y=('y', y), 
                                       x=('x', x)),
                           name=arr1.name,
                           attrs=dict(name=arr1.attrs['name']))
    del arr1
    return new_arr