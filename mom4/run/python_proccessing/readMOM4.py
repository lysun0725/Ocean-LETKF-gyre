
# coding: utf-8

# In[2]:


from netCDF4 import Dataset
import numpy as np

def read_netcdf_drifters(dfile):
    ocn_dr_infile = Dataset(dfile,mode='r')
    dr_ids = ocn_dr_infile.variables['ids'][:]
    dr_pos = ocn_dr_infile.variables['positions'][:,0:3]
    return dr_ids, dr_pos

# This function is for reading the andr_me, gsdr_me
def read_txt_drifters(dfile):
    ocn_dr_infile = np.loadtxt(dfile,usecols=range(0,4)) #num_dr -by- 4
    dr_ids = ocn_dr_infile[:,0]
    dr_pos = ocn_dr_infile[:,1:4]
    return dr_ids, dr_pos

# This function is for reading the ensemble drifter files.
def read_txt_gs_drifters(dfile):
    ocn_dr_infile = np.loadtxt(dfile,usecols=range(0,8),skiprows=2)
    n = np.shape(ocn_dr_infile)[0]
    
    dr_ids = np.zeros((n,))
    dr_pos = np.zeros((n,3))
    counter = 0
    for i in range(0,n):
        if (ocn_dr_infile[i,6] == 24.):
            counter += 1
            dr_ids[counter] = ocn_dr_infile[i,0]
            dr_pos[counter,0:2] = ocn_dr_infile[i,1:3]
    dr_ids = dr_ids[np.nonzero(dr_ids)]
    dr_pos = dr_pos[~np.all(dr_pos==0,axis=1)]
    return dr_ids, dr_pos

def read_netcdf_grids(gfile):
    ocn_grid_file = Dataset(gfile,mode='r') 
    ocn_lons = ocn_grid_file.variables['grid_x_T'][:]
    ocn_lats = ocn_grid_file.variables['grid_y_T'][:]
    ocn_levs = ocn_grid_file.variables['zt'][:]
    return ocn_lons, ocn_lats, ocn_levs

def read_netcdf_ts(tsfile):
    ocn_ts_file = Dataset(tsfile,mode='r')
    ocn_temp = ocn_ts_file.variables['temp'][0,:,:,:]
    ocn_salt = ocn_ts_file.variables['salt'][0,:,:,:]
    return ocn_temp, ocn_salt

def read_netcdf_uv(uvfile):
    ocn_uv_file = Dataset(uvfile,mode='r')
    ocn_uvel = ocn_uv_file.variables['u'][0,:,:,:,]
    ocn_vvel = ocn_uv_file.variables['v'][0,:,:,:,]
    return ocn_uvel,ocn_vvel

def read_netcdf_ss(ssfile,varname):
    ocn_ss_file = Dataset(ssfile,mode='r')
    ocn_ss_var = ocn_ss_file.variables[varname][0,:,:]
    return ocn_ss_var

def read_netcdf_3d(infile,varname):
    ocn_infile = Dataset(infile,mode='r')
    ocn_var = ocn_infile.variables[varname][0,:,:,:]
    return ocn_var

# find the NE corner of the cell
def find_a_cell(x,y,ocn_lons,ocn_lats):
    ocn_nlon = np.size(ocn_lons)
    ocn_nlat = np.size(ocn_lats)
    for i in range(1,ocn_nlon):
   	 if x>=ocn_lons[i-1] and x<ocn_lons[i]:
		break

    for j in range(1,ocn_nlat):
    	if y>=ocn_lats[j-1] and y<ocn_lats[j]:
		break
    return i,j 

# In[ ]:




