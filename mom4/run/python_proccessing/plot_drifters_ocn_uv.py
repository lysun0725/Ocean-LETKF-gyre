
# coding: utf-8

# In[40]:


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import readMOM4 as MOM
import os

#os.system('cd ../../OUTPUT/EXP2/1981043000/letkf/')
#os.system('rm -rf anme.*.nc')
#os.system('ncea anal0??.ocean_velocity.res.nc anme.ocean_velocity.res.nc')


file_dir = os.path.dirname('/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4/OUTPUT/1981010100/letkf/')
true_dir = os.path.dirname('/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4/OBS/1981010100/041/RESTART/')

ocn_grid_file = os.path.join(file_dir,'grid_spec.nc')
ocn_uv_infile = os.path.join(true_dir,'ocean_velocity.res.nc')# ocean velocity file
#ocn_dr_infile = 'drifters_inp.nc'
#ocn_dr_anfile = os.path.join(file_dir,'andr_me.txt')
ocn_dr_anfile = os.path.join(file_dir,'andr_me.nc')

ocn_lons,ocn_lats,ocn_levs = MOM.read_netcdf_grids(ocn_grid_file)
lon, lat = np.meshgrid(ocn_lons, ocn_lats)
ocn_uvel,ocn_vvel = MOM.read_netcdf_uv(ocn_uv_infile)

#dr_ids,dr_pos = MOM.read_txt_drifters(ocn_dr_anfile)
dr_ids,dr_pos = MOM.read_netcdf_drifters(ocn_dr_anfile)

plt.quiver(lon,lat,ocn_uvel[1,:,:],ocn_vvel[1,:,:])
plt.scatter(dr_pos[:,0],dr_pos[:,1])

# label the scatters
labels = ['{:2}'.format(int(id)) for id in dr_ids]
i=0
for xy in zip(dr_pos[:,0],dr_pos[:,1]):
    plt.annotate(labels[i],xy)
    i += 1
    
plt.grid()
plt.rcParams["figure.figsize"] = [16.0,32.0]
plt.axis('equal')
plt.show()







