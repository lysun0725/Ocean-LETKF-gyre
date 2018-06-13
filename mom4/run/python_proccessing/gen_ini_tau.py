
# coding: utf-8

# In[18]:


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

tau = Dataset('tau.nc',mode='r')

#lons = tau.variables['grid_x_C'][:]
#lats = tau.variables['grid_y_C'][:]
taux = tau.variables['taux'][0,:,:]
tauy = tau.variables['tauy'][0,:,:]

# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
#lon, lat = np.meshgrid(lons, lats)

# Plot Data
#plt.quiver(lon,lat,taux,tauy)

#plt.show()


# In[21]:


# initialize ensemble members
num_ens = 80
x_rand = np.random.rand(num_ens,1)*.2-.1
y_rand = np.random.rand(num_ens,1)*.2-.1
toexclude = ['taux','tauy']

for i in range(0,num_ens):
    tau_ens = Dataset('tau.%0.3d.nc' % (i+1),'w',format=tau.file_format)

    # create dimension
    for dimname in tau.dimensions.keys():
        dim = tau_ens.createDimension(dimname,tau.dimensions[dimname].size)
    
    # create variables
    for varname in tau.variables.keys():
        var = tau_ens.createVariable(varname,tau.variables[varname].dtype,tau.variables[varname].dimensions)
    
        # adding attributes
        for attrname in tau.variables[varname].ncattrs():
            var.setncattr(attrname,getattr(tau.variables[varname],attrname))
    
        # Passing variables value
        if varname not in toexclude:
            tau_ens.variables[varname][:] = tau.variables[varname][:]
        else:
            tau_ens.variables[varname][0,:,:] = tau.variables[varname][0,:,:]+(np.random.rand(1,1)*.2-.1)*np.ones(tau.variables[varname][0,:,:].shape)    
    # create global attributes
    for gattrname in tau.ncattrs():
        tau_ens.setncattr(gattrname,getattr(tau,gattrname))
    
    tau_ens.close()
    
tau.close()
    


# In[8]:





# In[ ]:




