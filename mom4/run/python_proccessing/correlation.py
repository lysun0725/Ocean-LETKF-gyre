
from netCDF4 import Dataset
import numpy as np
import readMOM4 as MOM
import os
import matplotlib.pyplot as plt
from datetime import datetime


root_dir = '/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4'
num_ens = 40
fin_date = datetime(1981,1,16,0,0)
date_str = fin_date.strftime('%Y%m%d%H')
letkf_dir = os.path.dirname(root_dir + '/OUTPUT/' + date_str + '/letkf/') 
true_dir  = os.path.dirname(root_dir + '/OBS/' + date_str + '/041/RESTART/')

gfile = os.path.join(letkf_dir,'grid_spec.nc')
lons,lats,levs = MOM.read_netcdf_grids(gfile)
nlon = np.size(lons)
nlat = np.size(lats)
nlev = np.size(levs)


dr_id = 44 
dpos_ens = np.zeros([2,num_ens])
dpos_ens_g = np.zeros([2,num_ens])

# read drifter ensemble data based on the given drifter id
for k in range(1,num_ens+1):
    edrfile = os.path.join(letkf_dir,'anal'+'%03.d'%(k)+'.drifters_inp.nc')
    edrfile_g = os.path.join(letkf_dir,'gs01'+'%03.d'%(k)+'.drifters_inp.txt')
    
    tmp1,tmp2 = MOM.read_netcdf_drifters(edrfile)
    dpos_ens[:,k-1] = tmp2[dr_id-1,0:2]
    tmp1,tmp2 = MOM.read_txt_gs_drifters(edrfile_g)  
    dpos_ens_g[:,k-1] = tmp2[dr_id-1,0:2]

# read true drifter position as mean
tdrfile = os.path.join(true_dir,'drifters_inp.nc')
tmp1,tmp2 = MOM.read_netcdf_drifters(tdrfile)
dpos_mean = tmp2[dr_id-1,0:2]
#dpos_mean = np.mean(dpos_ens,axis=1)
dpos_ens = dpos_ens - np.matrix(dpos_mean).transpose() * np.ones([1,num_ens])

dpos_ens_g = dpos_ens_g - np.matrix(dpos_mean).transpose() * np.ones([1,num_ens])

################# Collect Data in the Vertical Direction #####################
#start searching the indices of the NE grid points
xlon, ylat = MOM.find_a_cell(dpos_mean[0],dpos_mean[1],lons,lats)

#read ocean ensemble data in the vertical direction
temp_v_ens = np.zeros([nlev,num_ens])
salt_v_ens = np.zeros([nlev,num_ens])
uvel_v_ens = np.zeros([nlev,num_ens])
vvel_v_ens = np.zeros([nlev,num_ens])

temp_v_ens_g = np.zeros([nlev,num_ens])
salt_v_ens_g = np.zeros([nlev,num_ens])
uvel_v_ens_g = np.zeros([nlev,num_ens])
vvel_v_ens_g = np.zeros([nlev,num_ens])

for k in range(1,num_ens+1):
    etsfile = os.path.join(letkf_dir,'anal'+'%03.d'%(k)+'.ocean_temp_salt.res.nc')
    euvfile = os.path.join(letkf_dir,'anal'+'%03.d'%(k)+'.ocean_velocity.res.nc')
    etsfile_g = os.path.join(letkf_dir,'gs01'+'%03.d'%(k)+'.ocean_temp_salt.res.nc')
    euvfile_g = os.path.join(letkf_dir,'gs01'+'%03.d'%(k)+'.ocean_velocity.res.nc')

    tmp1,tmp2 = MOM.read_netcdf_ts(etsfile) # tmp1.shape = (nlev,nlat,nlon)
    temp_v_ens[:,k-1] = tmp1[:,ylat,xlon]
    salt_v_ens[:,k-1] = tmp2[:,ylat,xlon]
    
    tmp1,tmp2 = MOM.read_netcdf_uv(euvfile)
    uvel_v_ens[:,k-1] = tmp1[:,ylat,xlon]
    vvel_v_ens[:,k-1] = tmp2[:,ylat,xlon]

    tmp1,tmp2 = MOM.read_netcdf_ts(etsfile_g) # tmp1.shape = (nlev,nlat,nlon)
    temp_v_ens_g[:,k-1] = tmp1[:,ylat,xlon]
    salt_v_ens_g[:,k-1] = tmp2[:,ylat,xlon]
    
    tmp1,tmp2 = MOM.read_netcdf_uv(euvfile_g)
    uvel_v_ens_g[:,k-1] = tmp1[:,ylat,xlon]
    vvel_v_ens_g[:,k-1] = tmp2[:,ylat,xlon]

#read true ocean data as mean
ttsfile = os.path.join(true_dir,'ocean_temp_salt.res.nc')
tuvfile = os.path.join(true_dir,'ocean_velocity.res.nc')

tmp1,tmp2 = MOM.read_netcdf_ts(ttsfile) # tmp1.shape = (nlev,nlat,nlon)2
temp_v_mean = tmp1[:,ylat,xlon]
salt_v_mean = tmp2[:,ylat,xlon]

tmp1,tmp2 = MOM.read_netcdf_uv(tuvfile)
uvel_v_mean = tmp1[:,ylat,xlon]
vvel_v_mean = tmp2[:,ylat,xlon]

temp_v_ens = temp_v_ens - np.matrix(temp_v_mean).transpose() * np.ones([1,num_ens])
salt_v_ens = salt_v_ens - np.matrix(salt_v_mean).transpose() * np.ones([1,num_ens])
uvel_v_ens = uvel_v_ens - np.matrix(uvel_v_mean).transpose() * np.ones([1,num_ens])
vvel_v_ens = vvel_v_ens - np.matrix(vvel_v_mean).transpose() * np.ones([1,num_ens])

temp_v_ens_g = temp_v_ens_g - np.matrix(temp_v_mean).transpose() * np.ones([1,num_ens])
salt_v_ens_g = salt_v_ens_g - np.matrix(salt_v_mean).transpose() * np.ones([1,num_ens])
uvel_v_ens_g = uvel_v_ens_g - np.matrix(uvel_v_mean).transpose() * np.ones([1,num_ens])
vvel_v_ens_g = vvel_v_ens_g - np.matrix(vvel_v_mean).transpose() * np.ones([1,num_ens])

dr_v_temp_corr = np.corrcoef(dpos_ens,temp_v_ens)[0:2,2:(nlev+2)]
dr_v_salt_corr = np.corrcoef(dpos_ens,salt_v_ens)[0:2,2:(nlev+2)]
dr_v_uvel_corr = np.corrcoef(dpos_ens,uvel_v_ens)[0:2,2:(nlev+2)]
dr_v_vvel_corr = np.corrcoef(dpos_ens,vvel_v_ens)[0:2,2:(nlev+2)]

dr_v_temp_corr_g = np.corrcoef(dpos_ens_g,temp_v_ens_g)[0:2,2:(nlev+2)]
dr_v_salt_corr_g = np.corrcoef(dpos_ens_g,salt_v_ens_g)[0:2,2:(nlev+2)]
dr_v_uvel_corr_g = np.corrcoef(dpos_ens_g,uvel_v_ens_g)[0:2,2:(nlev+2)]
dr_v_vvel_corr_g = np.corrcoef(dpos_ens_g,vvel_v_ens_g)[0:2,2:(nlev+2)]

plt.figure(1)
nz=31
plt.subplot(241)
plt.plot(levs[0:(nz)],dr_v_temp_corr[0,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_temp_corr_g[0,0:nz],c='r',marker='*')
#plt.title('Temperature')
#plt.ylabel('Drifter_Lon',rotation=0)
plt.tight_layout()
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.ylim([-1,1])

plt.subplot(242)
plt.plot(levs[0:(nz)],dr_v_salt_corr[0,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_salt_corr_g[0,0:nz],c='r',marker='*')
#plt.title('Salinity')
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.ylim([-1,1])

plt.subplot(243)
plt.plot(levs[0:(nz)],dr_v_uvel_corr[0,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_uvel_corr_g[0,0:nz],c='r',marker='*')
#plt.title('Uvel')
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.ylim([-1,1])

plt.subplot(244)
plt.plot(levs[0:(nz)],dr_v_vvel_corr[0,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_vvel_corr_g[0,0:nz],c='r',marker='*')
#plt.title('Vvel')
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.ylim([-1,1])

plt.subplot(245)
plt.plot(levs[0:(nz)],dr_v_temp_corr[1,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_temp_corr_g[1,0:nz],c='r',marker='*')
#plt.ylabel('Drifter_Lat',rotation=0)
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.tight_layout()
plt.ylim([-1,1])

plt.subplot(246)
plt.plot(levs[0:(nz)],dr_v_salt_corr[1,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_salt_corr_g[1,0:nz],c='r',marker='*')
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.ylim([-1,1])

plt.subplot(247)
plt.plot(levs[0:(nz)],dr_v_uvel_corr[1,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_uvel_corr_g[1,0:nz],c='r',marker='*')
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.ylim([-1,1])

plt.subplot(248)
plt.plot(levs[0:(nz)],dr_v_vvel_corr[1,0:nz],c='b',marker='*')
plt.plot(levs[0:(nz)],dr_v_vvel_corr_g[1,0:nz],c='r',marker='*')
plt.xlabel('Depth (m)')
plt.ylabel('Correlation')
plt.xscale('log')
plt.axhline(y=0.0, color='black', linestyle='--')
plt.axvline(x=15.0, color='black', linestyle='--')
plt.ylim([-1,1])

plt.show()



