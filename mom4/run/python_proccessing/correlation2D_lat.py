
from netCDF4 import Dataset
import numpy as np
import readMOM4 as MOM
import os
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.transforms import offset_copy


root_dir = '/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4'
num_ens = 40
fin_date = datetime(1981,1,16,0,0)
date_str = fin_date.strftime('%Y%m%d%H')
letkf_dir = os.path.dirname(root_dir + '/OUTPUT/ORIG_PERTUB_DR/' + date_str + '/letkf/') 
true_dir  = os.path.dirname(root_dir + '/OBS/' + date_str + '/041/RESTART/')

gfile = os.path.join(letkf_dir,'grid_spec.nc')
lons,lats,levs = MOM.read_netcdf_grids(gfile)
nlon = np.size(lons)
nlat = np.size(lats)
nlev = np.size(levs)


dr_id = 11
dpos_ens_g = np.zeros([2,num_ens])

# read drifter ensemble data based on the given drifter id
for k in range(1,num_ens+1):
    edrfile_g = os.path.join(letkf_dir,'gs01'+'%03.d'%(k)+'.drifters_inp.txt')

    tmp1,tmp2 = MOM.read_txt_gs_drifters(edrfile_g)  
    dpos_ens_g[:,k-1] = tmp2[dr_id-1,0:2]

# read true drifter position as mean
tdrfile = os.path.join(true_dir,'drifters_inp.nc')
tmp1,tmp2 = MOM.read_netcdf_drifters(tdrfile)
dpos_mean = tmp2[dr_id-1,0:2]
dpos_ens_g = dpos_ens_g - np.matrix(dpos_mean).transpose() * np.ones([1,num_ens])

################# Collect Data in the Vertical and Longitude Direction #####################
#start searching the indices of the NE grid points
xlon, ylat = MOM.find_a_cell(dpos_mean[0],dpos_mean[1],lons,lats)

temp_v_ens_g = np.zeros([nlat,nlev,num_ens])
salt_v_ens_g = np.zeros([nlat,nlev,num_ens])
uvel_v_ens_g = np.zeros([nlat,nlev,num_ens])
vvel_v_ens_g = np.zeros([nlat,nlev,num_ens])
ssh_v_ens_g = np.zeros([1,num_ens])

for k in range(1,num_ens+1):
    etsfile_g = os.path.join(letkf_dir,'gs01' + '%03.d'%(k) + '.ocean_temp_salt.res.nc')
    euvfile_g = os.path.join(letkf_dir,'gs01' + '%03.d'%(k) + '.ocean_velocity.res.nc')

    # read ssh file
    emodel_dir = os.path.dirname(root_dir + '/OUTPUT/ORIG_PERTUB_DR/' + date_str + '/model/' + '%03.d'%(k) + '/RESTART/') 
    essfile_g = os.path.join(emodel_dir,'ocean_sbc.res.nc')

    tmp1,tmp2 = MOM.read_netcdf_ts(etsfile_g) # tmp1.shape = (nlev,nlat,nlon)
    temp_v_ens_g[:,:,k-1] = tmp1[:,0:nlat,xlon].transpose()
    salt_v_ens_g[:,:,k-1] = tmp2[:,0:nlat,xlon].transpose()
    
    tmp1,tmp2 = MOM.read_netcdf_uv(euvfile_g)
    uvel_v_ens_g[:,:,k-1] = tmp1[:,0:nlat,xlon].transpose()
    vvel_v_ens_g[:,:,k-1] = tmp2[:,0:nlat,xlon].transpose()

    tmp3 = MOM.read_netcdf_ss(essfile_g,'sea_lev') # tmp3.shape=(nlat,nlon)
    ssh_v_ens_g[:,k-1] = tmp3[ylat,xlon]

#read true ocean data as mean
ttsfile = os.path.join(true_dir,'ocean_temp_salt.res.nc')
tuvfile = os.path.join(true_dir,'ocean_velocity.res.nc')
tssfile = os.path.join(true_dir,'ocean_sbc.res.nc')

tmp1,tmp2 = MOM.read_netcdf_ts(ttsfile) # tmp1.shape = (nlev,nlat,nlon)2
temp_v_mean = tmp1[:,0:nlat,xlon].transpose()
salt_v_mean = tmp2[:,0:nlat,xlon].transpose()

tmp1,tmp2 = MOM.read_netcdf_uv(tuvfile)
uvel_v_mean = tmp1[:,0:nlat,xlon].transpose()
vvel_v_mean = tmp2[:,0:nlat,xlon].transpose()

tmp3 = MOM.read_netcdf_ss(tssfile,'sea_lev')
ssh_v_mean = tmp3[ylat,xlon]

for k in range(1,num_ens+1):
    temp_v_ens_g[:,:,k-1] = temp_v_ens_g[:,:,k-1] - temp_v_mean
    salt_v_ens_g[:,:,k-1] = salt_v_ens_g[:,:,k-1] - salt_v_mean
    uvel_v_ens_g[:,:,k-1] = uvel_v_ens_g[:,:,k-1] - uvel_v_mean
    vvel_v_ens_g[:,:,k-1] = vvel_v_ens_g[:,:,k-1] - vvel_v_mean
    ssh_v_ens_g[:,k-1] = ssh_v_ens_g[:,k-1] - ssh_v_mean

dr_v_temp_corr_g = np.zeros([2,nlat,nlev])
dr_v_salt_corr_g = np.zeros([2,nlat,nlev])
dr_v_uvel_corr_g = np.zeros([2,nlat,nlev])
dr_v_vvel_corr_g = np.zeros([2,nlat,nlev])

ssh_v_temp_corr_g = np.zeros([nlat,nlev])
ssh_v_salt_corr_g = np.zeros([nlat,nlev])
ssh_v_uvel_corr_g = np.zeros([nlat,nlev])
ssh_v_vvel_corr_g = np.zeros([nlat,nlev])


for i in range(1,nlat+1):
    dr_v_temp_corr_g[:,i-1,:] = np.corrcoef(dpos_ens_g,temp_v_ens_g[i-1,:,:])[0:2,2:(nlev+2)]
    dr_v_salt_corr_g[:,i-1,:] = np.corrcoef(dpos_ens_g,salt_v_ens_g[i-1,:,:])[0:2,2:(nlev+2)]
    dr_v_uvel_corr_g[:,i-1,:] = np.corrcoef(dpos_ens_g,uvel_v_ens_g[i-1,:,:])[0:2,2:(nlev+2)]
    dr_v_vvel_corr_g[:,i-1,:] = np.corrcoef(dpos_ens_g,vvel_v_ens_g[i-1,:,:])[0:2,2:(nlev+2)]

    ssh_v_temp_corr_g[i-1,:] = np.corrcoef(ssh_v_ens_g,temp_v_ens_g[i-1,:,:])[0:1,1:(nlev+1)]
    ssh_v_salt_corr_g[i-1,:] = np.corrcoef(ssh_v_ens_g,salt_v_ens_g[i-1,:,:])[0:1,1:(nlev+1)]
    ssh_v_uvel_corr_g[i-1,:] = np.corrcoef(ssh_v_ens_g,uvel_v_ens_g[i-1,:,:])[0:1,1:(nlev+1)]
    ssh_v_vvel_corr_g[i-1,:] = np.corrcoef(ssh_v_ens_g,vvel_v_ens_g[i-1,:,:])[0:1,1:(nlev+1)]

plt.figure(figsize=(14, 8))
nz=nlev
X,Y = np.meshgrid(lats,levs[0:(nz)])

ax=plt.subplot(3,4,1)
plt.pcolor(X,Y,dr_v_temp_corr_g[0,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.title('Temp')
ax.annotate('dr_lon',xy=(0, 0.5), xytext=(-ax.yaxis.labelpad-5,0),xycoords=ax.yaxis.label, textcoords='offset points',size='large', ha='right', va='center')
plt.colorbar()

plt.subplot(3,4,2)
plt.pcolor(X,Y,dr_v_salt_corr_g[0,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.title('Salt')
plt.colorbar()

plt.subplot(3,4,3)
plt.pcolor(X,Y,dr_v_uvel_corr_g[0,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.title('Uvel')
plt.colorbar()

plt.subplot(3,4,4)
plt.pcolor(X,Y,dr_v_vvel_corr_g[0,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.title('Vvel')
plt.colorbar()

ax=plt.subplot(3,4,5)
plt.pcolor(X,Y,dr_v_temp_corr_g[1,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
ax.annotate('dr_lat',xy=(0, 0.5), xytext=(-ax.yaxis.labelpad-5,0),xycoords=ax.yaxis.label, textcoords='offset points',size='large', ha='right', va='center')
plt.colorbar()

plt.subplot(3,4,6)
plt.pcolor(X,Y,dr_v_salt_corr_g[1,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.colorbar()

plt.subplot(3,4,7)
plt.pcolor(X,Y,dr_v_uvel_corr_g[1,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.colorbar()

plt.subplot(3,4,8)
plt.pcolor(X,Y,dr_v_vvel_corr_g[1,:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.colorbar()

ax=plt.subplot(3,4,9)
plt.pcolor(X,Y,ssh_v_temp_corr_g[:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
ax.annotate('SSH',xy=(0, 0.5), xytext=(-ax.yaxis.labelpad-5,0),xycoords=ax.yaxis.label, textcoords='offset points',size='large', ha='right', va='center')
plt.colorbar()

plt.subplot(3,4,10)
plt.pcolor(X,Y,ssh_v_salt_corr_g[:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.colorbar()

plt.subplot(3,4,11)
plt.pcolor(X,Y,ssh_v_uvel_corr_g[:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.colorbar()

plt.subplot(3,4,12)
plt.pcolor(X,Y,ssh_v_vvel_corr_g[:,0:nz].transpose(),cmap='RdBu',vmin=-1, vmax=1)
plt.tight_layout()
plt.xlabel('Lat (degree)')
plt.ylabel('Depth (m)')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.axhline(y=15.0, color='black', linestyle='--')
plt.axvline(x=lats[ylat], color='black', linestyle='--')
plt.colorbar()


plt.show()



