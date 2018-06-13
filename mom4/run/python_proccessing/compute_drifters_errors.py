
# coding: utf-8

# In[35]:


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import readMOM4 as MOM
from geopy.distance import great_circle
import os
from datetime import datetime, timedelta

def txt_netcdf(txtfile,ncfile):
    txt_ids,txt_pos = MOM.read_txt_drifters(txtfile)
    nc_ids,nc_pos = MOM.read_netcdf_drifters(ncfile)

    dist = np.zeros(txt_ids.shape)
    for i in range(0,txt_ids.shape[0]):
        for j in range(0,nc_ids.shape[0]):
            if txt_ids[i] == nc_ids[j]:
                dist[i] = great_circle(txt_pos[i,0:2],nc_pos[j,0:2]).km
                
            
    return txt_ids, dist

def txt_gstxt(txtfile1,txtfile2):
    txt1_ids,txt1_pos = MOM.read_txt_drifters(txtfile1)
    txt2_ids,txt2_pos = MOM.read_txt_gs_drifters(txtfile2)
   
    dist = np.zeros(txt1_ids.shape)
    for i in range(0,txt1_ids.shape[0]):
        for j in range(0,txt2_ids.shape[0]):
            if txt1_ids[i] == txt2_ids[j]:
                dist[i] = great_circle(txt1_pos[i,0:2],txt2_pos[j,0:2]).km    
   
    return txt1_ids, dist

def netcdf_netcdf(ncfile1,ncfile2):
    nc1_ids, nc1_pos = MOM.read_netcdf_drifters(ncfile1)
    nc2_ids, nc2_pos = MOM.read_netcdf_drifters(ncfile2)

    dist = np.zeros(nc1_ids.shape)
    for i in range(0,nc1_ids.shape[0]):
        for j in range(0,nc2_ids.shape[0]):
            if nc1_ids[i] == nc2_ids[j]:
                dist[i] = great_circle(nc1_pos[i,0:2],nc2_pos[i,0:2]).km

    return nc1_ids,dist

root_dir = '/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4'
days = 16 
num_dr = 50
num_ens = 40
iday = 1
tma_dist = np.zeros((num_dr,days))
tma2_dist = np.zeros((num_dr,days))
tma3_dist = np.zeros((num_dr,days))
tmg_dist = np.zeros((num_dr,days))
tmc_dist = np.zeros((num_dr,days))
spa_dist = np.zeros((num_dr,days))
spg_dist = np.zeros((num_dr,days))
ini_date = datetime(1981,1,15,0,0)
time = np.arange(1,days+1,iday)
tpos = np.zeros((num_dr,3,days))
gpos = np.zeros((num_dr,3,days))
apos = np.zeros((num_dr,3,days))
cpos = np.zeros((num_dr,3,days))

f = open(os.path.join(os.path.dirname(root_dir + '/OUTPUT/EXP5/'),"drif_rmse.txt"),"w")
f.write("TMG        TMA        SPG        SPA        CON\n")

for i in range(1,days+1):
  date = ini_date + timedelta(days=i)
  date_str = date.strftime('%Y%m%d%H')

  letkf_dir = os.path.dirname(root_dir + '/OUTPUT/EXP5/' + date_str + '/letkf/')
  #letkf2_dir = os.path.dirname(root_dir + '/OUTPUT/EXP2/' + date_str + '/letkf/') 
  #letkf3_dir = os.path.dirname(root_dir + '/OUTPUT/EXP3/' + date_str + '/letkf/')
  true_dir = os.path.dirname(root_dir + '/OBS2/' + date_str + '/041/RESTART/')
  ctrl_dir = os.path.dirname(root_dir + '/CONTROL/' + date_str + '/CON/RESTART/')

  anal = os.path.join(letkf_dir,'andr_me.txt') 
  #anal2 = os.path.join(letkf2_dir,'gsdr_me.txt')
  #anal3 = os.path.join(letkf3_dir,'gsdr_me.txt')
  gues = os.path.join(letkf_dir,'gsdr_me.txt')
  true = os.path.join(true_dir,'drifters_inp.nc')
  ctrl = os.path.join(ctrl_dir,'drifters_inp.nc')

  tma_ids,tma_dist[:,i-1] = txt_netcdf(anal,true)
  #tma_ids,tma2_dist[:,i-1] = txt_netcdf(anal2,true)
  #tma_ids,tma3_dist[:,i-1] = txt_netcdf(anal3,true)
  tmg_ids,tmg_dist[:,i-1] = txt_netcdf(gues,true)
  tmc_ids,tmc_dist[:,i-1] = netcdf_netcdf(ctrl,true)

  tmp_id = np.zeros((num_dr,1))
  tmp_an_sp = np.zeros((num_dr,1))
  tmp_gs_sp = np.zeros((num_dr,1))

  # code the andr_sp and gsdr_sp later 
  for j in range(1,num_ens+1):
    an_ens = os.path.join(letkf_dir,'anal' + '{:03d}'.format(j) + '.drifters_inp.nc')
    gs_ens = os.path.join(letkf_dir,'gs'+ '{:02d}'.format(iday)  + '{:03d}'.format(j) + '.drifters_inp.txt')   

    tmp_id,tmp_an_dist = txt_netcdf(anal,an_ens) 
    tmp_id,tmp_gs_dist = txt_gstxt(gues,gs_ens)

    tmp_an_sp[:,0] = tmp_an_sp[:,0] + np.square(tmp_an_dist[:])
    tmp_gs_sp[:,0] = tmp_gs_sp[:,0] + np.square(tmp_gs_dist[:])

  spa_dist[:,i-1] = np.sqrt(tmp_an_sp[:,0]/(num_ens-1))
  spg_dist[:,i-1] = np.sqrt(tmp_gs_sp[:,0]/(num_ens-1))

  # code for plotting drifter trajectories: true, gues_me and anal_me
  tid, tpos[:,:,i-1] = MOM.read_netcdf_drifters(true)
  gid, gpos[:,:,i-1] = MOM.read_txt_drifters(gues)
  aid, apos[:,:,i-1] = MOM.read_txt_drifters(anal)
  cid, cpos[:,:,i-1] = MOM.read_netcdf_drifters(ctrl)

  f.write(str(np.mean(tmg_dist[:,i-1])) + " " +  str(np.mean(tma_dist[:,i-1])) + " " + str(np.mean(spg_dist[:,i-1])) + " " + str(np.mean(spa_dist[:,i-1])) + " " + str(np.mean(tmc_dist[:,i-1])) + "\n")

f.close()

  

# plotting the average drifters distance error in the testing time interval
plt.figure(1)
plt.plot(time,np.mean(tmg_dist[:,:],axis=0),label='true - forecast',c='b')
plt.plot(time,np.mean(tma_dist[:,:],axis=0),label='true - analysis',c='r')
plt.plot(time,np.mean(tmc_dist[:,:],axis=0),label='true - control', c='g')
plt.plot(time,np.mean(spg_dist[:,:],axis=0),label='forecast spread',linestyle='--',c='b')
plt.plot(time,np.mean(spa_dist[:,:],axis=0),label='analysis spread',linestyle='--',c='r')
plt.yscale('log')
plt.xlabel('time (days)')
plt.ylabel('log distance (km)')
plt.title('Drifter Errors in Distance (km)')
plt.legend()
#plt.savefig('RMSE_dr.png')

# plotting the trajectories of a certain drifter
ind = 30 
plt.figure(2)
plt.plot(gpos[ind-1,0,:],gpos[ind-1,1,:],linestyle='--',marker='s',label='forecast traj',c='b')
plt.plot(apos[ind-1,0,:],apos[ind-1,1,:],linestyle='--',marker='v',label='analysis traj',c='r')
plt.plot(cpos[ind-1,0,:],cpos[ind-1,1,:],linestyle='--',marker='p',label='control traj',c='g')
plt.plot(tpos[ind-1,0,:],tpos[ind-1,1,:],linestyle='--',marker='*',label='true traj',c='k')
plt.xlabel('lon')
plt.ylabel('lat')
plt.title('Drifter Trajectories Comparison No.'+'{:02d}'.format(ind))
plt.legend()

ind = 47
plt.figure(3)
plt.plot(gpos[ind-1,0,:],gpos[ind-1,1,:],linestyle='--',marker='s',label='forecast traj',c='b')
plt.plot(apos[ind-1,0,:],apos[ind-1,1,:],linestyle='--',marker='v',label='analysis traj',c='r')
plt.plot(cpos[ind-1,0,:],cpos[ind-1,1,:],linestyle='--',marker='p',label='control traj',c='g')
plt.plot(tpos[ind-1,0,:],tpos[ind-1,1,:],linestyle='--',marker='*',label='true traj',c='k')
plt.xlabel('lon')
plt.ylabel('lat')
plt.title('Drifter Trajectories Comparison No.'+'{:02d}'.format(ind))
plt.legend()

plt.figure(4)
std_dist=great_circle([5,25],[5.1,25.1]).km

plt.plot(time,np.mean(tmc_dist[:,:],axis=0)/std_dist,label='Control Run',c='b')
#plt.plot(time,np.mean(tma3_dist[:,:],axis=0)/std_dist,label='EXP1',c='orange',linestyle='--')
plt.plot(time,np.mean(tma_dist[:,:],axis=0)/std_dist,label='EXP2', c='r',linestyle='--')
#plt.plot(time,np.mean(tma2_dist[:,:],axis=0)/std_dist,label='EXP3',linestyle='--',c='m')

plt.yscale('log')
plt.xlabel('time (days)')
plt.ylabel('|x_D|')
plt.title('Drifter Norm')
plt.legend()
plt.savefig('NORM_dr_gs.png')
#plt.show()

#print 'tma_dist = ', np.mean(tma_dist[:,:],axis=0)
#print 'tmg_dist = ', np.mean(tmg_dist[:,:],axis=0)
#print 'spa_dist = ', np.mean(spa_dist[:,:],axis=0)
#print 'spg_dist = ', np.mean(spg_dist[:,:],axis=0)




