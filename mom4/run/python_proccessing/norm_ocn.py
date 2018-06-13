import numpy as np
import matplotlib.pyplot as plt
import os

output_dir = os.path.dirname('/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4/OUTPUT/')
num_days = 105
ncol = 4 

label1 = 'DA w/ only T,S; 3L_R'
label2 = 'LaDA w/o T,S; 3L_R'
label3 = 'LaDA w/ T,S; 3L_R'

temp = np.zeros([num_days,ncol])
salt = np.zeros([num_days,ncol])
kine = np.zeros([num_days,ncol])
da_cycle = np.arange(1,num_days+1,1)

tfile = os.path.join(output_dir,'temp_norm_3d_1.txt')
sfile = os.path.join(output_dir,'salt_norm_3d_1.txt')
kfile = os.path.join(output_dir,'KE_norm_3d_1.txt')

temp = np.loadtxt(tfile,usecols=range(0,ncol),skiprows=1)
salt = np.loadtxt(sfile,usecols=range(0,ncol),skiprows=1)
kine = np.loadtxt(kfile,usecols=range(0,ncol),skiprows=1)

plt.figure(1)
plt.plot(da_cycle,temp[:,3]*100,label='Control Run',c='b')
plt.plot(da_cycle,temp[:,1]*100,label=label1,c='orange',linestyle='--')
plt.plot(da_cycle,temp[:,0]*100,label=label2,c='r',linestyle='--')
plt.plot(da_cycle,temp[:,2]*100,label=label3,c='m',linestyle='--')
#plt.plot(da_cycle,temp[:,2]*100,label='w/o T,S; 2L_R',linestyle='--',c='m')

plt.xlabel('time (days)')
plt.ylabel('|T| %')
plt.title('Temperature')
plt.legend()
plt.tight_layout()
plt.savefig('NORM_percent_t_ts.png')

plt.figure(2)
plt.plot(da_cycle,salt[:,3]*100,label='Control Run',c='b')
plt.plot(da_cycle,salt[:,1]*100,label=label1,c='orange',linestyle='--')
plt.plot(da_cycle,salt[:,0]*100,label=label2,c='r',linestyle='--')
plt.plot(da_cycle,salt[:,2]*100,label=label3,c='m',linestyle='--')
#plt.plot(da_cycle,salt[:,2]*100,label='w/o T,S; 2L_R',linestyle='--',c='m')

plt.xlabel('time (days)')
plt.ylabel('|S| %')
plt.title('Salinity')
#plt.yscale('log')
plt.legend()
plt.savefig('NORM_percent_s_ts.png')

plt.figure(3)
plt.plot(da_cycle,kine[:,3]*100,label='Control Run',c='b')
plt.plot(da_cycle,kine[:,1]*100,label=label1,c='orange',linestyle='--')
plt.plot(da_cycle,kine[:,0]*100,label=label2,c='r',linestyle='--')
plt.plot(da_cycle,kine[:,2]*100,label=label3,c='m',linestyle='--')
#plt.plot(da_cycle,kine[:,2]*100,label='w/o T,S; 2_R',linestyle='--',c='m')

plt.xlabel('time (days)')
plt.ylabel('|KE| %')
plt.title('KE')
plt.legend()
plt.savefig('NORM_percent_KE_ts.png')

plt.show()

