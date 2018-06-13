import numpy as np
import matplotlib.pyplot as plt
import os

zlev=np.array([5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 235.5619, 250.0461, 269.52, 297.2319, 334.083, 383.1031, 444.9251, 522.265, 615.3997])

output_dir=os.path.dirname('/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4/OUTPUT/')

label1 = 'DA w/ only T,S; 3L_R'
label2 = 'LaDA w/o T,S; 3L_R'
label3 = 'LaDA w/ T,S; 3L_R'

tfile = os.path.join(output_dir,'temp_norm_3d_z.txt')
sfile = os.path.join(output_dir,'salt_norm_3d_z.txt')
kefile = os.path.join(output_dir,'KE_norm_3d_z.txt')

temp = np.loadtxt(tfile,usecols=range(0,4),skiprows=1)
salt = np.loadtxt(sfile,usecols=range(0,4),skiprows=1)
ke = np.loadtxt(kefile,usecols=range(0,4),skiprows=1)

plt.figure(1,figsize=(14,8))

plt.subplot(1,3,1)
plt.plot(temp[:,3]*100,zlev,label='Control')
plt.plot(temp[:,1]*100,zlev,label=label1,c='orange',linestyle='--')
plt.plot(temp[:,0]*100,zlev,label=label2,c='r',linestyle='--')
plt.plot(temp[:,2]*100,zlev,label=label3,c='m',linestyle='--')
plt.legend()
plt.title('Temperature')
plt.xlabel('|T|%')
plt.ylabel('depth')
plt.tight_layout()
plt.xlim([0,1])
plt.gca().invert_yaxis()

plt.subplot(1,3,2)
plt.plot(salt[:,3]*100,zlev,label='Control')
plt.plot(salt[:,1]*100,zlev,label=label1,c='orange',linestyle='--')
plt.plot(salt[:,0]*100,zlev,label=label2,c='r',linestyle='--')
plt.plot(salt[:,2]*100,zlev,label=label3,c='m',linestyle='--')
plt.legend()
plt.title('Salinity')
plt.xlabel('|S|%')
plt.ylabel('depth')
plt.tight_layout()
plt.gca().invert_yaxis()

plt.subplot(1,3,3)
plt.plot(ke[:,3]*100,zlev,label='Control')
plt.plot(ke[:,1]*100,zlev,label=label1,c='orange',linestyle='--')
plt.plot(ke[:,0]*100,zlev,label=label2,c='r',linestyle='--')
plt.plot(ke[:,2]*100,zlev,label=label3,c='m',linestyle='--')
plt.legend()
plt.title('KE')
plt.xlabel('|KE|%')
plt.ylabel('depth')
plt.tight_layout()
plt.xlim([0,30])
plt.gca().invert_yaxis()

plt.savefig('NORM_z_105_ts.png')

plt.show()
