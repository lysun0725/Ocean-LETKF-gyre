*
* This script is for computing the rmse and output the corresponding data as a txt file
*
'open gues_me.ctl'
'set lev 15'
'open anal_me.ctl'
'open ../OBS/ocean_temp_salt_1.res.ctl'
'open ../OBS/ocean_velocity_1.res.ctl'
'open gues_sp.ctl'
'open anal_sp.ctl'
'open ../CONTROL/ocean_temp_salt_1.res.ctl'
'open ../CONTROL/ocean_velocity_1.res.ctl'
'set gxout print'
'set prnopts %1.6e 1 1'    ;*%6.2f: c format, 1: values to plot on each line, ;space between values
tmax = 75

********** TEMPERATURE ***********
write('temp.txt','TMG          TMA          SPG          SPA          CON')
t=1
while(t<=tmax)
  'set t 't 
  'define RMSE=sqrt(aave((t.1-t.3)*(t.1-t.3),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE'
  tmgline = sublin(result,2)
  tmg = subwrd(tmgline,1)

  'define RMSE2=sqrt(aave((t.2-t.3)*(t.2-t.3),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE2'
  tmaline = sublin(result,2)
  tma = subwrd(tmaline,1)
 
  'define RMSE3=sqrt(aave(t.5*t.5,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE3'
  spgline = sublin(result,2)
  spg = subwrd(spgline,1)

  'define RMSE4=sqrt(aave(t.6*t.6,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE4'
  spaline = sublin(result,2)
  spa = subwrd(spaline,1)

  'define RMSE5=sqrt(aave((t.7-t.3)*(t.7-t.3),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE5'
  conline = sublin(result,2)
  con = subwrd(conline,1)

  write('temp.txt', tmg' 'tma' 'spg' 'spa' 'con,append )
  t=t+1
endwhile

********** SALINITY ***********
write('salt.txt','TMG          TMA          SPG          SPA          CON')
t=1
while(t<=tmax)
  'set t 't 
  'define RMSE=sqrt(aave((s.1-s.3)*(s.1-s.3),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE'
  tmgline = sublin(result,2)
  tmg = subwrd(tmgline,1)

  'define RMSE2=sqrt(aave((s.2-s.3)*(s.2-s.3),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE2'
  tmaline = sublin(result,2)
  tma = subwrd(tmaline,1)
 
  'define RMSE3=sqrt(aave(s.5*s.5,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE3'
  spgline = sublin(result,2)
  spg = subwrd(spgline,1)

  'define RMSE4=sqrt(aave(s.6*s.6,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE4'
  spaline = sublin(result,2)
  spa = subwrd(spaline,1)

  'define RMSE5=sqrt(aave((s.7-s.3)*(s.7-s.3),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE5'
  conline = sublin(result,2)
  con = subwrd(conline,1)

  write('salt.txt', tmg' 'tma' 'spg' 'spa' 'con,append )
  t=t+1
endwhile

********** U ***********
write('uvel.txt','TMG          TMA          SPG          SPA          CON')
t=1
while(t<=tmax)
  'set t 't 
  'define RMSE=sqrt(aave((u.1-u.4)*(u.1-u.4),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE'
  tmgline = sublin(result,2)
  tmg = subwrd(tmgline,1)

  'define RMSE2=sqrt(aave((u.2-u.4)*(u.2-u.4),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE2'
  tmaline = sublin(result,2)
  tma = subwrd(tmaline,1)
 
  'define RMSE3=sqrt(aave(u.5*u.5,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE3'
  spgline = sublin(result,2)
  spg = subwrd(spgline,1)

  'define RMSE4=sqrt(aave(u.6*u.6,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE4'
  spaline = sublin(result,2)
  spa = subwrd(spaline,1)

  'define RMSE5=sqrt(aave((u.8-u.4)*(u.8-u.4),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE5'
  conline = sublin(result,2)
  con = subwrd(conline,1)

  write('uvel.txt', tmg' 'tma' 'spg' 'spa' 'con,append)
  t=t+1
endwhile

********** V ***********
write('vvel.txt','TMG          TMA          SPG          SPA          CON')
t=1
while(t<=tmax)
  'set t 't 
  'define RMSE=sqrt(aave((v.1-v.4)*(v.1-v.4),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE'
  tmgline = sublin(result,2)
  tmg = subwrd(tmgline,1)

  'define RMSE2=sqrt(aave((v.2-v.4)*(v.2-v.4),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE2'
  tmaline = sublin(result,2)
  tma = subwrd(tmaline,1)
 
  'define RMSE3=sqrt(aave(v.5*v.5,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE3'
  spgline = sublin(result,2)
  spg = subwrd(spgline,1)

  'define RMSE4=sqrt(aave(v.6*v.6,lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE4'
  spaline = sublin(result,2)
  spa = subwrd(spaline,1)

  'define RMSE5=sqrt(aave((v.8-v.4)*(v.8-v.4),lon=0.1666667,lon=9.833333,lat=15.1667,lat=34.8333))'
  'd RMSE5'
  conline = sublin(result,2)
  con = subwrd(conline,1)

  write('vvel.txt', tmg' 'tma' 'spg' 'spa' 'con,append)
  t=t+1
endwhile
