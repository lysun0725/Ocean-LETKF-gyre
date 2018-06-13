DSET ocean_velocity.res.nc
*DSET ^%y4%m2%d2.000000.ocean_velocity.res.nc
*DTYPE netcdf
TITLE Ocean U/V Velocity restart
OPTIONS template
*UNDEF 99999.0

XDEF 30 LEVELS 0.1666667, 0.5, 0.8333333, 1.166667, 1.5, 1.833333, 2.166667, 
    2.5, 2.833333, 3.166667, 3.5, 3.833333, 4.166667, 4.5, 4.833333, 
    5.166667, 5.5, 5.833333, 6.166667, 6.5, 6.833333, 7.166667, 7.5, 
    7.833333, 8.166667, 8.5, 8.833333, 9.166667, 9.5, 9.833333 
YDEF 60 LEVELS 15.16667, 15.5, 15.83333, 16.16667, 16.5, 16.83333, 17.16667, 
    17.5, 17.83333, 18.16667, 18.5, 18.83333, 19.16667, 19.5, 19.83333, 
    20.16667, 20.5, 20.83333, 21.16667, 21.5, 21.83333, 22.16667, 22.5, 
    22.83333, 23.16667, 23.5, 23.83333, 24.16667, 24.5, 24.83333, 25.16667, 
    25.5, 25.83333, 26.16667, 26.5, 26.83333, 27.16667, 27.5, 27.83333, 
    28.16667, 28.5, 28.83333, 29.16667, 29.5, 29.83333, 30.16667, 30.5, 
    30.83333, 31.16667, 31.5, 31.83333, 32.16667, 32.5, 32.83333, 33.16667, 
    33.5, 33.83333, 34.16667, 34.5, 34.83333
ZDEF 50 LEVELS 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 
    165, 175, 185, 195, 205, 215, 225, 235.5619, 250.0461, 269.52, 297.2319, 
    334.083, 383.1031, 444.9251, 522.265, 615.3997, 726.6511, 855.8699, 
    1004.924, 1173.189, 1362.042, 1570.357, 1799.008, 2046.37, 2312.824, 
    2596.269, 2896.635, 3211.393, 3540.079, 3879.808, 4229.801, 4586.908, 
    4950.129, 5316.14

TDEF Time 1000 linear 0z2jan1980 1dy

@ xaxis_1 String units degree_east 
@ xaxis_1 String long_name Nominal Longitude of T-cell center
@ yaxis_1 String units degree_north
@ yaxis_1 String long_name Nominal Latitude of T-cell center
@ zaxis_1 String units meters
@ zaxis_1 String long_name tz
@ Time String units days since 0001-01-01 00:00:00
@ Time String long_name Time

@ u String long_name initial u
@ u String units m/s
@ v String long_name initial v
@ v String units m/s
