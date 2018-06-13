* Dr. Stephen G. Penny, University of Maryland, 30-Sep-2011
* Visiting Scientist, NOAA/NCEP
*DSET ^ocean_temp_salt.res.nc
*DSET ^%y4%m2%d2.000000.ocean_temp_salt.res.nc
DSET ^%y4%m2%d2%h2.grd
OPTIONS template big_endian
*DTYPE netcdf
TITLE MOM4 MODEL restart data in letkf.grd format
*TITLE Ocean Temperature/Salinity restart
UNDEF 0.0
*UNDEF 0

XDEF 360 LEVELS -279.5, -278.5, -277.5, -276.5, -275.5, -274.5, -273.5, -272.5,
    -271.5, -270.5, -269.5, -268.5, -267.5, -266.5, -265.5, -264.5, -263.5,
    -262.5, -261.5, -260.5, -259.5, -258.5, -257.5, -256.5, -255.5, -254.5,
    -253.5, -252.5, -251.5, -250.5, -249.5, -248.5, -247.5, -246.5, -245.5,
    -244.5, -243.5, -242.5, -241.5, -240.5, -239.5, -238.5, -237.5, -236.5,
    -235.5, -234.5, -233.5, -232.5, -231.5, -230.5, -229.5, -228.5, -227.5,
    -226.5, -225.5, -224.5, -223.5, -222.5, -221.5, -220.5, -219.5, -218.5,
    -217.5, -216.5, -215.5, -214.5, -213.5, -212.5, -211.5, -210.5, -209.5,
    -208.5, -207.5, -206.5, -205.5, -204.5, -203.5, -202.5, -201.5, -200.5,
    -199.5, -198.5, -197.5, -196.5, -195.5, -194.5, -193.5, -192.5, -191.5,
    -190.5, -189.5, -188.5, -187.5, -186.5, -185.5, -184.5, -183.5, -182.5,
    -181.5, -180.5, -179.5, -178.5, -177.5, -176.5, -175.5, -174.5, -173.5,
    -172.5, -171.5, -170.5, -169.5, -168.5, -167.5, -166.5, -165.5, -164.5,
    -163.5, -162.5, -161.5, -160.5, -159.5, -158.5, -157.5, -156.5, -155.5,
    -154.5, -153.5, -152.5, -151.5, -150.5, -149.5, -148.5, -147.5, -146.5,
    -145.5, -144.5, -143.5, -142.5, -141.5, -140.5, -139.5, -138.5, -137.5,
    -136.5, -135.5, -134.5, -133.5, -132.5, -131.5, -130.5, -129.5, -128.5,
    -127.5, -126.5, -125.5, -124.5, -123.5, -122.5, -121.5, -120.5, -119.5,
    -118.5, -117.5, -116.5, -115.5, -114.5, -113.5, -112.5, -111.5, -110.5,
    -109.5, -108.5, -107.5, -106.5, -105.5, -104.5, -103.5, -102.5, -101.5,
    -100.5, -99.5, -98.5, -97.5, -96.5, -95.5, -94.5, -93.5, -92.5, -91.5,
    -90.5, -89.5, -88.5, -87.5, -86.5, -85.5, -84.5, -83.5, -82.5, -81.5,
    -80.5, -79.5, -78.5, -77.5, -76.5, -75.5, -74.5, -73.5, -72.5, -71.5,
    -70.5, -69.5, -68.5, -67.5, -66.5, -65.5, -64.5, -63.5, -62.5, -61.5,
    -60.5, -59.5, -58.5, -57.5, -56.5, -55.5, -54.5, -53.5, -52.5, -51.5,
    -50.5, -49.5, -48.5, -47.5, -46.5, -45.5, -44.5, -43.5, -42.5, -41.5,
    -40.5, -39.5, -38.5, -37.5, -36.5, -35.5, -34.5, -33.5, -32.5, -31.5,
    -30.5, -29.5, -28.5, -27.5, -26.5, -25.5, -24.5, -23.5, -22.5, -21.5,
    -20.5, -19.5, -18.5, -17.5, -16.5, -15.5, -14.5, -13.5, -12.5, -11.5,
    -10.5, -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5,
    1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5,
    14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5,
    26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5,
    38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5,
    50.5, 51.5, 52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5, 59.5, 60.5, 61.5,
    62.5, 63.5, 64.5, 65.5, 66.5, 67.5, 68.5, 69.5, 70.5, 71.5, 72.5, 73.5,
    74.5, 75.5, 76.5, 77.5, 78.5, 79.5
YDEF 200 LEVELS -81.5, -80.5, -79.5, -78.5, -77.5, -76.5, -75.5, -74.5, -73.5,
    -72.5, -71.5, -70.5, -69.5, -68.5, -67.5, -66.5, -65.5, -64.5, -63.5,
    -62.5, -61.5, -60.5, -59.5, -58.5, -57.5, -56.5, -55.5, -54.5, -53.5,
    -52.5, -51.5, -50.5, -49.5, -48.5, -47.5, -46.5, -45.5, -44.5, -43.5,
    -42.5, -41.5, -40.5, -39.5, -38.5, -37.5, -36.5, -35.5, -34.5, -33.5,
    -32.5, -31.5, -30.5, -29.5, -28.5014266967773, -27.5071048736572,
    -26.5197925567627, -25.5421199798584, -24.5765609741211,
    -23.6253776550293, -22.6905841827393, -21.7739162445068,
    -20.876802444458, -20.0003337860107, -19.1452465057373,
    -18.3119125366211, -17.5003337860107, -16.7101364135742,
    -15.9405832290649, -15.1905832290649, -14.4587106704712,
    -13.7432279586792, -13.0421209335327, -12.3531246185303,
    -11.673770904541, -11.001425743103, -10.3333330154419, -9.66666603088379,
    -9.00205135345459, -8.34354209899902, -7.69504070281982,
    -7.06020450592041, -6.44235372543335, -5.84438943862915,
    -5.26872444152832, -4.7172212600708, -4.19114875793457,
    -3.69114899635315, -3.21722149848938, -2.76872420310974,
    -2.34438943862915, -1.94235360622406, -1.56020474433899,
    -1.1950409412384, -0.843542039394379, -0.502051472663879,
    -0.166666254401207, 0.16666704416275, 0.502052307128906,
    0.843542814254761, 1.19504177570343, 1.56020557880402, 1.9423543214798,
    2.34439015388489, 2.76872515678406, 3.21722221374512, 3.69114971160889,
    4.19114971160889, 4.71722221374512, 5.26872491836548, 5.84438991546631,
    6.44235420227051, 7.06020545959473, 7.69504165649414, 8.34354305267334,
    9.00205230712891, 9.66666698455811, 10.3333339691162, 11.0014266967773,
    11.6737718582153, 12.3531255722046, 13.0421209335327, 13.7432289123535,
    14.4587116241455, 15.1905841827393, 15.9405841827393, 16.7101364135742,
    17.5003337860107, 18.3119125366211, 19.1452465057373, 20.0003337860107,
    20.8768043518066, 21.7739181518555, 22.6905841827393, 23.6253776550293,
    24.5765628814697, 25.542121887207, 26.5197925567627, 27.5071048736572,
    28.5014266967773, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5,
    38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5,
    50.5, 51.5, 52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5, 59.5, 60.5, 61.5,
    62.5, 63.5, 64.5, 65.5, 66.5, 67.5, 68.5, 69.5, 70.5, 71.5, 72.5, 73.5,
    74.5, 75.5, 76.5, 77.5, 78.5, 79.5, 80.5, 81.5, 82.5, 83.5, 84.5, 85.5,
    86.5, 87.5, 88.5, 89.5
ZDEF 50 LEVELS 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155,
    165, 175, 185, 195, 205, 215, 225, 236.122817993164, 250.599975585938,
    270.620819091797, 298.304931640625, 335.675628662109, 384.63427734375,
    446.936645507812, 524.170593261719, 617.736328125, 728.828491210938,
    858.421508789062, 1007.25708007812, 1175.83483886719, 1364.40625,
    1572.97131347656, 1801.27868652344, 2048.82861328125, 2314.87915039062,
    2598.45629882812, 2898.365234375, 3213.20581054688, 3541.38989257812,
    3881.162109375, 4230.62060546875, 4587.74267578125, 4950.40869140625,
    5316.4287109375

TDEF 1000 LINEAR 2jan1990 1dy

VARS 7
U 31 99 U-wind [m/s]
V 31 99 V-wind [m/s]
T 31 99 Temperature [C]
S 31 99 Specific Humidity [psu]
SSH 0 99 Surface Height [cm]
SST 0 99 Surface Temperature [C]
SSS 0 99 Surface Salinity [psu]
ENDVARS

