DSET ^%y4%m2%d2%h2/model_prep/%e/SFR_0%e_daily_dsw.nc
*DSET ^%y4%m2%d2%h2/model/%e/INPUT/ocean_temp_salt.res.nc
* Dr. Stephen G. Penny, University of Maryland, 21-May-2011
DTYPE netcdf
TITLE GEFS Surface Forcing
UNDEF 0.0
OPTIONS template

XDEF 576 LEVELS 0, 0.625, 1.25, 1.875, 2.5, 3.125, 3.75, 4.375, 5, 5.625, 6.25, 
    6.875, 7.5, 8.125, 8.75, 9.375, 10, 10.625, 11.25, 11.875, 12.5, 13.125, 
    13.75, 14.375, 15, 15.625, 16.25, 16.875, 17.5, 18.125, 18.75, 19.375, 
    20, 20.625, 21.25, 21.875, 22.5, 23.125, 23.75, 24.375, 25, 25.625, 
    26.25, 26.875, 27.5, 28.125, 28.75, 29.375, 30, 30.625, 31.25, 31.875, 
    32.5, 33.125, 33.75, 34.375, 35, 35.625, 36.25, 36.875, 37.5, 38.125, 
    38.75, 39.375, 40, 40.625, 41.25, 41.875, 42.5, 43.125, 43.75, 44.375, 
    45, 45.625, 46.25, 46.875, 47.5, 48.125, 48.75, 49.375, 50, 50.625, 
    51.25, 51.875, 52.5, 53.125, 53.75, 54.375, 55, 55.625, 56.25, 56.875, 
    57.5, 58.125, 58.75, 59.375, 60, 60.625, 61.25, 61.875, 62.5, 63.125, 
    63.75, 64.375, 65, 65.625, 66.25, 66.875, 67.5, 68.125, 68.75, 69.375, 
    70, 70.625, 71.25, 71.875, 72.5, 73.125, 73.75, 74.375, 75, 75.625, 
    76.25, 76.875, 77.5, 78.125, 78.75, 79.375, 80, 80.625, 81.25, 81.875, 
    82.5, 83.125, 83.75, 84.375, 85, 85.625, 86.25, 86.875, 87.5, 88.125, 
    88.75, 89.375, 90, 90.625, 91.25, 91.875, 92.5, 93.125, 93.75, 94.375, 
    95, 95.625, 96.25, 96.875, 97.5, 98.125, 98.75, 99.375, 100, 100.625, 
    101.25, 101.875, 102.5, 103.125, 103.75, 104.375, 105, 105.625, 106.25, 
    106.875, 107.5, 108.125, 108.75, 109.375, 110, 110.625, 111.25, 111.875, 
    112.5, 113.125, 113.75, 114.375, 115, 115.625, 116.25, 116.875, 117.5, 
    118.125, 118.75, 119.375, 120, 120.625, 121.25, 121.875, 122.5, 123.125, 
    123.75, 124.375, 125, 125.625, 126.25, 126.875, 127.5, 128.125, 128.75, 
    129.375, 130, 130.625, 131.25, 131.875, 132.5, 133.125, 133.75, 134.375, 
    135, 135.625, 136.25, 136.875, 137.5, 138.125, 138.75, 139.375, 140, 
    140.625, 141.25, 141.875, 142.5, 143.125, 143.75, 144.375, 145, 145.625, 
    146.25, 146.875, 147.5, 148.125, 148.75, 149.375, 150, 150.625, 151.25, 
    151.875, 152.5, 153.125, 153.75, 154.375, 155, 155.625, 156.25, 156.875, 
    157.5, 158.125, 158.75, 159.375, 160, 160.625, 161.25, 161.875, 162.5, 
    163.125, 163.75, 164.375, 165, 165.625, 166.25, 166.875, 167.5, 168.125, 
    168.75, 169.375, 170, 170.625, 171.25, 171.875, 172.5, 173.125, 173.75, 
    174.375, 175, 175.625, 176.25, 176.875, 177.5, 178.125, 178.75, 179.375, 
    180, 180.625, 181.25, 181.875, 182.5, 183.125, 183.75, 184.375, 185, 
    185.625, 186.25, 186.875, 187.5, 188.125, 188.75, 189.375, 190, 190.625, 
    191.25, 191.875, 192.5, 193.125, 193.75, 194.375, 195, 195.625, 196.25, 
    196.875, 197.5, 198.125, 198.75, 199.375, 200, 200.625, 201.25, 201.875, 
    202.5, 203.125, 203.75, 204.375, 205, 205.625, 206.25, 206.875, 207.5, 
    208.125, 208.75, 209.375, 210, 210.625, 211.25, 211.875, 212.5, 213.125, 
    213.75, 214.375, 215, 215.625, 216.25, 216.875, 217.5, 218.125, 218.75, 
    219.375, 220, 220.625, 221.25, 221.875, 222.5, 223.125, 223.75, 224.375, 
    225, 225.625, 226.25, 226.875, 227.5, 228.125, 228.75, 229.375, 230, 
    230.625, 231.25, 231.875, 232.5, 233.125, 233.75, 234.375, 235, 235.625, 
    236.25, 236.875, 237.5, 238.125, 238.75, 239.375, 240, 240.625, 241.25, 
    241.875, 242.5, 243.125, 243.75, 244.375, 245, 245.625, 246.25, 246.875, 
    247.5, 248.125, 248.75, 249.375, 250, 250.625, 251.25, 251.875, 252.5, 
    253.125, 253.75, 254.375, 255, 255.625, 256.25, 256.875, 257.5, 258.125, 
    258.75, 259.375, 260, 260.625, 261.25, 261.875, 262.5, 263.125, 263.75, 
    264.375, 265, 265.625, 266.25, 266.875, 267.5, 268.125, 268.75, 269.375, 
    270, 270.625, 271.25, 271.875, 272.5, 273.125, 273.75, 274.375, 275, 
    275.625, 276.25, 276.875, 277.5, 278.125, 278.75, 279.375, 280, 280.625, 
    281.25, 281.875, 282.5, 283.125, 283.75, 284.375, 285, 285.625, 286.25, 
    286.875, 287.5, 288.125, 288.75, 289.375, 290, 290.625, 291.25, 291.875, 
    292.5, 293.125, 293.75, 294.375, 295, 295.625, 296.25, 296.875, 297.5, 
    298.125, 298.75, 299.375, 300, 300.625, 301.25, 301.875, 302.5, 303.125, 
    303.75, 304.375, 305, 305.625, 306.25, 306.875, 307.5, 308.125, 308.75, 
    309.375, 310, 310.625, 311.25, 311.875, 312.5, 313.125, 313.75, 314.375, 
    315, 315.625, 316.25, 316.875, 317.5, 318.125, 318.75, 319.375, 320, 
    320.625, 321.25, 321.875, 322.5, 323.125, 323.75, 324.375, 325, 325.625, 
    326.25, 326.875, 327.5, 328.125, 328.75, 329.375, 330, 330.625, 331.25, 
    331.875, 332.5, 333.125, 333.75, 334.375, 335, 335.625, 336.25, 336.875, 
    337.5, 338.125, 338.75, 339.375, 340, 340.625, 341.25, 341.875, 342.5, 
    343.125, 343.75, 344.375, 345, 345.625, 346.25, 346.875, 347.5, 348.125, 
    348.75, 349.375, 350, 350.625, 351.25, 351.875, 352.5, 353.125, 353.75, 
    354.375, 355, 355.625, 356.25, 356.875, 357.5, 358.125, 358.75, 359.375

YDEF 288 LEVELS -89.522, -88.904, -88.281, -87.658, -87.035, -86.411, -85.787, 
    -85.164, -84.54, -83.916, -83.292, -82.668, -82.044, -81.421, -80.797, 
    -80.173, -79.549, -78.925, -78.301, -77.677, -77.053, -76.429, -75.806, 
    -75.182, -74.558, -73.934, -73.31, -72.686, -72.062, -71.438, -70.814, 
    -70.19, -69.566, -68.943, -68.319, -67.695, -67.071, -66.447, -65.823, 
    -65.199, -64.575, -63.951, -63.327, -62.703, -62.08, -61.456, -60.832, 
    -60.208, -59.584, -58.96, -58.336, -57.712, -57.088, -56.464, -55.84, 
    -55.217, -54.593, -53.969, -53.345, -52.721, -52.097, -51.473, -50.849, 
    -50.225, -49.601, -48.977, -48.353, -47.73, -47.106, -46.482, -45.858, 
    -45.234, -44.61, -43.986, -43.362, -42.738, -42.114, -41.49, -40.866, 
    -40.243, -39.619, -38.995, -38.371, -37.747, -37.123, -36.499, -35.875, 
    -35.251, -34.627, -34.003, -33.379, -32.756, -32.132, -31.508, -30.884, 
    -30.26, -29.636, -29.012, -28.388, -27.764, -27.14, -26.516, -25.893, 
    -25.269, -24.645, -24.021, -23.397, -22.773, -22.149, -21.525, -20.901, 
    -20.277, -19.653, -19.029, -18.406, -17.782, -17.158, -16.534, -15.91, 
    -15.286, -14.662, -14.038, -13.414, -12.79, -12.166, -11.542, -10.919, 
    -10.295, -9.671, -9.047, -8.423, -7.799, -7.175, -6.551, -5.927, -5.303, 
    -4.679, -4.055, -3.432, -2.808, -2.184, -1.56, -0.936, -0.312, 0.312, 
    0.936, 1.56, 2.184, 2.808, 3.432, 4.055, 4.679, 5.303, 5.927, 6.551, 
    7.175, 7.799, 8.423, 9.047, 9.671, 10.295, 10.919, 11.542, 12.166, 12.79, 
    13.414, 14.038, 14.662, 15.286, 15.91, 16.534, 17.158, 17.782, 18.406, 
    19.029, 19.653, 20.277, 20.901, 21.525, 22.149, 22.773, 23.397, 24.021, 
    24.645, 25.269, 25.893, 26.516, 27.14, 27.764, 28.388, 29.012, 29.636, 
    30.26, 30.884, 31.508, 32.132, 32.756, 33.379, 34.003, 34.627, 35.251, 
    35.875, 36.499, 37.123, 37.747, 38.371, 38.995, 39.619, 40.243, 40.866, 
    41.49, 42.114, 42.738, 43.362, 43.986, 44.61, 45.234, 45.858, 46.482, 
    47.106, 47.73, 48.353, 48.977, 49.601, 50.225, 50.849, 51.473, 52.097, 
    52.721, 53.345, 53.969, 54.593, 55.217, 55.84, 56.464, 57.088, 57.712, 
    58.336, 58.96, 59.584, 60.208, 60.832, 61.456, 62.08, 62.703, 63.327, 
    63.951, 64.575, 65.199, 65.823, 66.447, 67.071, 67.695, 68.319, 68.943, 
    69.566, 70.19, 70.814, 71.438, 72.062, 72.686, 73.31, 73.934, 74.558, 
    75.182, 75.806, 76.429, 77.053, 77.677, 78.301, 78.925, 79.549, 80.173, 
    80.797, 81.421, 82.044, 82.668, 83.292, 83.916, 84.54, 85.164, 85.787, 
    86.411, 87.035, 87.658, 88.281, 88.904, 89.522


ZDEF 1 LEVELS 1

TDEF 365 LINEAR 1jan2011 1dy

EDEF 20 names 01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20

VARS 1
dsw 0 t,y,x dsw [?]
ENDVARS

