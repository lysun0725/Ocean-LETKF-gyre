#!/bin/sh
# Special compile script just for UMD's DT2

f90file=read_drifters.f90
f90file1=read_drifters_1.f90
exefile=rd.x
exefile1=rd1.x

a=0
while [ $a -lt 10 ]
do
  file1="drifters_out.nc.0000 $a"
  file2=drifters_out.nc.00000$a
  if [ -f "$file1" ]; then
    mv "$file1" $file2
  fi
  a=`expr $a + 1`
done

ifort -I$NETCDF_FORTRAN_INCDIR $f90file -o $exefile -L$NETCDF_FORTRAN_LIBDIR -lnetcdff -Wl,-rpath,$NETCDF_FORTRAN_LIBDIR

txtfile=drifters_inp.txt
outfile=drifters_inpc.txt

if [ -f $txtfile ]; then
  rm -rf $txtfile
fi
./rd.x -np 80 -o $txtfile > rd.out # "20" depends on how many processors you have used.

ifort -I$NETCDF_FORTRAN_INCDIR $f90file1 -o $exefile1 -L$NETCDF_FORTRAN_LIBDIR -lnetcdff -Wl,-rpath,$NETCDF_FORTRAN_LIBDIR

if [ -f $outfile ]; then
  rm -rf $outfile
fi

./rd1.x -np 80 -o $outfile > rd1.out





 
