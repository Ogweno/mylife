'''
Read netcdf (GMT .grd) file and convert to stupid ESRI integer ASCII format
The order of the headers is native, you still need to flip them if the file
is going to be used by GeoClaw.

If you don't have the netCDF module jsut do $pip install netCDF4 at your terminal

'''

import sys
from netCDF4 import Dataset
from numpy import where,nan,flipud,isnan

# If you want input from the shell
#infile=sys.argv[1]
#outfile=sys.argv[2]

#Otherwise define your files here
infile='/Users/dmelgar/Lefkada2015/scripts/slip.grd'
outfile='/Users/dmelgar/Lefkada2015/scripts/slip_esri.txt'

print 'Converting from grd to ESRI ASCII...'
print '... input file is '+infile
print '... output file is '+outfile

grd = Dataset(infile, 'r', format='NETCDF4')

#Sometimes x is x, sometimes x is lon, deal with both
try:
    x=grd.variables['x'][:]
    y=grd.variables['y'][:]
    z=grd.variables['z'][:]
except:
    x=grd.variables['lon'][:]
    y=grd.variables['lat'][:]
    z=grd.variables['z'][:]

z=flipud(z) #Because GMT stores things from lower left corner

#get header crap
f=open(outfile,'w')
f.write('NCOLS '+str(len(x))+'\n')
f.write('NROWS '+str(len(y))+'\n')
f.write('XLLCENTER %.12f\n' %(x.min()))
f.write('YLLCENTER %.12f\n' %(y.min()))
f.write('CELLSIZE %.12f\n' %(x[1]-x[0]))
f.write('NO_DATA_VALUE -9999\n')

#Traverse grid
for k in range(len(y)):
    line=z[k,:].data
    i=where(isnan(line)==1)[0] #inf Nans
    line[i]=-9999 #Write nans as -9999
    a=line[0]
    lineout=str(float(a))
    for kw in range(1,len(line)):
        lineout=lineout+'\t'+str(float(line[kw]))
    f.write(lineout+'\n')
f.close()