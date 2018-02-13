from numpy import genfromtxt,unique,fliplr

f='/Users/dmelgar/GlarmS/socal/stations/chanfile_loc20.dat'
outfile='/Users/dmelgar/GlarmS/socal/stations/socal_gps.txt'

sta=genfromtxt(f,usecols=1,dtype='S')
lonlat=genfromtxt(f,usecols=[4,5])
lonlat=fliplr(lonlat)
a,i=unique(sta,return_index=True)
sta=sta[i]
lonlat=lonlat[i,:]

fout=open(outfile,'w')
for k in range(len(sta)):
    line='%s\t%12.6f\t%12.6f\n'%(sta[k],lonlat[k,0],lonlat[k,1])
    fout.write(line)
fout.close()