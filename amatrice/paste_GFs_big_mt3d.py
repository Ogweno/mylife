from numpy import genfromtxt,r_,zeros,unique,where,diff,c_
from matplotlib import pyplot as plt

#variables=['Mxx','Mxy','Mxz','Myy','Myz','Mzz']
#variables=['Mxy','Mxz','Myy','Myz','Mzz']
variables=['Mxx']
limits=[10421,11661]

for kvar in range(5):
    variable=variables[kvar]
    print variable
    fout=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.los'
    
    f0=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished1.los')
    f1=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished2.los')
    f2=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished3.los')
    
    print '... finished reading data'
    
    f0sta=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished1.los',usecols=0,dtype='S')
    f1sta=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished2.los',usecols=0,dtype='S')
    f2sta=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished3.los',usecols=0,dtype='S')
    
    print '... finished reading station names'
    
    i1=where(f0[:,1]<limits[0])[0]
    i2=where((f1[:,1]>=limits[0]) & (f1[:,1]<limits[1]))[0]
    i3=where(f2[:,1]>=limits[1])[0]
    
    out=r_[f0[i1,1:],f1[i2,1:],f2[i3,1:]]
    outsta=r_[f0sta[i1],f1sta[i2],f2sta[i3]]
    
    
    
    outfile=open(fout, 'w')
    outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
    
    print '... writing output file'
    
    for k in range(len(out)):
        # T40349	11692	13.3149	42.8287	13.2621	43.1000	13.0000	3.9962e-05
        line='%s\t%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (outsta[k],out[k,0],out[k,1],out[k,2],out[k,3],out[k,4],out[k,5],out[k,6])
        outfile.write(line)
        
    outfile.close()