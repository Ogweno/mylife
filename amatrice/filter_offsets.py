from numpy import genfromtxt,argmin,savetxt

data=genfromtxt('/Users/dmelgar/code/GMT/Amatrice/offsets.txt')
synth=genfromtxt('/Users/dmelgar/code/GMT/Amatrice/gps_synth.txt')
out=synth.copy()
for k in range(len(synth)):
    dist=((synth[k,0]-data[:,0])**2+(synth[k,1]-data[:,1])**2)**0.5
    i=argmin(dist)
    out[k,:]=data[i,:]


savetxt('/Users/dmelgar/code/GMT/Amatrice/offsets_filt.txt',out,fmt='%.6f')
