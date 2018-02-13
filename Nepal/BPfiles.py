from numpy import genfromtxt,savetxt,where,c_,zeros,ones
from string import rjust

bpfile='/Users/dmelgar/code/GMT/Nepal/HFdotsAU'
outpath='/Users/dmelgar/code/GMT/Nepal/BPfiles/'
bp=genfromtxt(bpfile,usecols=[0,1,2,3])
#Write zero file
tempout=zeros((1,3))
tempout[0,0]=0
tempout[0,1]=90
tempout[0,2]=0
savetxt(outpath+'BP0000.txt',tempout,fmt='%.6f\t%.6f\t%.6f')
#Do all others
for k in range(1,145):
    t=0.5*k
    i=where((bp[:,0]-10<=t) & (bp[:,0]>0))[0]
    if len(i)<1:
        out=tempout.copy()
    else:
        out=c_[bp[i,2],bp[i,1],bp[i,0]-10]
    outfile=outpath+'BP'+rjust(str(k),4,'0')+'.txt'
    savetxt(outfile,out,fmt='%.6f\t%.6f\t%.6f')
    