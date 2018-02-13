from numpy import savetxt,c_,r_
from obspy import read

sta='SYBC'
outfile='/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/data/waveforms/'+sta+'.txt'
e=read('/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/data/waveforms/'+sta+'.disp.e')
n=read('/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/data/waveforms/'+sta+'.disp.n')
u=read('/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/data/waveforms/'+sta+'.disp.u')
t1=n[0].stats.starttime

header='''Data provided by the Nepal Geodetic Array run by the California Institute of Technology
and the Departement of Mines and Geology (Nepal).

Raw data processed by the Scripps Permanent Array and Observation Center (SOPAC).

If you use this data please cite:
    
Galetzka et al.(2015), Slip pulse and resonance of Kathmandu basin during the 2015 
Mw 7.8 Gorkha earthquake, Nepal imaged with space geodesy, Science.

'''
header+='GPS station '+sta+'\nFirst sample is at '+str(t1)
header+='\n\nTime(s)  North(m)     East(m)      Up(m)'


out=c_[n[0].times(),n[0].data,e[0].data,u[0].data]
savetxt(outfile,out,fmt='%6.2f\t%10.4f\t%10.4f\t%10.4f',header=header)
