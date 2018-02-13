from numpy import genfromtxt,arange,r_,sort,array,savetxt
from matplotlib import pyplot as plt

f=genfromtxt('/Users/dmelgar/Slip_inv/Lefkada65/output/inverse_models/models/gps_sm_insar_5win_v2.6_refine.0003.inv')

i1=arange(0,40,3)
i2=arange(1,40,3)
i3=arange(2,40,3)

irow=array([ 0,  1,  2,  6,  7,  8, 12, 13, 14, 18, 19, 
                20, 24, 25, 26, 30, 31, 32, 36, 37, 38])
                
icol1=r_[irow,irow+40,irow+80]
icol3=r_[irow+240,irow+280,irow+320]

irow=array([ 0,  1,  2,  6,  7,  8, 12, 13, 14, 18, 19, 
                20, 24, 25, 26, 30, 31, 32,36])
icol2=r_[irow+123,irow+163,irow+203]
icol4=r_[irow+363,irow+403,irow+443]

i=r_[icol1,icol2,icol3,icol4]

f[:,9]=0
f[:,8]=0
#f[i,8]=-0.2
#f[i+1*480,8]=-0.2
#f[i+2*480,8]=-0.2
#f[i+3*480,8]=-0.2
#f[i+4*480,8]=-0.2

#one window
f[i,8]=-1.0
f=f[0:480,:]

fmt='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
savetxt(u'/Users/dmelgar/Slip_inv/Lefkada_fwd/forward_models/checkerboard_1win.rupt',f,fmt)

iplot=arange(0,400)
plt.scatter(f[iplot,1],f[iplot,2],c=f[iplot,8],s=50,lw=0.5,cmap=plt.cm.magma_r)
plt.colorbar()
plt.show()

