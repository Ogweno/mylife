from numpy import genfromtxt,savetxt,log10,c_,zeros
from obspy.core import UTCDateTime

#YEAR	MON	DAY	HOUR	MIN	SEC	LAT	LON	DEP	MW	Stk1	Dip1	Rake1	Stk2	Dip2	Rake2	Ptr	Ppl	Ttr	Tpl	Ntr	Npl	Regime	Fit
#
#	 (c) Focal mechanisms in Harvard CMT convention
#	     X, Y, depth, strike1, dip1, rake1, strike2, dip2, rake2, moment, newX, newY, event_title
#	     with moment in 2 columns : mantissa and exponent corresponding to seismic moment in dynes-cm

c=genfromtxt('/Users/dmelgar/Amatrice2016/afters/sluitcatalog.csv')
out='/Users/dmelgar/Amatrice2016/afters/sluitcatalog.psmeca'
#c=genfromtxt('/Users/dmelgar/Amatrice2016/afters/slu_2016_large.csv')
#out='/Users/dmelgar/Amatrice2016/afters/slu_2016_large.psmeca'

def frexp_10(decimal):
   logdecimal = log10(decimal)
   return 10 ** (logdecimal - int(logdecimal)), int(logdecimal)

moment=10**(1.5*c[:,9]+9.1)
exponent=zeros(len(c))
mantissa=zeros(len(c))
for k in range(len(c)):
    mantissa[k],exponent[k]=frexp_10(moment[k])


#Days  since origin
days=zeros(len(c))
t0=UTCDateTime('%s-%s-%s' % (int(c[0,0]),int(c[0,1]),int(c[0,2])))
for k in range(len(c)-1):
    t1=UTCDateTime('%s-%s-%s' % (int(c[k+1,0]),int(c[k+1,1]),int(c[k+1,2])))
    days[k+1]=(t1-t0)/86400
    

c_out=c_[c[:,7],c[:,6],days,c[:,10:13],c[:,13:16],mantissa,exponent]
savetxt(out,c_out,fmt='%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%d')