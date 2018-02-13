from numpy import genfromtxt,ones,tile,sqrt,sum,log10

mts=genfromtxt('/Users/dmelgar/code/GMT/SStsunami/cmt_catalogue.txt')
ids=genfromtxt('/Users/dmelgar/code/GMT/SStsunami/cmt_catalogue.txt',usecols=12,dtype='S')

exponent=10**mts[:,9]
exponent=tile(exponent,(6,1)).T
mts[:,3:9]=mts[:,3:9]*exponent
M0=(1./sqrt(2))*sqrt(sum(mts[:,3:9]**2,axis=1))
M0=M0/1e7
Mw=(2./3)*(log10(M0)-9.1)

f=open(u'/Users/dmelgar/SStsunami/CMT_summary.txt','w')
f.write('#       ID          lon          lat     depth(km)     Mw\n')
for k in range(len(mts)):
    line='%14s\t%10.4f\t%10.4f\t%8.2f\t%6.2f\n' % (ids[k],mts[k,0],mts[k,1],mts[k,2],Mw[k])
    f.write(line)
f.close()