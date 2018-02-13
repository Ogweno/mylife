from numpy import genfromtxt

g=genfromtxt('/Users/dmelgar/Puebla2017/GPS/static/unr.txt')
sta=genfromtxt('/Users/dmelgar/Puebla2017/GPS/static/unr.txt',usecols=0,dtype='S')


for k in range(len(sta)):
    f=open("/Users/dmelgar/Slip_inv/puebla/data/statics/"+sta[k]+'.neu','w')
    f.write('%.6f\n%.6f\n%.6f' % (g[k,4]/1000,g[k,3]/1000,g[k,5]/1000))
    f.close()