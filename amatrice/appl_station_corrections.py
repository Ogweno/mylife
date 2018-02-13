from obspy import read

stations=['LSS','RM33','SPD','TERO','ASP','PTI','MNF','FEMA','FOS','TRE','SPM','NRC','AMT','GSA','PZI1']
corrections=[0.90,0.25,0.57,2.96,2.52,0.37,-0.19,-0.51,0.49,1.92,2.25,-0.1,-0.47,2.04,0.47]
rootpath='/Users/dmelgar/Amatrice2016/strong_motion/sac/'
outpath='/Users/dmelgar/Amatrice2016/strong_motion/sac_tshifted/'

for k in range(len(stations)):
    n=read(rootpath+stations[k]+'.HNN.sac')
    e=read(rootpath+stations[k]+'.HNE.sac')
    z=read(rootpath+stations[k]+'.HNZ.sac')
    n[0].stats.starttime=n[0].stats.starttime+corrections[k]
    e[0].stats.starttime=e[0].stats.starttime+corrections[k]
    z[0].stats.starttime=z[0].stats.starttime+corrections[k]
    
    #n[0].stats['sac']['nzsec']=n[0].stats['sac']['nzsec']+corrections[k]
    #e[0].stats['sac']['nzsec']=e[0].stats['sac']['nzsec']+corrections[k]
    #z[0].stats['sac']['nzsec']=z[0].stats['sac']['nzsec']+corrections[k]
    
    n.write(outpath+stations[k]+'.HNN.sac',format='SAC')
    e.write(outpath+stations[k]+'.HNE.sac',format='SAC')
    z.write(outpath+stations[k]+'.HNZ.sac',format='SAC')
