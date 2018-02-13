from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.event import read_events,Catalog
from numpy import arange
from string import rjust


#Read the catalogue
print 'Loading humongous catalog...'
cat=read_events(u'/Users/dmelgar/Amatrice2016/afters/INGV_withpicks.xml')

out_picks='/Users/dmelgar/Amatrice2016/afters/catalogs/with_picks/'
out_nopicks='/Users/dmelgar/Amatrice2016/afters/catalogs/no_picks/'

#Break the catalogue up into many more because it takes too long to read.
Nevents=len(cat)
Nbreak=40
Nperbreak=Nevents/Nbreak

for k in range(Nbreak):
    
    event_range=arange(k*Nperbreak,(k+1)*Nperbreak)
    
    print event_range
    
    if event_range[-1]>Nevents-1:
        cat_picks=cat[k*Nperbreak:Nevents]
        cat_picks.write(out_picks+'INGV_withpicks.'+rjust(str(k),4,'0')+'.xml',format='QUAKEML')
        
        cat_nopicks=cat[k*Nperbreak:Nevents]
        for kevent in range(len(cat_nopicks)):
            cat_nopicks[kevent].picks=[]
        cat_nopicks.write(out_nopicks+'INGV_nopicks.'+rjust(str(k),4,'0')+'.xml',format='QUAKEML')
    else:
        cat_picks=cat[k*Nperbreak:(k+1)*Nperbreak]
        cat_picks.write(out_picks+'INGV_withpicks.'+rjust(str(k),4,'0')+'.xml',format='QUAKEML')
        
        cat_nopicks=cat[k*Nperbreak:(k+1)*Nperbreak]
        for kevent in range(len(cat_nopicks)):
            cat_nopicks[kevent].picks=[]
        cat_nopicks.write(out_nopicks+'INGV_nopicks.'+rjust(str(k),4,'0')+'.xml',format='QUAKEML')
