from obspy.taup import TauPyModel
from numpy import genfromtxt,arange,meshgrid,ones,zeros,save,array
from pyproj import Geod
from obspy.geodetics import kilometer2degrees
from scipy.linalg import norm
from obspy import read,UTCDateTime

model = TauPyModel(model="spica")
lonlat=genfromtxt('/Users/dmelgar/Slip_inv/puebla/data/station_info/sm.gflist',usecols=[1,2])
sta=genfromtxt('/Users/dmelgar/Slip_inv/puebla/data/station_info/sm.gflist',usecols=[0],dtype='S')

# coarse ettings
#time_epi=[UTCDateTime('2017-09-19T18:14:37'),UTCDateTime('2017-09-19T18:14:38'),UTCDateTime('2017-09-19T18:14:39'),UTCDateTime('2017-09-19T18:14:40'),UTCDateTime('2017-09-19T18:14:41'),UTCDateTime('2017-09-19T18:14:42')]
#xsource=arange(-99,-98.3,0.05)
#ysource=arange(18.2,18.6,0.05)
#zsource=arange(30,60,2)

# fine ettings
#time_epi=[UTCDateTime('2017-09-19T18:14:35.5'),UTCDateTime('2017-09-19T18:14:36'),UTCDateTime('2017-09-19T18:14:36.5'),UTCDateTime('2017-09-19T18:14:37'),UTCDateTime('2017-09-19T18:14:37.5'),UTCDateTime('2017-09-19T18:14:38')]
#xsource=arange(-98.7,-98.6,0.01)
#ysource=arange(18.15,18.35,0.01)
#zsource=arange(55,58,1)

#time_epi=[UTCDateTime('2017-09-19T18:14:32'),UTCDateTime('2017-09-19T18:14:33'),UTCDateTime('2017-09-19T18:14:34'),UTCDateTime('2017-09-19T18:14:35'),UTCDateTime('2017-09-19T18:14:36'),UTCDateTime('2017-09-19T18:14:37')]
#xsource=array([-98.65])
#ysource=array([18.22])
#zsource=arange(55,70,2)

time_epi=[UTCDateTime('2017-09-19T18:14:35.5'),UTCDateTime('2017-09-19T18:14:36'),UTCDateTime('2017-09-19T18:14:36.5'),UTCDateTime('2017-09-19T18:14:37'),UTCDateTime('2017-09-19T18:14:37.5'),UTCDateTime('2017-09-19T18:14:38')]
xsource=arange(-98.7,-98.6,0.01)
ysource=arange(18.15,18.35,0.01)
zsource=array([50])



[xs,ys,zs]=meshgrid(xsource,ysource,zsource)
xs=xs.ravel()
ys=ys.ravel()
zs=zs.ravel()


p=Geod(ellps='WGS84')

L2=zeros((len(zs),len(time_epi)))

#read all stations once
for ksta in range(len(sta)):
    if ksta==0:
        st=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/_withPicks/'+sta[ksta]+'.HLZ.sac')
    else:
        st+=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/_withPicks/'+sta[ksta]+'.HLZ.sac')


for kepi in range(len(time_epi)):
    for k in range(len(xs)):
        
        print 'epi %d, source %d/%d' % (kepi,k,len(zs))
        
        #What are the ptimes??
        tobs=zeros(len(sta))
        
        #get distances from all sites to current source point
        az,baz,d_in_m=p.inv(lonlat[:,0],lonlat[:,1],xs[k]*ones(len(lonlat)),ys[k]*ones(len(lonlat)))
        d_in_deg=kilometer2degrees(d_in_m/1000.)
        
        tsynth=zeros(len(lonlat))
        for ksta in range(len(lonlat)):
            tp=st[ksta].stats.starttime+st[0].stats['sac']['a']
            tobs[ksta]=tp-time_epi[kepi]
            arrivals = model.get_travel_times(source_depth_in_km=zs[k],distance_in_degree=d_in_deg[ksta],phase_list=["P", "p"])
            tsynth[ksta]=arrivals[0].time
        L2[k,kepi]=norm(tsynth-tobs)
        
save('/Users/dmelgar/Puebla2017/velocity_model/fine_depth_50km',L2)
