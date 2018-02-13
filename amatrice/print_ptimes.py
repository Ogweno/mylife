from obspy.core import UTCDateTime
from obspy import read
from numpy import genfromtxt,where,arange,ones,zeros
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees

time_epi=UTCDateTime('2016-08-24T01:36:32')

stations=['LSS','RM33','SPD','TERO','ASP','PTI','MNF','FEMA','FOS','TRE','SPM','NRC','AMT']
Nsta=len(stations)
relative_ptimes=[5.831,65.7839,5.4050,67.8967,5.3254,6.5378,5.7208,63.1667,5.8395,5.9275,5.8311,5.825,5.66964]
rootpath='/Users/dmelgar/Amatrice2016/strong_motion/sac/'


for ksta in range(len(stations)):
    print stations[ksta]
    st=read(rootpath+stations[ksta]+'.HNZ.sac')
    print st[0].stats.starttime+relative_ptimes[ksta]