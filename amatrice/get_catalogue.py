from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.event import read_events,Catalog

client = Client('INGV')
lat_center=42.81
lon_center=13.17
radius=0.6 #In degrees
starttime = UTCDateTime("2016-08-24")
endtime = UTCDateTime("2016-11-02")

get_catalogue=False
get_stations=True

#This gets the whole catalog
if get_catalogue==True:
    cat = client.get_events(starttime=starttime, endtime=endtime,includearrivals=False,
        latitude=lat_center,longitude=lon_center,maxradius=radius,minmagnitude=2)
    cat.write(u'/Users/dmelgar/Amatrice2016/afters/INGV_nopicks.xml',format='QUAKEML')

else: #read it
    cat=read_events(u'/Users/dmelgar/Amatrice2016/afters/INGV_nopicks.xml')
    
#loop over all events and build new catalogue, this time with picks  
Nevents=len(cat)
cat2=Catalog()
for kevent in range(Nevents):
    print 'Working on event %d of %d' % (kevent,Nevents)
    ev=cat[kevent]
    evID=int(ev.resource_id.id.split('=')[-1])
    try:
        cat_temp=client.get_events(eventid=evID,includearrivals=True)
        ev=cat_temp[0]
        cat2+=ev
    except:
        print '... Could not get event'
    
    
cat2.write(u'/Users/dmelgar/Amatrice2016/afters/INGV_withpicks.xml',format='QUAKEML')

if get_stations==True:
    inventory = client.get_stations(network="IV",starttime=starttime,endtime=endtime)
    inventory += client.get_stations(network="BA",starttime=starttime,endtime=endtime)
    inventory += client.get_stations(network="MN",starttime=starttime,endtime=endtime)
    inventory += client.get_stations(network="TV",starttime=starttime,endtime=endtime)
    inventory += client.get_stations(network="XO",starttime=starttime,endtime=endtime)
