from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.event import read_events,Catalog
from glob import glob

#Path to write waveforms to:
out='/Users/dmelgar/Amatrice2016/afters/SEED/'

#Path to all catalogs
cat_path=u'/Users/dmelgar/Amatrice2016/afters/catalogs/with_picks/'

pwindow=5. #in s.
swindow=10. #in s.

#Init the client
client = Client('INGV')

#How many catalogs are there?
catalogs=glob(cat_path+'*.xml')

for kcatalog in range(22,len(catalogs)):
    
    print 'workin on catalog '+catalogs[kcatalog]
    cat=read_events(catalogs[kcatalog])
 
    #Loop through events in this catalog
    for kevent in range(len(cat)):
        
        if kevent % 10 == 0:
            print '... processing event %d of %d' % (kevent,len(cat))
        
        #Current event
        ev=cat[kevent]
        
        #loop through picks
        for kpick in range(len(ev.picks)):
        
            if kpick % 40 == 0:
                print '... .. pick %d of %d' % (kpick,len(ev.picks))
        
            #Current pick
            pick=ev.picks[kpick]
            
            #Extract time and station info
            t=pick.time
            net=pick.waveform_id.network_code
            sta=pick.waveform_id.station_code
            chan=pick.waveform_id.channel_code
            loc=pick.waveform_id.location_code
            phase=pick.phase_hint
            
            try:
                #Go fetch that data
                if phase=='P':
                    st = client.get_waveforms(net,sta,loc,chan,t-pwindow/2,t+pwindow/2)
                    st[0].trim(starttime=t-pwindow/2,endtime=t+pwindow/2)
                elif phase=='S':
                    st = client.get_waveforms(net,sta,loc,chan,t-swindow/2,t+swindow/2)
                    st[0].trim(starttime=t-swindow/2,endtime=t+swindow/2)
                
                timestring=t.strftime('%Y%m%d%H%M%S')
                filename=out+net+'.'+sta+'.'+chan+'.'+timestring+'.mseed'
                st.write(filename,format='MSEED')
            except:
                print '... ... ... woah, no data!'
