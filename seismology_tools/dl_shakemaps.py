#!/usr/bin/env python
import argparse
import sys
import httplib
from bs4 import BeautifulSoup
from urllib2 import urlopen,HTTPError

SHAKEMAP_URL="http://earthquake.usgs.gov/earthquakes/shakemap/list.php"
REGULAR="x=1"
SCENARIO="s=1"
YEAR="y=%d"
NETWORK="n=%s"
NETWORKS=["sc","nc","pn","nn","ut","ak","hv","nm","global"]
NETWORK_NAMES=["S. California","N. California","Pacific NW","Nevade","Utah","Alaska","Hawaii","New Madrid","Global"]

XML_URL='http://earthquake.usgs.gov/earthquakes/shakemap/%(network)s/shake/%(eventid)s/download/grid.xml'
TEXT_URL='http://earthquake.usgs.gov/earthquakes/shakemap/%(network)s/shake/%(eventid)s/download/grid.xyz.zip'

def download(url,filename, chunk_size=8192, report_hook=None):
    print("%s -> %s"%(url,filename))
    try:
        response = urlopen(url);
        total_size=None
        if 'content-length' in response.info():
            total_size = response.info().getheader('Content-Length')
            total_size = int(total_size.strip())
        bytes_so_far = 0

        with open(filename,'w') as f:
            while True:
                try:
                    chunk = response.read(chunk_size)
                except httplib.IncompleteRead, e:
                    return e.partial

                f.write(chunk)
                bytes_so_far += len(chunk)

                if not chunk: 
                    sys.stdout.write('\n')
                    break

                if total_size:
                    percent = round(100*float(bytes_so_far)/total_size,2)
                    sys.stdout.write("Downloading %d of %d bytes (%0.2f%%)\r" % (bytes_so_far, total_size, percent))
                else:
                    sys.stdout.write("Downloading %d bytes (size unknown)\r" % (bytes_so_far))
    except HTTPError as err:
        print('%d : %s'%(err.code,err.reason))
        
def extend(l,m):
    l.extend(m)
    return l

def getEvents(args):
    url="http://earthquake.usgs.gov/earthquakes/shakemap/list.php?"

    url_arg=[]
    events=[]
    if args.network: url_arg.append(NETWORK%args.network)
    if args.year: url_arg.append(YEAR%args.year)
    if args.real:
        events.extend(downloadEvents(url+'&'.join(extend([REGULAR],url_arg))))
    if args.scenario:
        events.extend(downloadEvents(url+'&'.join(extend([SCENARIO],url_arg))))
    return events

def downloadEvents(url):
    events=[]
    print('Download %s'%url)
    soup=BeautifulSoup(urlopen(url))
    soup.find("div",{"id":"main"})
    main_div=soup.find("div",{"id":"main"})
    tbody=main_div.find('tbody')
    for a in tbody.findAll('a'):
        href=a.attrs['href']
        parts=href.split('/')
        eventid=parts[-2]
        network=parts[-4]
        print(href)
        events.append({'eventid':eventid,'network':network})
    return events

if __name__=="__main__":
    parser = argparse.ArgumentParser('Batch Downloads USGS shakemaps data')
    parser.add_argument('-v',dest='verbose',action='store_true',help='Verbose output')
    parser.add_argument('-n',dest='network',nargs="?",choices=NETWORKS)
    parser.add_argument('-s',dest='scenario',action='store_true',help='Download scenario events')
    parser.add_argument('-r',dest='real',action='store_true',help='Download real events')
    parser.add_argument('-x',dest='xml',action='store_true',default=True,help='Download XML formatted data')
    parser.add_argument('-t',dest='text',action='store_true',default=True,help='Download text formatted data')
    parser.add_argument('-y',dest='year',nargs="?",help='Download events from a cetain year')

    args=parser.parse_args()

    events=getEvents(args)

    for event in events:
        if args.xml:
            download(XML_URL%event,event['eventid']+'.xml')
        if args.text:
            download(TEXT_URL%event,event['eventid']+'.xyz.zip')
