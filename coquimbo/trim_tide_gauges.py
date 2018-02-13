from obspy import read
from glob import glob
from datetime import timedelta

path='/Users/dmelgar/Coquimbo2015/tsunami/sac/'
outpath='/Users/dmelgar/Coquimbo2015/tsunami/proc/chao/'
files=glob(path+'*.sac')
td=timedelta(seconds=7200)
for k in range(len(files)):
    sta=files[k].split('/')[-1].split('.')[0]
    t=read(files[k])
    #if sta=='buca':
    #    td=timedelta(seconds=4395)
    #if sta=='cald':
    #    td=timedelta(seconds=5951)
    #if sta=='cons':
    #    td=timedelta(seconds=6000)
    #if sta=='coqu':
    #    td=timedelta(seconds=3951)
    #if sta=='huas':
    #    td=timedelta(seconds=3801) #BAD
    #if sta=='juan':
    #    td=timedelta(seconds=4347)
    #if sta=='pich':
    #    td=timedelta(seconds=2370)
    #if sta=='quin':
    #    td=timedelta(seconds=2975)
    #if sta=='sanf':
    #    td=timedelta(seconds=5612)
    #if sta=='valp':
    #    td=timedelta(seconds=2967)
    if sta=='buca':
        td=timedelta(seconds=54*60)
    if sta=='cald':
        td=timedelta(seconds=57*60)
    if sta=='cons':
        td=timedelta(seconds=76*60)
    if sta=='coqu':
        td=timedelta(seconds=38*60)
    if sta=='huas':
        td=timedelta(seconds=42*60) #BAD
    if sta=='juan':
        td=timedelta(seconds=4347)
    if sta=='pich':
        td=timedelta(seconds=26*60)
    if sta=='quin':
        td=timedelta(seconds=100*60)
    if sta=='sanf':
        td=timedelta(seconds=90)
    if sta=='valp':
        td=timedelta(seconds=42*60)
    
    t[0].trim(endtime=t[0].stats.starttime+td)
    t.write(outpath+'c'+sta+'.tsun',format='SAC')
    