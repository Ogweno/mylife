from mudpy import analysis
from glob import glob

folders=glob('/Users/dmelgar/FakeQuakes/Cascadia_final1/output/dreger_4tau/cascadia*')
vel_model='/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia'
gf_list=u'/Users/dmelgar/FakeQuakes/Cascadia_final1/data/station_info/cascadia_small.gflist'

for k in range(163,len(folders)):
    event=folders[k].split('/')[-1]
    print event
    event_log=folders[k]+'/_'+event+'.log'
    out_file='/Users/dmelgar/FakeQuakes/Cascadia_final1/output/dreger_4tau/_picks/'+event+'.picks'
    analysis.dump_picks(event_log,vel_model,gf_list,out_file)