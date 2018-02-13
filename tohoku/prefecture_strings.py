from glob import glob
prefs=glob(u'/Users/dmelgar/KMLs/japan/prefectures/JPN*')

#anto.lon,anto.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_antofagasta.txt')
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    print  p+'.lon,'+p+'.lat=read_region('+"'"+prefs[k]+"')"
    
#anto=region()
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    print  p+'=region()'
    
#regions=[anto,arec,atac,icax,moqu,tacn,tara]
#out='['
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    out=out+p+','
#print out

#if argmin(d)==0:
#    anto.eta.append(eta_out[k])
#f=open(u'/Users/dmelgar/KMLs/japan/prefectures/temp_string.txt','w')
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    f.write('if argmin(d)=='+str(k)+':\n')
#    f.write('\t'+p+'.eta.append(eta_out[k])\n')
#f.close()

#anto.color=warn_color(anto)
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    print  p+'.color=warn_color('+p+')'
    
#f.write('Antofagasta='+anto.color+'\n')
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    plong=prefs[k].split('_')[1].split('.')[0]
#    print  'f.write('+"'"+plong+'='+"'"+'+'+p+'.color+'+')'




# FOR GMT

#ctarapaca=`grep "Tarapaca" $warnfile | awk -F "=" '{print $2}'`
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    plong=prefs[k].split('_')[1].split('.')[0]
#    print  'c'+plong.lower()+'=`grep '+'"'+plong+'"'+' $warnfile | awk -F '+'"'+'='+'"'+' '+"'"+'{print $2}'+"'`"
    
#pol_tarapaca=/Users/dmelgar/KMLs/chile/regions/CHL_tarapaca.txt
#for k in range(len(prefs)):
#    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
#    plong=prefs[k].split('_')[1].split('.')[0]
#    print  'pol_'+plong.lower()+'=/Users/dmelgar/KMLs/japan/prefectures/JPN_'+plong+'.txt'
    
#gmt psxy $pol_tarapaca -W0.4p,black -G$ctarapaca -t$transp -R -J -O -K >> $nNAME.eps
for k in range(len(prefs)):
    p=prefs[k].split('_')[1].split('.')[0].lower()[0:4]
    plong=prefs[k].split('_')[1].split('.')[0]
    print  'psxy $pol'+plong.lower()+' -W0.4p,black -G$c'+plong.lower()+' -t$transp -R -J -O -K >> $nNAME.eps'