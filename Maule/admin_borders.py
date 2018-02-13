#Count number of regions
#from numpy import unique
#
#f=open('/Users/dmelgar/KMLs/chile/chile/CHL_adm1.txt','r')
#lineout=[]
#for line in f:
#    if line[0]=='>':
#        lineout.append(line)
#f.close()
#print unique(lineout)
        
##Separate into individual files
#['> -L"Antofagasta"\n' '> -L"Araucania"\n' '> -L"Atacama"\n'
# '> -L"Aysen"\n' '> -L"Bio-Bio"\n' '> -L"Coquimbo"\n' '> -L"Los Lagos"\n'
# '> -L"Magallanes"\n' '> -L"Maule"\n' '> -L"O\'Higgins"\n'
# '> -L"Santiago"\n' '> -L"Tarapaca"\n' '> -L"Valparaiso"\n']

test='Valpa'
#g=open('/Users/dmelgar/KMLs/chile/CHL_antofagasta.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_aracuania.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_atacama.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_aysen.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_biobio.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_coquimbo.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_loslagos.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_magallanes.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_maule.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_ohigins.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_santiago.txt','w')
#g=open('/Users/dmelgar/KMLs/chile/CHL_tarapaca.txt','w')
g=open('/Users/dmelgar/KMLs/chile/CHL_valparaiso.txt','w')

f=open('/Users/dmelgar/KMLs/chile/CHL_adm1.txt','r')
main_file=True
write_line=False 
k=0
while main_file:
    line=f.readline()
    if '>' in line: #Make a decision
        if test in line:
            write_line=True
        else:
            write_line=False
    if write_line:
        g.write(line)
    if line=='': #End of file
        main_file=False
    k+=1
f.close()
g.close()