fin=open('/Users/dmelgar/Cascadia/mesh/cascadia_2km_25deep.gmsh','r')

zs=-10

points=[]
while True:
    line=fin.readline()
    if 'Point' in line:
        z=int(float(line.split()[4].replace(',','')))
        print z
        if z==zs:
            points.append(int(line.split()[0].replace('Point(','').replace(')','')))
    if line=='':
        break
 
fin.close()       
print points
    