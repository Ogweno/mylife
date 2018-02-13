#strike line coordinates
s1=[360-73.407557,-37.045842,360-72.136641,-34.167922]
s2=[360-73.96121,-36.903269,360-72.737874,-34.065193]
#dip lines
d1=[360-73.272363,-36.905475,360-74.415519,-36.579790]
d2=[360-72.898813,-36.107714,360-74.06815,-35.750574]
d3=[360-72.512036,-35.295536,360-73.655134,-34.970736]
d4=[360-72.146905,-34.656404,360-73.301863,-34.281861]

#plot([s1[0],s1[2]],[s1[1],s1[3]])
#plot([s2[0],s2[2]],[s2[1],s2[3]])
#plot([d1[0],d1[2]],[d1[1],d1[3]])
#plot([d2[0],d2[2]],[d2[1],d2[3]])
#plot([d3[0],d3[2]],[d3[1],d3[3]])
#plot([d4[0],d4[2]],[d4[1],d4[3]])

def make_line(p,npoints):
    from numpy import linspace
    m=(p[3]-p[1])/(p[2]-p[0])
    b=p[1]-m*p[0]
    x=linspace(p[2],p[0],npoints)
    y=m*x+b
    return x,y
    
xs1,ys1=make_line(s1,10)
xs2,ys2=make_line(s2,10)
xd1,yd1=make_line(d1,10)
xd2,yd2=make_line(d2,10)
xd3,yd3=make_line(d3,10)
xd4,yd4=make_line(d4,10)

def print_gauge(x,y,indices):
    from string import rjust
    for k in range(len(x)):
        print u'rundata.gaugedata.gauges.append([%d, %.6f,%.6f, 0., 1.e10])' %(indices[k],x[k],y[k])

indices=range(10,20)
print_gauge(xs1,ys1,indices)
indices=range(20,30)
print_gauge(xs2,ys2,indices)
indices=range(30,40)
print_gauge(xd1,yd1,indices)
indices=range(40,50)
print_gauge(xd2,yd2,indices)
indices=range(50,60)
print_gauge(xd3,yd3,indices)
indices=range(70,80)
print_gauge(xd4,yd4,indices)