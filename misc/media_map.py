'''
D. Melgar 6/2016

Create simple map for media inquiries
'''

def run_map(out_file,epicenter,path_to_aftershocks,path_to_shakemap,path_to_faults,zoom_start_level=10,
            plot_afters=True,plot_shakemap=True,plot_faults=True,plot_mmi=True):

    '''
    Use folium to make the map
    
    IN:
        out_file: String, path and name of the html file
        epicenter: List or tuple, lat,lon of epicenter
        path_to_aftershocks: String, path and filename of file with aftershocks
        path_to_shakemap: String, path and file name of xyz shakemap
        path_to_faults: String, path and file name of fault trace coordinates file
        zoom_start_level: Int (optional), zoom level when map is first opened
        
    RETURNS:
        Nothing
    '''
    

    import folium
    from numpy import genfromtxt,fliplr,zeros,unique
    from matplotlib import cm
    from matplotlib.pyplot import tricontour
    
    
    #Map attributes
    tileset=r'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png'
    attribution='Poop'
    
    #Aftershocks
    afters=genfromtxt(path_to_aftershocks)
    afters=fliplr(afters)
    
    #Shakemap
    shakemap=genfromtxt(path_to_shakemap,usecols=[1,0,4])
    lat1=shakemap[:,0].min()
    lat2=shakemap[:,0].max()
    lon1=shakemap[:,1].min()
    lon2=shakemap[:,1].max()
    shake_bounds=[lat1,lat2,lon1,lon2]
    
    
    #Process shakemap
    k=0
    MMI=zeros((len(unique(shakemap[:,0])),len(unique(shakemap[:,1]))))
    for klat in range(len(unique(shakemap[:,0]))):
        for klon in range(len(unique(shakemap[:,1]))):
                MMI[klat,klon]=shakemap[k,2]
                k+=1
    MMIrgb= cm.jet(MMI/MMI.max())
    
    
    #Make map
    map_bsl=folium.Map(location=epicenter,tiles=tileset,attr=attribution,zoom_start=zoom_start_level,control_scale=True)
    
    #Epicenter
    folium.Marker(location=epicenter,icon=folium.Icon(color='red',icon='glyphicon glyphicon-asterisk'),popup='August 28, Mw 6.1').add_to(map_bsl)
    
    
    #Shakemap
    if plot_shakemap:
        map_bsl.image_overlay(MMIrgb,min_lat=shake_bounds[0], max_lat=shake_bounds[1], min_lon=shake_bounds[2], max_lon=shake_bounds[3],opacity=0.5,mercator_project=True)
        #map_bsl.add_children(plugins.ImageOverlay(MMIrgb,bounds=shake_bounds,opacity=0.5,mercator_project=True))
    
    #Add ucerf faults
    if plot_faults:
        ucerf=open(path_to_faults,'r')
        ucerf.readline()
        lonlat=[]
        all_lines=[]
        while True:
            line=ucerf.readline()
            if '>' in line:
                all_lines.append(lonlat)
                lonlat=[]
            elif line=='':
                break
            else:
                xy=(float(line.split()[1]),float(line.split()[0]))
                lonlat.append(xy)
        ucerf_lines = folium.MultiPolyLine(locations=all_lines, color='#800000', weight=2, opacity=0.5)
        map_bsl.add_children(ucerf_lines)    

    #Now shakemap MMI4 and MMI7 contours
    if plot_mmi:
        mmi_contours=tricontour(shakemap[:,0],shakemap[:,1],shakemap[:,2],levels=[6,7,8])
        mmi_lines=[mmi_contours.collections[0].get_paths()[0].vertices.tolist()]
        mmi = folium.MultiPolyLine(locations=mmi_lines, color='#DC143C', weight=6, opacity=1)
        map_bsl.add_children(mmi)
        
        mmi_lines=[mmi_contours.collections[1].get_paths()[0].vertices.tolist()]
        mmi = folium.MultiPolyLine(locations=mmi_lines, color='#DC143C', weight=6, opacity=1)
        map_bsl.add_children(mmi)
        
        mmi_lines=[mmi_contours.collections[2].get_paths()[0].vertices.tolist()]
        mmi = folium.MultiPolyLine(locations=mmi_lines, color='#DC143C', weight=6, opacity=1)
        map_bsl.add_children(mmi)
    
    #Aftershocks
    if plot_afters:
        for quake in afters:
            map_bsl.add_children(folium.CircleMarker(location=quake, radius=200, color='#191970',fill_color='#191970'))
    
    #Done, save map
    map_bsl.save(out_file)
