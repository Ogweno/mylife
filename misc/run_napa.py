import media_map

path_to_shakemap='/Users/dmelgar/BSL/media_map/napa/napa.shake.xyz'
path_to_aftershocks='/Users/dmelgar/code/GMT/Napa/afters.xy'
path_to_faults='/Users/dmelgar/code/GMT/Napa/UCERF3.xy'
out_file='/Users/dmelgar/BSL/media_map/napa_3levels.html'
epicenter=[38.215, -122.312]
zoom_start_level=9

media_map.run_map(out_file,epicenter,path_to_aftershocks,path_to_shakemap,path_to_faults,zoom_start_level,
            plot_shakemap=True,plot_faults=False)

