from obspy.taup import taup_create,TauPyModel


tvel=u'/Users/dmelgar/FakeQuakes/M6_validation/structure/test.tvel'
out='/Users/dmelgar/FakeQuakes/M6_validation/structure/'

print 'Building model'

taup_create.build_taup_model(tvel,output_folder=out)

print 'Loading model'
velmod=TauPyModel(model='/Users/dmelgar/FakeQuakes/M6_validation/structure/test.npz',verbose=True)

print 'Getting phases'
#paths=velmod.get_ray_paths(10.1,1.0,phase_list=['S','s','Sn'])
#paths2=velmod.get_ray_paths(10.1,3.0,phase_list=['S','s','Sn'])
paths3=velmod.get_ray_paths(18,0.4,phase_list=['S','s','Sn','SmS','Sm'])

paths3.plot(plot_type='cartesian')
#paths.plot(plot_type='cartesian')
#paths2.plot(plot_type='cartesian')


#from obspy.taup import TauPyModel
#
#velmod=TauPyModel(model='iasp91')
#paths=velmod.get_ray_paths(10.1,3.0,phase_list=['p','P','Pn','S','s','Sn'])
#paths.plot(plot_type='cartesian')
