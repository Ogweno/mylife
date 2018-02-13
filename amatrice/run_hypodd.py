import glob
from hypoddpy import HypoDDRelocator

working_dir=u'/Users/dmelgar/Amatrice2016/afters/hypodd/'
output_cross_correlation_file='/Users/dmelgar/Amatrice2016/afters/hypodd/output_files/crosscorr.txt'
relocator = HypoDDRelocator(working_dir=working_dir,
    cc_time_before=0.05,
    cc_time_after=0.2,
    
    cc_maxlag=0.1,
    cc_filter_min_freq=1.0,
    cc_filter_max_freq=20.0,
    cc_p_phase_weighting={"Z": 1.0},
    cc_s_phase_weighting={"Z": 1.0, "E": 1.0, "N": 1.0},
    cc_min_allowed_cross_corr_coeff=0.6)
    
    
relocator.add_event_files(u'/Users/dmelgar/Amatrice2016/afters/catalogs/with_picks/INGV_withpicks.0000.xml')
relocator.add_waveform_files(glob.glob("/Users/dmelgar/Amatrice2016/afters/SEED/*.mseed"))
relocator.add_station_files(u'/Users/dmelgar/Amatrice2016/afters/station_inventory.xml')

relocator.setup_velocity_model(
    model_type="layered_p_velocity_with_constant_vp_vs_ratio",
    layer_tops=[(-10000, 5.8)],
    vp_vs_ratio=1.73)


relocator.log("Starting relocator...")
relocator._parse_station_files()
relocator._write_station_input_file()
relocator._read_event_information()
relocator._write_ph2dt_inp_file()
relocator._create_event_id_map()
relocator._write_catalog_input_file()
relocator._compile_hypodd()
relocator._run_ph2dt()
relocator._parse_waveform_files()
relocator._cross_correlate_picks(outfile=output_cross_correlation_file)
relocator._write_hypoDD_inp_file()
relocator._run_hypodd()
relocator._create_output_event_file()
relocator._create_plots()

#relocator.start_relocation(output_event_file="relocated_events.xml")
