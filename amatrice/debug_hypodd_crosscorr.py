import copy
import fnmatch
import glob
import json
import logging
import math
from obspy.core import read, Stream, UTCDateTime
from obspy.core.event import Catalog, Comment, Origin, read_events, \
    ResourceIdentifier
from obspy.signal import xcorrPickCorrection
from obspy.xseed import Parser
from obspy.core.inventory import read_inventory
import os
import progressbar
import shutil
import subprocess
import sys
import warnings
import glob
from hypoddpy import HypoDDRelocator


from hypodd_compiler import HypoDDCompiler





class HypoDDException(Exception):
    pass


# Set up relocator

working_dir=u'/Users/dmelgar/Amatrice2016/afters/hypodd/'
outfile='/Users/dmelgar/Amatrice2016/afters/hypodd/output_files/crosscorr.txt'
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

#relocator.setup_velocity_model(
#    model_type="layered_p_velocity_with_constant_vp_vs_ratio",
#    layer_tops=[(-10000, 5.8)],
#    vp_vs_ratio=1.73)
#
#
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


################################################################################


ct_file_path = os.path.join(relocator.paths["input_files"], "dt.cc")
if os.path.exists(ct_file_path):
    relocator.log("dt.cc input file already exists")
    #return
# This is by far the lengthiest operation and will be broken up in
# smaller steps
cc_dir = os.path.join(relocator.paths["working_files"], "cc_files")
if not os.path.exists(cc_dir):
    os.makedirs(cc_dir)
# Read the dt.ct file and get all event pairs.
dt_ct_path = os.path.join(relocator.paths["input_files"], "dt.ct")
if not os.path.exists(dt_ct_path):
    msg = "dt.ct does not exists. Did ph2dt run successfully?"
    raise HypoDDException(msg)
event_id_pairs = []
with open(dt_ct_path, "r") as open_file:
    for line in open_file:
        line = line.strip()
        if not line.startswith("#"):
            continue
        # Remove leading hashtag.
        line = line[1:]
        event_id_1, event_id_2 = map(int, line.split())
        event_id_pairs.append((event_id_1, event_id_2))
# Now for every event pair, calculate cross correlated differential
# travel times for every pick.
# Setup a progress bar.
for event_1, event_2 in event_id_pairs:
    # Update the progress bar.
    # filename for event_pair
    event_pair_file = os.path.join(cc_dir, "%i_%i.txt" %
                                    (event_1, event_2))
    current_pair_strings = []
    # Find the corresponding events.
    event_id_1 = relocator.event_map[event_1]
    event_id_2 = relocator.event_map[event_2]
    event_1_dict = event_2_dict = None
    for event in relocator.events:
        if event["event_id"] == event_id_1:
            event_1_dict = event
        if event["event_id"] == event_id_2:
            event_2_dict = event
        if event_1_dict is not None and event_2_dict is not None:
            break

    # Write the leading string in the dt.cc file.
    current_pair_strings.append(
        "# {event_id_1}  {event_id_2} 0.0".format(
            event_id_1=event_1, event_id_2=event_2))
    # Now try to cross-correlate as many picks as possible.
    for pick_1 in event_1_dict["picks"]:
        pick_1_station_id = pick_1["station_id"]
        pick_1_phase = pick_1["phase"]
        # Try to find the corresponding pick for the second event.
        pick_2 = None
        for pick in event_2_dict["picks"]:
            if pick["station_id"] == pick_1_station_id and pick["phase"] == pick_1_phase:
                pick_2 = pick
                break
        # No corresponding pick could be found.
        if pick_2 is None:
            continue

        else:
            station_id = pick_1["station_id"]
            # Try to find data for both picks.
            data_files_1 = relocator._find_data(station_id,
                                        pick_1["pick_time"] -
                                        relocator.cc_param["cc_time_before"],
                                        relocator.cc_param["cc_time_before"] +
                                        relocator.cc_param["cc_time_after"])
            data_files_2 = relocator._find_data(station_id,
                                        pick_2["pick_time"] -
                                        relocator.cc_param["cc_time_before"],
                                        relocator.cc_param["cc_time_before"] +
                                        relocator.cc_param["cc_time_after"])
            # If any pick has no data, skip this pick pair.
            if data_files_1 is False or data_files_2 is False:
                continue
            # Read all files.
            stream_1 = Stream()
            stream_2 = Stream()
            for waveform_file in data_files_1:
                stream_1 += read(waveform_file)
            for waveform_file in data_files_2:
                stream_2 += read(waveform_file)
            # Get the corresponing pick weighting dictionary.
            if pick_1_phase == "P":
                pick_weight_dict = relocator.cc_param[
                    "cc_p_phase_weighting"]
            elif pick_1_phase == "S":
                pick_weight_dict = relocator.cc_param[
                    "cc_s_phase_weighting"]
            all_cross_correlations = []
            # Loop over all picks and weight them.
            for channel, channel_weight in pick_weight_dict.iteritems():
                pick2_corr=None
                if channel_weight == 0.0:
                    continue
                # Filter the files to obtain the correct trace.
                network, station = station_id.split(".")
                st_1 = stream_1.select(network=network, station=station,
                                        channel="*%s" % channel)
                st_2 = stream_2.select(network=network, station=station,
                                        channel="*%s" % channel)
                max_starttime_st_1 = pick_1["pick_time"] - \
                    relocator.cc_param["cc_time_before"]
                min_endtime_st_1 = pick_1["pick_time"] + \
                    relocator.cc_param["cc_time_after"]
                max_starttime_st_2 = pick_2["pick_time"] - \
                    relocator.cc_param["cc_time_before"]
                min_endtime_st_2 = pick_2["pick_time"] + \
                    relocator.cc_param["cc_time_after"]
                # Attempt to find the correct trace.
                for trace in st_1:
                    if trace.stats.starttime > max_starttime_st_1 or \
                        trace.stats.endtime < min_endtime_st_1:
                        st_1.remove(trace)
                for trace in st_2:
                    if trace.stats.starttime > max_starttime_st_2 or \
                        trace.stats.endtime < min_endtime_st_2:
                        st_2.remove(trace)

                # cleanup merges, in case the event is included in
                # multiple traces (happens for events with very close
                # origin times)
                st_1.merge(-1)
                st_2.merge(-1)

                if len(st_1) > 1:
                    msg = "More than one matching trace found for {pick}"
                    relocator.log(msg.format(pick=str(pick_1)), level="warning")
                    relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = msg
                    continue
                elif len(st_1) == 0:
                    msg = "No matching trace found for {pick}"
                    relocator.log(msg.format(pick=str(pick_1)), level="warning")
                    relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = msg
                    continue
                trace_1 = st_1[0]

                if len(st_2) > 1:
                    msg = "More than one matching trace found for {pick}"
                    relocator.log(msg.format(pick=str(pick_1)), level="warning")
                    relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = msg
                    continue
                elif len(st_2) == 0:
                    msg = "No matching trace found for {pick}"
                    relocator.log(msg.format(pick=str(pick_1)), level="warning")
                    relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = msg
                    continue
                trace_2 = st_2[0]

                if trace_1.id != trace_2.id:
                    msg = "Non matching ids during cross correlation. "
                    msg += "(%s and %s)" % (trace_1.id, trace_2.id)
                    relocator.log(msg, level="warning")
                    relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = msg
                    continue
                if trace_1.stats.sampling_rate != \
                        trace_2.stats.sampling_rate:
                    msg = ("Non matching sampling rates during cross "
                            "correlation. ")
                    msg += "(%s and %s)" % (trace_1.id, trace_2.id)
                    relocator.log(msg, level="warning")
                    relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = msg
                    continue

                # Call the cross correlation function.
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    try:
                        pick2_corr, cross_corr_coeff = \
                            xcorrPickCorrection(
                                pick_1["pick_time"], trace_1,
                                pick_2["pick_time"], trace_2,
                                t_before=relocator.cc_param["cc_time_before"],
                                t_after=relocator.cc_param["cc_time_after"],
                                cc_maxlag=relocator.cc_param["cc_maxlag"],
                                filter="bandpass",
                                filter_options={
                                    "freqmin":
                                    relocator.cc_param["cc_filter_min_freq"],
                                    "freqmax":
                                    relocator.cc_param["cc_filter_max_freq"]},
                                plot=False)
                    except Exception, err:
                        # XXX: Maybe maxlag is too short?
                        if not err.message.startswith("Less than 3"):
                            msg = "Error during cross correlating: "
                            msg += err.message
                            relocator.log(msg, level="error")
                            relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = msg
                            continue
                    if pick2_corr!=None:
                        all_cross_correlations.append((pick2_corr,
                                                cross_corr_coeff, channel_weight))
            if len(all_cross_correlations) == 0:
                relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = "No cross correlations performed"
                continue
            # Now combine all of them based upon their weight.
            pick2_corr = sum([_i[0] * _i[2] for _i in
                                all_cross_correlations])
            cross_corr_coeff = sum([_i[1] * _i[2] for _i in
                                    all_cross_correlations])
            weight = sum([_i[2] for _i in all_cross_correlations])
            pick2_corr /= weight
            cross_corr_coeff /= weight
            relocator.cc_results.setdefault(pick_1['id'], {})[pick_2['id']] = (pick2_corr, cross_corr_coeff)
        # If the cross_corr_coeff is under the allowed limit, discard
        # it.
        if cross_corr_coeff < \
                relocator.cc_param["cc_min_allowed_cross_corr_coeff"]:
            continue
        # Otherwise calculate the corrected differential travel time.
        diff_travel_time = (pick_2["pick_time"] + pick2_corr -
            event_2_dict["origin_time"]) - (pick_1["pick_time"] -
            event_1_dict["origin_time"])
        string = "{station_id} {travel_time:.6f} {weight:.4f} {phase}"
        string = string.format(
            station_id=pick_1["station_id"],
            travel_time=diff_travel_time,
            weight=cross_corr_coeff,
            phase=pick_1["phase"])
        current_pair_strings.append(string)
    # Write the file.
    with open(event_pair_file, "w") as open_file:
        open_file.write("\n".join(current_pair_strings))
pbar.finish()
relocator.log("Finished calculating cross correlations.")
if outfile:
    relocator.save_cross_correlation_results(outfile)
# Assemble final file.
final_string = []
for cc_file in glob.iglob(os.path.join(cc_dir, "*.txt")):
    with open(cc_file, "r") as open_file:
        final_string.append(open_file.read().strip())
final_string = "\n".join(final_string)
with open(ct_file_path, "w") as open_file:
    open_file.write(final_string)