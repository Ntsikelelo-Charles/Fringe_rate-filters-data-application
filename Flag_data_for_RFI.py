import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import h5py
import hdf5plugin  # REQUIRED to have the compression plugins available
import numpy as np
from scipy import constants, interpolate
import copy
import glob
import re
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_rows', 1000)
from uvtools.plot import plot_antpos, plot_antclass
from hera_qm import ant_metrics, ant_class, xrfi
from hera_qm import ant_metrics, xrfi
from hera_cal import io, utils, redcal, apply_cal, datacontainer, abscal
from hera_filters import dspec
from IPython.display import display, HTML
import linsolve
import time
tstart = time.time()
import glob
import re
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run calibration script with custom parameters.")
    parser.add_argument( "--N", type=int, default=275, help="Number of files to calibrate.")
    parser.add_argument( "--lst", type=str, default="0h", help="lst center to be filtered")
    
    
    args = parser.parse_args()
    lsts_center=args.lst
    tstart = time.time()
    print(f"Number of files: {args.N}")

    
    ## load file names to calibrate
    N=args.N
    path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
    path_data2="/net/sinatra/vault2/ntsikelelo/Data_H6C/"
    
    def numerical_sort(value):
        # Find all digit groups and convert them to integers
        numbers = re.findall(r'\d+', value)
        return list(map(int, numbers)) if numbers else [0]
    
    uvh5name = sorted(glob.glob(path_data2+"*sum.uvh5"), key=numerical_sort)
    uvh5name_diff = sorted(glob.glob(path_data2+"*diff.uvh5"), key=numerical_sort)
    
    
    print("done loading data file names")
    
    
    # parse omnical settings
    OC_MAX_DIMS = int(os.environ.get("OC_MAX_DIMS", 4))
    OC_MIN_DIM_SIZE = int(os.environ.get("OC_MIN_DIM_SIZE", 8))
    OC_SKIP_OUTRIGGERS = os.environ.get("OC_SKIP_OUTRIGGERS", "FALSE").upper() == "TRUE"
    OC_MIN_BL_LEN = float(os.environ.get("OC_MIN_BL_LEN", 1))
    OC_MAX_BL_LEN = float(os.environ.get("OC_MAX_BL_LEN", 1e100))
    OC_MAXITER = int(os.environ.get("OC_MAXITER", 50))
    OC_MAX_RERUN = int(os.environ.get("OC_MAX_RERUN", 4))
    OC_RERUN_MAXITER = int(os.environ.get("OC_MAXITER", 25))
    OC_MAX_CHISQ_FLAGGING_DYNAMIC_RANGE = float(os.environ.get("OC_MAX_CHISQ_FLAGGING_DYNAMIC_RANGE", 1))
    OC_USE_PRIOR_SOL = os.environ.get("OC_USE_PRIOR_SOL", "FALSE").upper() == "TRUE"
    OC_PRIOR_SOL_FLAG_THRESH = float(os.environ.get("OC_PRIOR_SOL_FLAG_THRESH", .95))
    OC_USE_GPU = os.environ.get("SAVE_RESULTS", "FALSE").upper() == "TRUE"
    
    # parse RFI settings
    RFI_DPSS_HALFWIDTH = float(os.environ.get("RFI_DPSS_HALFWIDTH", 300e-9))
    RFI_NSIG = float(os.environ.get("RFI_NSIG", 4))
    
    
    
    
    #set flag parameters, same as standard hera calibration pipepline 
        
    
    # ant_metrics bounds for low correlation / dead antennas
    am_corr_bad = (0, float(os.environ.get("AM_CORR_BAD", 0.3)))
    am_corr_suspect = (float(os.environ.get("AM_CORR_BAD", 0.3)), float(os.environ.get("AM_CORR_SUSPECT", 0.5)))
    
    # ant_metrics bounds for cross-polarized antennas
    am_xpol_bad = (-1, float(os.environ.get("AM_XPOL_BAD", -0.1)))
    am_xpol_suspect = (float(os.environ.get("AM_XPOL_BAD", -0.1)), float(os.environ.get("AM_XPOL_SUSPECT", 0)))
    
    # bounds on solar altitude (in degrees)
    good_solar_altitude = (-90, float(os.environ.get("SUSPECT_SOLAR_ALTITUDE", 0)))
    suspect_solar_altitude = (float(os.environ.get("SUSPECT_SOLAR_ALTITUDE", 0)), 90)
    
    # bounds on zeros in spectra
    good_zeros_per_eo_spectrum = (0, int(os.environ.get("MAX_ZEROS_PER_EO_SPEC_GOOD", 2)))
    suspect_zeros_per_eo_spectrum = (0, int(os.environ.get("MAX_ZEROS_PER_EO_SPEC_SUSPECT", 8)))
    
    # bounds on autocorrelation power
    auto_power_good = (float(os.environ.get("AUTO_POWER_GOOD_LOW", 5)), float(os.environ.get("AUTO_POWER_GOOD_HIGH", 30)))
    auto_power_suspect = (float(os.environ.get("AUTO_POWER_SUSPECT_LOW", 1)), float(os.environ.get("AUTO_POWER_SUSPECT_HIGH", 60)))
    
    # bounds on autocorrelation slope
    auto_slope_good = (float(os.environ.get("AUTO_SLOPE_GOOD_LOW", -0.4)), float(os.environ.get("AUTO_SLOPE_GOOD_HIGH", 0.4)))
    auto_slope_suspect = (float(os.environ.get("AUTO_SLOPE_SUSPECT_LOW", -0.6)), float(os.environ.get("AUTO_SLOPE_SUSPECT_HIGH", 0.6)))
    
    # bounds on autocorrelation RFI
    auto_rfi_good = (0, float(os.environ.get("AUTO_RFI_GOOD", 1.5)))
    auto_rfi_suspect = (0, float(os.environ.get("AUTO_RFI_SUSPECT", 2)))
    
    # bounds on autocorrelation shape
    auto_shape_good = (0, float(os.environ.get("AUTO_SHAPE_GOOD", 0.1)))
    auto_shape_suspect = (0, float(os.environ.get("AUTO_SHAPE_SUSPECT", 0.2)))
    
    # bound on per-xengine non-noiselike power in diff
    bad_xengine_zcut = float(os.environ.get("BAD_XENGINE_ZCUT", 10.0))
    
    # bounds on chi^2 per antenna in omnical
    oc_cspa_good = (0, float(os.environ.get("OC_CSPA_GOOD", 2)))
    oc_cspa_suspect = (0, float(os.environ.get("OC_CSPA_SUSPECT", 3)))
    
    # print bounds
    for bound in ['am_corr_bad', 'am_corr_suspect', 'am_xpol_bad', 'am_xpol_suspect', 
                  'good_solar_altitude', 'suspect_solar_altitude',
                  'good_zeros_per_eo_spectrum', 'suspect_zeros_per_eo_spectrum',
                  'auto_power_good', 'auto_power_suspect', 'auto_slope_good', 'auto_slope_suspect',
                  'auto_rfi_good', 'auto_rfi_suspect', 'auto_shape_good', 'auto_shape_suspect',
                  'bad_xengine_zcut', 'oc_cspa_good', 'oc_cspa_suspect']:
        print(f'{bound} = {eval(bound)}')
    
    # print settings
    for setting in [ 'OC_MAX_DIMS', 'OC_MIN_DIM_SIZE', 'OC_SKIP_OUTRIGGERS', 
                'OC_MIN_BL_LEN', 'OC_MAX_BL_LEN', 'OC_MAXITER', 'OC_MAX_RERUN', 'OC_RERUN_MAXITER', 
                'OC_MAX_CHISQ_FLAGGING_DYNAMIC_RANGE', 'OC_USE_PRIOR_SOL', 'OC_PRIOR_SOL_FLAG_THRESH', 
                'OC_USE_GPU', 'RFI_DPSS_HALFWIDTH', 'RFI_NSIG']:
        print(f'{setting} = {eval(setting)}')
    
    
    
    n=0

    if lsts_center=="0h":
            n=0
       

    if lsts_center=="1h":
            n=1 
       

    if lsts_center=="2h":
            n=2

    if lsts_center=="3h":
            n=3    
            
    rfi_flags_all=[]  
    for x in range (N): 
        i=4*x+n
        print(i)
        read_start = time.time()
        
        # load data
        ##filter data when filter are used
        hd = io.HERAData(uvh5name[i])
        
        SUM_FILE=uvh5name[i]
      
        DIFF_FILE=SUM_FILE[0:58]+"diff.uvh5"
        hd_unfil = io.HERAData(SUM_FILE)
        hd_diff = io.HERAData(DIFF_FILE)
        times=hd_unfil.times[0:2]
        freq_range=hd_unfil.freq_array
        
        data_unfil, _, _ = hd_unfil.read(frequencies=freq_range, times=times)
        diff_data, _, _ = hd_diff.read(frequencies=freq_range, times=times)
        data, _, _ = hd.read(frequencies=freq_range,times=times)
        print(f'Finished loading data in {(time.time() - read_start) / 60:.2f} minutes.')
           
        
        # Antenna classification step 
        
        ants = sorted(set([ant for bl in hd.bls for ant in utils.split_bl(bl)]))
        auto_bls = [bl for bl in data if (bl[0] == bl[1]) and (utils.split_pol(bl[2])[0] == utils.split_pol(bl[2])[1])]
        antpols = sorted(set([ant[1] for ant in ants]))
        
        am = ant_metrics.AntennaMetrics(SUM_FILE, DIFF_FILE, sum_data=data_unfil, diff_data=diff_data)
        am.iterative_antenna_metrics_and_flagging(crossCut=am_xpol_bad[1], deadCut=am_corr_bad[1])
        am.all_metrics = {}  # this saves time and disk by getting rid of per-iteration information we never use
        
        # Turn ant metrics into classifications
        totally_dead_ants = [ant for ant, i in am.xants.items() if i == -1]
        am_totally_dead = ant_class.AntennaClassification(good=[ant for ant in ants if ant not in totally_dead_ants], bad=totally_dead_ants)
        am_corr = ant_class.antenna_bounds_checker(am.final_metrics['corr'], bad=[am_corr_bad], suspect=[am_corr_suspect], good=[(0, 1)])
        am_xpol = ant_class.antenna_bounds_checker(am.final_metrics['corrXPol'], bad=[am_xpol_bad], suspect=[am_xpol_suspect], good=[(-1, 1)])
        ant_metrics_class = am_totally_dead + am_corr + am_xpol
        if np.all([ant_metrics_class[utils.split_bl(bl)[0]] == 'bad' for bl in auto_bls]):
            raise ValueError('All antennas are flagged for ant_metrics.')
        
        min_sun_alt = np.min(utils.get_sun_alt(hd.times))
        solar_class = ant_class.antenna_bounds_checker({ant: min_sun_alt for ant in ants}, good=[good_solar_altitude], suspect=[suspect_solar_altitude])
        
        zeros_class = ant_class.even_odd_zeros_checker(data_unfil, diff_data, good=good_zeros_per_eo_spectrum, suspect=suspect_zeros_per_eo_spectrum)
        if np.all([zeros_class[utils.split_bl(bl)[0]] == 'bad' for bl in auto_bls]):
            raise ValueError('All antennas are flagged for too many even/odd zeros.')
            
        auto_power_class = ant_class.auto_power_checker(data_unfil, good=auto_power_good, suspect=auto_power_suspect)
        if np.all([(auto_power_class)[utils.split_bl(bl)[0]] == 'bad' for bl in auto_bls]):
            raise ValueError('All antennas are flagged for bad autocorrelation power/slope.')
        overall_class = auto_power_class  + zeros_class + ant_metrics_class + solar_class  
        
        print("done with antenna classification")
        
        ## RFI Flagging 
        
        def auto_bl_zscores(data, flag_array, cache={}):
            '''This function computes z-score arrays for each delay-filtered autocorrelation, normalized by the expected noise. 
            Flagged times/channels for the whole array are given 0 weight in filtering and are np.nan in the z-score.'''
            zscores = {}
            for bl in auto_bls:
                wgts = np.array(np.logical_not(flag_array), dtype=np.float64)
                model, _, _ = dspec.fourier_filter(data.freqs, data[bl], wgts, filter_centers=[0], filter_half_widths=[RFI_DPSS_HALFWIDTH], mode='dpss_solve',
                                                    suppression_factors=[1e-9], eigenval_cutoff=[1e-9], cache=cache)
                res = data[bl] - model
                int_time = 24 * 3600 * np.median(np.diff(data.times))
                chan_res = np.median(np.diff(data.freqs))
                int_count = int(int_time * chan_res)
                sigma = np.abs(model) / np.sqrt(int_count / 2)
                zscores[bl] = res / sigma    
                zscores[bl][flag_array] = np.nan
        
            return zscores
        
        def rfi_from_avg_autos(data, auto_bls_to_use, prior_flags=None, nsig=RFI_NSIG):
            '''Average together all baselines in auto_bls_to_use, then find an RFI mask by looking for outliers after DPSS filtering.'''
            
            # Compute int_count for all unflagged autocorrelations averaged together
            int_time = 24 * 3600 * np.median(np.diff(data.times_by_bl[auto_bls[0][0:2]]))
            chan_res = np.median(np.diff(data.freqs))
            int_count = int(int_time * chan_res) * len(auto_bls_to_use)
            avg_auto = {(-1, -1, 'ee'): np.mean([data[bl] for bl in auto_bls_to_use], axis=0)}
            
            # Flag RFI first with channel differences and then with DPSS
            antenna_flags, _ = xrfi.flag_autos(avg_auto, int_count=int_count, nsig=(RFI_NSIG * 5))
            if prior_flags is not None:
                antenna_flags[(-1, -1, 'ee')] = prior_flags
            _, rfi_flags = xrfi.flag_autos(avg_auto, int_count=int_count, flag_method='dpss_flagger',
                                           flags=antenna_flags, freqs=data.freqs, filter_centers=[0],
                                           filter_half_widths=[RFI_DPSS_HALFWIDTH], eigenval_cutoff=[1e-9], nsig=nsig)
        
            return rfi_flags    
        
        
        
        antenna_flags, array_flags = xrfi.flag_autos(data_unfil, flag_method="channel_diff_flagger", nsig=RFI_NSIG * 5, 
                                                     antenna_class=overall_class, flag_broadcast_thresh=.5)
        for key in antenna_flags:
            antenna_flags[key] = array_flags
        cache = {}
        _, array_flags = xrfi.flag_autos(data_unfil, freqs=data_unfil.freqs, flag_method="dpss_flagger",
                                         nsig=RFI_NSIG, antenna_class=overall_class,
                                         filter_centers=[0], filter_half_widths=[RFI_DPSS_HALFWIDTH],
                                         eigenval_cutoff=[1e-9], flags=antenna_flags, mode='dpss_matrix', 
                                         cache=cache, flag_broadcast_thresh=.5)
        
        xengine_diff_class = ant_class.non_noiselike_diff_by_xengine_checker(data_unfil, diff_data, flag_waterfall=array_flags, 
                                                                             antenna_class=overall_class, 
                                                                             xengine_chans=1, bad_xengine_zcut=bad_xengine_zcut)
        overall_class += xengine_diff_class
        if np.all([overall_class[utils.split_bl(bl)[0]] == 'bad' for bl in auto_bls]):
            raise ValueError('All antennas are flagged after flagging non-noiselike diffs.')
        
        
        # Iteratively develop RFI mask, excess RFI classification, and autocorrelation shape classification
        stage = 1
        rfi_flags = np.array(array_flags)
        prior_end_states = set()
        while True:
            print(stage)
            # compute DPSS-filtered z-scores with current array-wide RFI mask
            zscores = auto_bl_zscores(data_unfil, rfi_flags)
            rms = {bl: np.nanmean(zscores[bl]**2)**.5 if np.any(np.isfinite(zscores[bl])) else np.inf for bl in zscores}
            
            # figure out which autos to use for finding new set of flags
            candidate_autos = [bl for bl in auto_bls if overall_class[utils.split_bl(bl)[0]] != 'bad']
            if stage == 1:
                # use best half of the unflagged antennas
                med_rms = np.nanmedian([rms[bl] for bl in candidate_autos])
                autos_to_use = [bl for bl in candidate_autos if rms[bl] <= med_rms]
            elif stage == 2:
                # use all unflagged antennas which are auto RFI good, or the best half, whichever is larger
                med_rms = np.nanmedian([rms[bl] for bl in candidate_autos])
                best_half_autos = [bl for bl in candidate_autos if rms[bl] <= med_rms]
                good_autos = [bl for bl in candidate_autos if (overall_class[utils.split_bl(bl)[0]] != 'bad')
                              and (auto_rfi_class[utils.split_bl(bl)[0]] == 'good')]
                autos_to_use = (best_half_autos if len(best_half_autos) > len(good_autos) else good_autos)
            elif stage == 3:
                # use all unflagged antennas which are auto RFI good or suspect
                autos_to_use = [bl for bl in candidate_autos if (overall_class[utils.split_bl(bl)[0]] != 'bad')]
        
            # compute new RFI flags
            rfi_flags = rfi_from_avg_autos(data_unfil, autos_to_use)
        
            # perform auto shape and RFI classification
            overall_class = auto_power_class  + zeros_class + ant_metrics_class + solar_class + xengine_diff_class
            auto_rfi_class = ant_class.antenna_bounds_checker(rms, good=auto_rfi_good, suspect=auto_rfi_suspect, bad=(0, np.inf))
            overall_class += auto_rfi_class
            auto_shape_class = ant_class.auto_shape_checker(data_unfil, good=auto_shape_good, suspect=auto_shape_suspect,
                                                            flag_spectrum=np.sum(rfi_flags, axis=0).astype(bool), 
                                                            antenna_class=overall_class)
            overall_class += auto_shape_class
            
            # check for convergence by seeing whether we've previously gotten to this number of flagged antennas and channels
            if stage == 3:
                if (len(overall_class.bad_ants), np.sum(rfi_flags)) in prior_end_states:
                    break
                prior_end_states.add((len(overall_class.bad_ants), np.sum(rfi_flags)))
            else:
                stage += 1
        
        
        auto_class = auto_power_class + auto_rfi_class + auto_shape_class
        if np.all([overall_class[utils.split_bl(bl)[0]] == 'bad' for bl in auto_bls]):
            raise ValueError('All antennas are flagged after flagging for bad autos power/slope/rfi/shape.')        
        
        rfi_flags_all.append(rfi_flags)
        print("done with RFI flagging")
        
        
    combined_flags=[]
    for i in range (len(rfi_flags_all)):
        for ti in range (2):
            combined_flags.append(rfi_flags_all[i][ti,:])
    combined_flags_reduced = np.logical_or.reduce(np.array(combined_flags), axis=0)  
    np.save(path_data+"All_RFI_flags_"+lsts_center+".npy",rfi_flags_all) 
    np.save(path_data+"Reduced_combined_RFI_flags_"+lsts_center+".npy",combined_flags_reduced) 
    
        

if __name__ == "__main__":
    main()         
        
    
    print("done with file!")
