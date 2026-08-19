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
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run calibration script with custom parameters.")
    parser.add_argument("--filter_name", type=str, default="Notch_filter_40_mHz", help="Name of the filter to apply.")
    parser.add_argument( "--N", type=int, default=275, help="Number of files to calibrate.")
    parser.add_argument( "--lst", type=str, default="2h", help="lst center to be filtered")
    parser.add_argument( "--spw", type=str, default="low", help="spectral window to be filtered")
    
    
    args = parser.parse_args()
    lst=args.lst
    spw=args.spw
    tstart = time.time()
    print(f"Filter name: {args.filter_name}")
    print(f"Number of files: {args.N}")
    filter_name=args.filter_name
    
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
    
   
    #breakup the filtered data
    hd_full = io.HERAData(path_data+"filtered_single_pol_combined_"+filter_name+"_data_"+lst+"_"+spw+".uvh5")
    
    n=0

    if lst=="0h":
            n=0
       

    if lst=="1h":
            n=1 
       

    if lst=="2h":
            n=2

    if lst=="3h":
            n=3
    print("doing lst "+lst)  
    for x in range (N): 
        i=4*x+n
        hd_temp = io.HERAData(uvh5name[i])
        hd_full.read(times=hd_temp.times[0])        
        data,_,_=hd_full.build_datacontainers()
        hd_temp.read(frequencies=hd_full.freqs, times=hd_temp.times[0])
        hd_temp.update(data=data)
        hd_temp.write_uvh5("/net/sinatra/vault-ike/ntsikelelo/filtered_data/"+filter_name+"_filtered_data_"+uvh5name[i][43:69], clobber=True)
    
    # for x in range (N): 
    #     i=4*x+n
        print(i)
        read_start = time.time()
        
        # load data
        ##filter data when filter are used
        hd = io.HERAData("/net/sinatra/vault-ike/ntsikelelo/filtered_data/"+filter_name+"_filtered_data_"+uvh5name[i][43:69])
        
        SUM_FILE=uvh5name[i]
        DIFF_FILE=uvh5name_diff[i]
        hd_unfil = io.HERAData(SUM_FILE)
        hd_diff = io.HERAData(DIFF_FILE)
        times=hd_unfil.times[0:2]
        freq_range=hd_unfil.freq_array ### might need to look this again
        
        data_unfil, _, _ = hd_unfil.read(frequencies=freq_range, times=times)
        diff_data, _, _ = hd_diff.read(frequencies=freq_range, times=times)
        data, _, _ = hd.read(times=times[0]) #not setting the requency range yet
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
        
        
        print("done with RFI flagging")
        
        
         ## Perform redundant calibration    
        
        def classify_off_grid(reds, all_ants):
            '''Returns AntennaClassification of all_ants where good ants are in reds while bad ants are not.'''
            ants_in_reds = set([ant for red in reds for bl in red for ant in utils.split_bl(bl)])
            on_grid = [ant for ant in all_ants if ant in ants_in_reds]
            off_grid = [ant for ant in all_ants if ant not in ants_in_reds]
            return ant_class.AntennaClassification(good=on_grid, bad=off_grid)
        
        def per_pol_filter_reds(reds, pols=['nn', 'ee'], **kwargs):
            '''Performs redcal filtering separately on polarizations (which might have different min_dim_size issues).'''
            return [red for pol in pols for red in redcal.filter_reds(copy.deepcopy(reds), pols=[pol], **kwargs)]  
        
        def check_if_whole_pol_flagged(redcal_class, pols=['Jee', 'Jnn']):
            '''Checks if an entire polarization is flagged. If it is, returns True and marks all antennas as bad in redcal_class.'''
            if np.logical_or(*[np.all([redcal_class[ant] == 'bad' for ant in redcal_class.ants if ant[1] == pol]) for pol in pols]):
                print('An entire polarization has been flagged. Stopping redcal.')
                for ant in redcal_class:
                    redcal_class[ant] = 'bad'
                return True
            return False    
        
        
        def recheck_chisq(cspa, sol, cutoff, avg_alg, rfi_flags):
            '''Recompute chisq per ant without apparently bad antennas to see if any antennas get better.'''
            avg_cspa = {ant: avg_alg(np.where(rfi_flags, np.nan, cspa[ant])) for ant in cspa}
            sol2 = redcal.RedSol(sol.reds, gains={ant: sol[ant] for ant in avg_cspa if avg_cspa[ant] <= cutoff}, vis=sol.vis)
            new_chisq_per_ant = {ant: np.array(cspa[ant]) for ant in sol2.gains}
            if len(set([bl[2] for red in per_pol_filter_reds(sol2.reds, ants=sol2.gains.keys(), antpos=hd.data_antpos, **fr_settings) for bl in red])) >= 2:
                redcal.expand_omni_gains(sol2, sol2.reds, data, chisq_per_ant=new_chisq_per_ant)
            for ant in avg_cspa:
                if ant in new_chisq_per_ant:
                    if np.any(np.isfinite(new_chisq_per_ant[ant])):
                        if not np.all(np.isclose(new_chisq_per_ant[ant], 0)):
                            new_avg_cspa = avg_alg(np.where(rfi_flags, np.nan, cspa[ant]))
                            if new_avg_cspa > 0:
                                avg_cspa[ant] = np.min([avg_cspa[ant], new_avg_cspa])
            return avg_cspa
        
        import ctypes
        def malloc_trim():
            try:
                ctypes.CDLL('libc.so.6').malloc_trim(0) 
            except OSError:
                pass
        
        
        
        # this enables better memory management on linux
        def redundant_calibration (data, redcal_class, overall_class, reds, rfi_flags):
        
            redcal_start = time.time()
            rc_settings = {'max_dims': OC_MAX_DIMS, 'oc_conv_crit': 1e-10, 'gain': 0.4, 'run_logcal': False,
                           'oc_maxiter': OC_MAXITER, 'check_after': OC_MAXITER, 'use_gpu': OC_USE_GPU}
            
            if check_if_whole_pol_flagged(redcal_class):
                # skip redcal, initialize empty sol and meta 
                sol = redcal.RedSol(reds)
                meta = {'chisq': None, 'chisq_per_ant': None}
            else:    
                    # perform first stage of redundant calibration 
                    meta, sol = redcal.redundantly_calibrate(data, reds, **rc_settings)
                
                    max_dly = np.max(np.abs(list(meta['fc_meta']['dlys'].values())))  # Needed for RFI delay-slope cal
                    median_cspa = recheck_chisq(meta['chisq_per_ant'], sol, oc_cspa_suspect[1] * 5, np.nanmedian, rfi_flags)
                
                     # remove particularly bad antennas (5x the bound on median, not mean)
                    cspa_class = ant_class.antenna_bounds_checker(median_cspa, good=(oc_cspa_good[0], oc_cspa_suspect[1] * 5), bad=[(-np.inf, np.inf)])
                    redcal_class += cspa_class
                    print(f'Removing {cspa_class.bad_ants} for >5x high median chi^2.')
                    for ant in cspa_class.bad_ants:
                        print(f'\t{ant}: {median_cspa[ant]:.3f}')
            
            
                    # iteratively rerun redundant calibration
            redcal_done = False
            rc_settings['oc_maxiter'] = rc_settings['check_after'] = OC_RERUN_MAXITER
            for i in range(OC_MAX_RERUN + 1-OC_MAX_RERUN ):
                # refilter reds and update classification to reflect new off-grid ants, if any
                reds = per_pol_filter_reds(reds, ex_ants=(overall_class + redcal_class).bad_ants, antpos=hd.data_antpos, **fr_settings)
                reds = sorted(reds, key=len, reverse=True)
                redcal_class += classify_off_grid(reds, ants)
                ants_in_reds = set([ant for red in reds for bl in red for ant in utils.split_bl(bl)])
                
                # check to see whether we're done
                if check_if_whole_pol_flagged(redcal_class) or redcal_done or (i == OC_MAX_RERUN):
                    break
            
                # re-run redundant calibration using previous solution, updating bad and suspicious antennas
                meta, sol = redcal.redundantly_calibrate(data, reds, sol0=sol, **rc_settings)
                malloc_trim()
                
                # recompute chi^2 for bad antennas without bad antennas to make sure they are actually bad
                mean_cspa = recheck_chisq(meta['chisq_per_ant'], sol, oc_cspa_suspect[1], np.nanmean, rfi_flags)
                
                # remove bad antennas
                cspa_class = ant_class.antenna_bounds_checker(mean_cspa, good=oc_cspa_good, suspect=oc_cspa_suspect, bad=[(-np.inf, np.inf)])
                for ant in cspa_class.bad_ants:
                    if mean_cspa[ant] < np.max(list(mean_cspa.values())) / OC_MAX_CHISQ_FLAGGING_DYNAMIC_RANGE:
                        cspa_class[ant] = 'suspect'  # reclassify as suspect if they are much better than the worst antennas
                redcal_class += cspa_class
                print(f'Removing {cspa_class.bad_ants} for high mean unflagged chi^2.')
                for ant in cspa_class.bad_ants:
                    print(f'\t{ant}: {mean_cspa[ant]:.3f}')
            
                if len(cspa_class.bad_ants) == 0:
                    redcal_done = True  # no new antennas to flag
            
            print(f'Finished redcal in {(time.time() - redcal_start) / 60:.2f} minutes.')
            overall_class += redcal_class
            
            expanded_reds = redcal.get_reds(hd.data_antpos, pols=['ee', 'nn'], pol_mode='2pol')
            expanded_reds = per_pol_filter_reds(expanded_reds, ex_ants=(ant_metrics_class + solar_class + zeros_class + auto_class + xengine_diff_class).bad_ants,
                                                max_dims=OC_MAX_DIMS, min_dim_size=OC_MIN_DIM_SIZE)
            if OC_SKIP_OUTRIGGERS:
                expanded_reds = redcal.filter_reds(expanded_reds, ex_ants=[ant for ant in ants if ant[0] >= 320])
            if len(sol.gains) > 0:
                redcal.expand_omni_vis(sol, expanded_reds, data, chisq=meta['chisq'], chisq_per_ant=meta['chisq_per_ant'])
            
            # now figure out flags, nsamples etc.
            omni_flags = {ant: (~np.isfinite(g)) | (ant in overall_class.bad_ants) for ant, g in sol.gains.items()}
            vissol_flags = datacontainer.RedDataContainer({bl: ~np.isfinite(v) for bl, v in sol.vis.items()}, reds=sol.vis.reds)
            single_nsamples_array = np.ones((len(data.times), len(data.freqs)), dtype=float)
            nsamples = datacontainer.DataContainer({bl: single_nsamples_array for bl in data})
            vissol_nsamples = redcal.count_redundant_nsamples(nsamples, [red for red in expanded_reds if red[0] in vissol_flags], 
                                                              good_ants=[ant for ant in overall_class if ant not in overall_class.bad_ants])
            for bl in vissol_flags:
                vissol_flags[bl][vissol_nsamples[bl] == 0] = True
            sol.make_sol_finite()
            
            return sol, meta, overall_class, vissol_flags
        
        
        
        # figure out and filter reds and classify antennas based on whether or not they are on the main grid
        fr_settings = {'max_dims': OC_MAX_DIMS, 'min_dim_size': OC_MIN_DIM_SIZE, 'min_bl_cut': OC_MIN_BL_LEN, 'max_bl_cut': OC_MAX_BL_LEN}
        reds = redcal.get_reds(hd.data_antpos, pols=['ee', 'nn'], pol_mode='2pol')
        reds = per_pol_filter_reds(reds, ex_ants=overall_class.bad_ants, antpos=hd.data_antpos, **fr_settings)
        if OC_SKIP_OUTRIGGERS:
            reds = redcal.filter_reds(reds, ex_ants=[ant for ant in ants if ant[0] >= 320])
        redcal_class = classify_off_grid(reds, ants) 
    
        
        indexfreq=np.where( (hd.freqs[0]<=freq_range) & (freq_range <=hd.freqs[-1]))[0]    
        rfi_flags_freq_select=rfi_flags[0,indexfreq[0]-1:indexfreq[-1]]    
        sol, meta, overall_class, vissol_flags=redundant_calibration (data, redcal_class, overall_class, reds, rfi_flags_freq_select)  
        
        np.save(path_data+"gains/gains_"+filter_name+"_redcal_"+lst+"_"+spw+"_"+uvh5name[i][43:69]+".npy",sol.gains)
        
        chis_redcal=np.where(rfi_flags_freq_select, np.nan, meta['chisq']['Jnn']).ravel()
        np.save(path_data+"chisq/chisq_"+filter_name+"_redcal_"+uvh5name[i][43:69]+".npy",chis_redcal)
        hd.update(data=sol.vis)
        hd.write_uvh5(path_data+"redundant_"+filter_name+"_cal_data_temp.uvh5", clobber=True)
        
        print("done with redundant calibration")
        
        ### Perform abscal NOTE that flags along freqs are not dealt with 
        
        
        def produce_model_data(mode="full_sky_model", freqs=freq_range, times=times):
            filter_name_here=filter_name[0:19]                                       
            hdm_gleam = io.HERAData(path_data+"Filtered_model_"+filter_name_here+"_gleam_vis_expanded_"+lst+"_"+spw+".uvh5")
            hdm_gleam.read(times=times, frequencies=freqs)
            model_data_gleam,_,_=hdm_gleam.build_datacontainers()
            
            if mode=="full_sky_model":
                hdm_diffuse = io.HERAData(path_data+"Filtered_model_"+filter_name_here+"_diffuse_vis_expanded_"+lst+"_"+spw+".uvh5")
                hdm_diffuse.read(times=times, frequencies=freqs)
                model_data_diffuse,_,_=hdm_diffuse.build_datacontainers()
            
            if mode=="gleam_only":
                model_data=copy.deepcopy(model_data_gleam)
            if mode=="full_sky_model":
                model_data=copy.deepcopy(model_data_gleam)
                for bl in model_data:
                    model_data[bl]=model_data[bl]+model_data_diffuse[bl]
            return model_data 
        hdm_gleam = io.HERAData("/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Filtered_model_Notch_filter_40_mHz_gleam_vis_expanded_"+lst+"_"+spw+".uvh5")
        freq_range_reduced=freq_range[np.where( (hdm_gleam.freqs[0]<=freq_range) & (freq_range <=hdm_gleam.freqs[-1]))]
        indexfreq=np.where( (hdm_gleam.freqs[0]<=freq_range) & (freq_range <=hdm_gleam.freqs[-1]))[0]
        model_data_expanded_data=produce_model_data(mode="gleam_only", freqs=freq_range_reduced, times=times[0])
        
        reds_here = redcal.get_reds(hd.data_antpos, pols=['nn'], pol_mode='1pol')
        reds_filtered = per_pol_filter_reds(reds_here, ex_ants=overall_class.bad_ants, antpos=hd.data_antpos, **fr_settings)
        all_bls=[]
        for rg in reds_filtered:
                    for bl in rg:
                        all_bls.append(bl)
        
        matched_bls=[]
        for key in all_bls:
            if key in model_data_expanded_data.keys():
                matched_bls.append(key)
        
        hd_cal= io.HERAData(path_data+"redundant_"+filter_name+"_cal_data_temp.uvh5")
        hd_cal.read(bls=matched_bls, frequencies=freq_range_reduced, times=times[0])
        redcal_data,_,_=hd_cal.build_datacontainers()
        
        
        hd_data= io.HERAData(SUM_FILE)
        hd_data.read(bls=matched_bls, frequencies=freq_range_reduced, times=times[0])
        data_reduced,_,_=hd_data.build_datacontainers()
        
        def abscal_post_red_cal(data, model_data, redcal_data, redcal_gains, rfi_flags, hd):
            antpos, ants = hd.get_ENU_antpos(pick_data_ants=True)
            antpos_d = dict(zip(ants, antpos))
            #choose only baselines in data in model
            model_data_modi=copy.deepcopy(redcal_data)
            for k in redcal_data.keys():
                    model_data_modi[k]=np.array(model_data[k])
            
            gain_cal_pol={}
            for k in redcal_gains:
                if k[1]=='Jnn':
                    gain_cal_pol[k]=redcal_gains[k]
            
            rc_flags = {k: np.zeros_like(gain_cal_pol[k], dtype=bool) for k in gain_cal_pol}
        
            #get waterfall flags
            # there is a per frequency slove so rfi should not be a problem
            for ant in rc_flags:
                rc_flags[ant]=rfi_flags
                   
        
                
            noise_wgts = {k: np.ones_like(redcal_data[k], dtype=float)  for k in redcal_data}
            n_cut=0
            n_include=0
            if filter_name=="Main_lobe_20_mHz" or filter_name=="Main_lobe_baseline_dependent":
                filter_center_all=np.load("filter_center_mainlobe_full_baseline.npy",allow_pickle=True).item()
                for k in noise_wgts:
                    if k not in filter_center_all:
                        noise_wgts[k] = np.ones(noise_wgts[k].shape)*1e-40
                        n_include+=1
                print("number of baseline exclude is "+str(n_include)) 
            else:    
                for k in noise_wgts:
            
                    blvec = (antpos[np.where(ants==k[0])][0] - antpos[np.where(ants==k[1])][0])
            
                    bl_len_EW = np.abs(blvec[0])
            
                    if bl_len_EW <30: 
                        n_cut+=1
            
                        noise_wgts[k] = np.ones(noise_wgts[k].shape)*1e-40
                print("number of baseline cuts in "+str(n_cut))        
                # calibration with lincal gains


            
            
            calibrated_data = copy.deepcopy(redcal_data)
            calibrated_data_filtered=copy.deepcopy(redcal_data)
            
            abscal_gains = abscal.post_redcal_abscal(model_data_modi, calibrated_data, noise_wgts, rc_flags, verbose=False, use_abs_amp_lincal=False)
            total_gains_log_cal = abscal.merge_gains([gain_cal_pol, abscal_gains])
            calibrated_data_final = copy.deepcopy(data)
            abscal.calibrate_in_place(calibrated_data_final, total_gains_log_cal)

            print("perfoming lincal filter")
            abscal_gains_lincal = abscal.post_redcal_abscal(model_data_modi, calibrated_data, noise_wgts, rc_flags, verbose=False, use_abs_amp_logcal=False)
            total_gains_lin_cal = abscal.merge_gains([total_gains_log_cal, abscal_gains_lincal])
            
            
            
            abscal.calibrate_in_place(calibrated_data_filtered, abscal_gains)
        
            residual_data_filtered=copy.deepcopy(redcal_data)
            for bl in residual_data_filtered:
                residual_data_filtered[bl]=calibrated_data_filtered[bl]-model_data_modi[bl]
        
            
            return calibrated_data_final,  total_gains_lin_cal, calibrated_data_filtered
            
    
        redcal_gains={}
        for key in sol.gains:
            redcal_gains[key]=sol.gains[key][0,:]
            
       
        model_data_full=produce_model_data(mode="full_sky_model", freqs=freq_range_reduced, times=times[0])
        calibrated_data_final_full, total_gains_full, calibrated_data_final_full_filtered=abscal_post_red_cal(data_reduced, model_data_full, redcal_data, redcal_gains, rfi_flags[0,indexfreq[0]-1:indexfreq[-1]], hd)

        calibrated_data_final_gleam, total_gains_gleam, calibrated_data_final_gleam_filtered=abscal_post_red_cal(data_reduced, model_data_expanded_data, redcal_data, redcal_gains, rfi_flags[0,indexfreq[0]-1:indexfreq[-1]], hd)

    
        
        ## FLagged data
        
        
        
        def flag_data_and_gains(calibrated_data_final, total_gains, rfi_flags):
        
            data_flagged_redcal=copy.deepcopy(calibrated_data_final)
            for bls in data_flagged_redcal:
               
                if not overall_class[utils.split_bl(bls)[0]]=='bad':
                    data_flagged_redcal[bls]=np.array(calibrated_data_final[bls])
                    data_flagged_redcal[bls][0,:][np.where(rfi_flags==True)]=np.nan
                        
                    
                else:
                   
                    data_flagged_redcal[bls]=np.ones(calibrated_data_final[bls].shape)*np.nan
            total_gains_flagged=copy.deepcopy(total_gains)
            for ant in total_gains_flagged:
                
                if not overall_class[utils.split_bl(bls)[0]]=='bad':
                    total_gains_flagged[ant][0,:][np.where(rfi_flags==True)]=np.nan
                            
                        
                else:
                       
                    total_gains_flagged[ant]=np.ones(total_gains_flagged[ant].shape)*np.nan            
            return data_flagged_redcal, total_gains_flagged
        
        data_flagged_redcal_gleam, total_gains_flagged_gleam=flag_data_and_gains(calibrated_data_final_gleam, total_gains_gleam, rfi_flags[0,indexfreq[0]-1:indexfreq[-1]])
        
        data_flagged_redcal_gleam_filtered, total_gains_flagged_gleam=flag_data_and_gains(calibrated_data_final_gleam_filtered, total_gains_gleam, rfi_flags[0,indexfreq[0]-1:indexfreq[-1]])
        
        data_flagged_redcal_full_filtered, total_gains_flagged_full=flag_data_and_gains(calibrated_data_final_full_filtered, total_gains_full, rfi_flags[0,indexfreq[0]-1:indexfreq[-1]])

        data_flagged_redcal_full, total_gains_flagged_full=flag_data_and_gains(calibrated_data_final_full, total_gains_full, rfi_flags[0,indexfreq[0]-1:indexfreq[-1]])
        
        ## save data: gains and calibrated vis
        
        np.save(path_data+"gains/gains_"+filter_name+"_gleam_flagged_abscal"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy",total_gains_flagged_gleam)
        np.save(path_data+"gains/gains_"+filter_name+"_full_flagged_abscal"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy",total_gains_flagged_full)
        
        np.save(path_data+"gains/gains_"+filter_name+"_gleam_abscal"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy",total_gains_gleam)
        np.save(path_data+"gains/gains_"+filter_name+"_full_abscal"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy",total_gains_full)
        
        bls=matched_bls
        # antpos, ants = hd.get_ENU_antpos(pick_data_ants=True)
        # antpos_d = dict(zip(ants, antpos))
        # for bl in data_flagged_redcal_full:
        #     blvec = (antpos[np.where(ants==bl[0])] - antpos[np.where(ants==bl[1])])
        #     if np.linalg.norm(blvec)-14<2:
        #         bls.append(bl)
        print("total of "+str(len(bls))+" 14 m baselines in calibration as opposed to "+str(len(reds[0]))) 
        
        cal_full={}
        for bl in bls:
            if bl in data_flagged_redcal_full:
                cal_full[bl]=data_flagged_redcal_full[bl]
        

        cal_full_filtered={}
        for bl in bls:
            if bl in data_flagged_redcal_full_filtered:
                 cal_full_filtered[bl]=data_flagged_redcal_full_filtered[bl]
    
        cal_gleam_filtered={}
        for bl in bls:
            if bl in data_flagged_redcal_gleam_filtered:
                cal_gleam_filtered[bl]=data_flagged_redcal_gleam_filtered[bl] 

        cal_gleam={}
        for bl in bls:
            if bl in data_flagged_redcal_gleam:
                cal_gleam[bl]=data_flagged_redcal_gleam[bl] 

           
                
        np.save(path_data+"abs_calibrated_data_"+filter_name+"_full_"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy", cal_full)
        np.save(path_data+"abs_calibrated_data_"+filter_name+"_full_filtered_"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy", cal_full_filtered)
        
        np.save(path_data+"abs_calibrated_data_"+filter_name+"_gleam_"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy", cal_gleam)
        np.save(path_data+"abs_calibrated_data_"+filter_name+"_gleam_filtered_"+uvh5name[i][43:69]+"_"+lst+"_"+spw+".npy", cal_gleam_filtered)
        print("done with file! "+filter_name+"_"+lst+"_"+spw)

if __name__ == "__main__":
    main()        
        
        
