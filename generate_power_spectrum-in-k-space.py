from pyuvdata import UVData, UVBeam, utils as uvutils
import numpy as np
import hdf5plugin
import copy
import uvtools as uvt
import matplotlib.pyplot as plt
import hera_cal as hc
# import healvis
from scipy.signal import windows

import astropy.io.fits as fits
import pandas as pd
import h5py
from hera_cal import io, utils, redcal, apply_cal, datacontainer, abscal
import glob
import os
from scipy import stats
import hera_pspec as hp
import astropy
from astropy.time import Time
from astropy.coordinates import get_sun
import astropy.units as u
from hera_pspec import conversions


import glob
import re
path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
path_data2="/net/sinatra/vault2/ntsikelelo/Data_H6C/"
def numerical_sort(value):
    # Find all digit groups and convert them to integers
    numbers = re.findall(r'\d+', value)
    return list(map(int, numbers)) if numbers else [0]

uvh5name = sorted(glob.glob(path_data2+"*sum.uvh5"), key=numerical_sort)

Model_complete_file = uvh5name[0]
uvd1 = UVData()
uvd1.read(Model_complete_file, read_data=False)  
pol=uvd1.polarization_array[0]
uvd1.read(Model_complete_file, polarizations=pol) 
antpos=uvd1.antenna_positions
ants=uvd1.antenna_numbers
antpos_d = dict(zip(ants, antpos))

#make redundant groups model        
reds = hc.redcal.get_pos_reds(antpos_d, include_autos=True)


filter_name="Main_lobe_20_mHz"
filter_name_2="Notch_filter_40_mHz_interpolated"
filter_name_3="Notch_filter_80_mHz"
filter_name_4="Main_lobe_baseline_dependent"
filter_name_5="Notch_filter_20_mHz"
lst="0h"
spw="low"

uvd1 = UVData()
uvd1.read("/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Full_LST_gleam_model_final_"+lst+"_"+spw+".uvh5", read_data=False)
freqs = uvd1.freq_array*1e-6    
N=30
k_int=50

hd = io.HERAData(uvh5name[0])
def interpolate_vis_data(vis=None, freqs=freqs):
        data_no_nans=vis[~np.isnan(vis)]
        freq_no_nans=freqs[~np.isnan(vis)]
        data_fitted=[]
        if len(freq_no_nans)>150:
            coeffs = np.polyfit(freq_no_nans, data_no_nans, 6)
            poly = np.poly1d(coeffs)
            x_fit=freqs
            y_fit=poly(x_fit)
            
            for f in range (vis.shape[0]):
                if np.isnan(vis[f]):
                    data_fitted.append(poly(x_fit[f]))
                else:
                    data_fitted.append(vis[f])
        data_fitted=np.array(data_fitted) 
        return data_fitted


def cylindrical_power_spectra(vis_file="", N=25,kint=0):
    antpos, ants = hd.get_ENU_antpos(pick_data_ants=True)
    antpos_d = dict(zip(ants, antpos))
    
    c = 3e8
    
    # --- WINDOW + DELAY TRANSFORM SETUP ---
    freqs = uvd1.freq_array
    df = np.median(np.diff(freqs))
    window = np.blackman(len(freqs))
    
    f1=62.5e6
    taus = np.fft.fftshift(np.fft.fftfreq(len(freqs), d=df))
    tau_pos=taus[np.where(taus>0)]
    
    f0=1420.0e6
    lambda_0=c/f0
    lambda_1=c/f1
    z=(lambda_1-lambda_0)/lambda_0
    
    # --- COSMOLOGY ---
    
    cosmo=hp.conversions.Cosmo_Conversions(Om_l=0.7,H0=70,Om_c=0.3)
    factor=cosmo.bl_to_kperp(z=z)
    Beff=np.sum(window*df)
    Delta_D=cosmo.dRpara_df(z)*Beff
    

    # --- ACCUMULATORS ---
    all_k_d = []
    all_ps_d = []
    
    # --- LOOP OVER BASELINES ---
    
    
    for pi in range (N):
        n=pi+0
        vis_data={}
        vis_temp=np.load(vis_file+uvh5name[4*n][43:69]+"_0h_"+spw+".npy", allow_pickle=True).item()
        vis_temp2=np.load(vis_file+uvh5name[4*n+1][43:69]+"_1h_"+spw+".npy", allow_pickle=True).item()
        for bls_temp in vis_temp:
            if bls_temp in vis_temp2: 
                data_temp=[]
                data_temp.append(vis_temp[bls_temp][0,:])
                data_temp.append(vis_temp2[bls_temp][0,:])
                vis_data[bls_temp]=np.array(data_temp)

        vis_temp3=np.load(vis_file+uvh5name[4*n+2][43:69]+"_2h_"+spw+".npy", allow_pickle=True).item()
        vis_temp4=np.load(vis_file+uvh5name[4*n+3][43:69]+"_3h_"+spw+".npy", allow_pickle=True).item()
        for bls_temp in vis_temp3:
            if bls_temp in vis_temp4: 
                data_temp=[]
                data_temp.append(vis_temp3[bls_temp][0,:])
                data_temp.append(vis_temp4[bls_temp][0,:])
                vis_data[bls_temp]=np.array(data_temp)
            
        for k in vis_data.keys():
            # print(k)
            i, j = k[0],k[1]
        
            vis1 = interpolate_vis_data(vis=vis_data[k][0,:], freqs=freqs)
            vis2 = interpolate_vis_data(vis=vis_data[k][1,:], freqs=freqs)

            
            
           
            vis_win = vis1 * window
            vis_win2 = vis2 * window
            delay_vis = np.fft.fftshift(np.fft.fft(vis_win))[np.where(taus>=0)]
            delay_vis2 = np.fft.fftshift(np.fft.fft(vis_win2))[np.where(taus>=0)]
            k_par=cosmo.tau_to_kpara(z=z)*tau_pos
        
            ps_un = delay_vis*np.conjugate(delay_vis2)/np.max(np.abs(delay_vis*np.conjugate(delay_vis2))) # (Ndelay,)
        
            # --- BASELINE VECTOR ---
            blvec = (antpos[np.where(ants == i)][0]- antpos[np.where(ants == j)][0])
            bl_len_m = np.linalg.norm(blvec)
            k_per=factor*bl_len_m    
            k_mag = np.sqrt(k_per**2+k_par**2)
            for i in range (len(k_mag)):
                k_per=factor*bl_len_m
                if np.abs(bl_len_m-14)<2 or np.abs(bl_len_m-29)<2:
                    k_mag = np.sqrt(k_per**2+k_par**2)
                    all_k_d.append(k_mag[i])
                    all_ps_d.append(ps_un[i])
         
        
    all_k = np.array(all_k_d)
    all_ps = np.array(all_ps_d)
    print(len(all_ps),len(all_k))    
        # --- STACK ---
        
        
    # # --- SPHERICAL (VOLUME) AVERAGING ---
    k_bins = np.linspace(all_k.min(), all_k.max(), len(k_par))
    
    Pk_ave= []
    counts = []
    sigma=[]
    
    for i in range(len(k_bins)-1):
        mask = (all_k >= k_bins[i]) & (all_k < k_bins[i+1])
        # print(k_bins[i], len(all_ps[mask]))
        if np.any(mask):
            Pk=[]
            for k in range (len(all_ps[mask])-1):
                
                Pk.append(all_ps[mask][k])
            Pk_ave.append(np.mean(np.array(Pk)))
            counts.append(len(all_k [mask]))
        else:
            counts.append(0)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
    return k_centers, np.sqrt(Pk_ave)


k_centers, Pk_ave=cylindrical_power_spectra(vis_file=path_data+"abs_calibrated_data_gleam_", N=N,kint=50)
k_centers, Pk_ave_full=cylindrical_power_spectra(vis_file=path_data+"abs_calibrated_data_full_", N=N, kint=50)
k_centers, Pk_ave_filter_notch=cylindrical_power_spectra(vis_file= path_data+"abs_calibrated_data_"+filter_name_2+"_gleam_", N=N, kint=50)
k_centers, Pk_ave_filter=cylindrical_power_spectra(vis_file= path_data+"abs_calibrated_data_"+filter_name+"_gleam_", N=N, kint=50)

np.save("all_baseline_power_spectrum_cal_data_gleam.npy",  Pk_ave)
np.save("all_baseline_power_spectrum_cal_data_full.npy",  Pk_ave_full)
np.save("all_baseline_power_spectrum_cal_data_gleam_filter.npy",  Pk_ave_filter)
np.save("all_baseline_power_spectrum_cal_data_gleam_filter_notch.npy",  Pk_ave_filter_notch)
np.save("k_centers.npy", k_centers)
print("done main lobe filter")



