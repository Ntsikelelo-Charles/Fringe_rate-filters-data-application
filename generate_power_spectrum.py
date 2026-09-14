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

Model_complete_file = uvh5name[0]
uvd1 = UVData()
uvd1.read(Model_complete_file, read_data=False)  

antpos=uvd1.antenna_positions
ants=uvd1.antenna_numbers
antpos_d = dict(zip(ants, antpos))
times_total=[]
lsts_total=[]
n=0
#get metadata 
for i in range (250):
   
    # print (lsts_center, 4*i+n)
    Model_complete_file =uvh5name[4*i+n]
    uvd1 = UVData()
    uvd1.read(Model_complete_file, read_data=False)  
    times=np.unique(uvd1.time_array)
    times_total.append(times[0])
    lsts=np.rad2deg(np.unique(uvd1.lst_array))*(24/360)
    lsts_total.append(lsts[0])


uvd1 = UVData()
uvd1.read("/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Full_LST_gleam_model_final_"+lst+"_"+spw+".uvh5", read_data=False)
freqs = uvd1.freq_array*1e-6    
N=249
baseline=29
hd = io.HERAData(uvh5name[0])
hd.read()
reds = redcal.get_reds(hd.data_antpos, pols=['nn'], pol_mode='1pol')


data_bls,_,_=hd.build_datacontainers()
bls=[]
antpos, ants = hd.get_ENU_antpos(pick_data_ants=True)
antpos_d = dict(zip(ants, antpos))
for bl in data_bls:
        blvec = (antpos[np.where(ants==bl[0])] - antpos[np.where(ants==bl[1])])
        
        if np.abs(np.linalg.norm(blvec)-baseline)<2:
            
             bls.append(bl)

def JytomK(Jy, f0=163e6):
    """Convert Jansky to mK

    Args:
        Jy (float): brightness in Jy
        f0 (float): frequency in Hz

    Returns:
        float: brightness temperature in mK
    """

    # wavelength in m
    c=3.0e8
    lam = c / f0

    # flux density in W/m^2/Hz
    S = Jy * 1e-26
    kb=1.38064852e-23
    # temperature in K
    T = S * lam ** 2 / (2 * kb)

    # rettourn temperature in mK
    return T * 1e3

def power_spectrum(V, baseline, f1, df, delays):
    window = windows.blackmanharris(delays.shape[0])
    comso=hp.conversions.Cosmo_Conversions(Om_l=0.7,H0=70,Om_c=0.3)
    c=3e8
   
    f0=1420.0e6
    lambda_0=c/f0
    lambda_1=c/f1

    z=(lambda_1-lambda_0)/lambda_0
    factor=comso.bl_to_kperp(z=z)
    theta=(lambda_1/baseline)
    Dc=comso.dRperp_dtheta(z)*theta
    Beff=np.sum(window*df)
    Delta_D=comso.dRpara_df(z)*Beff
    k_par=comso.tau_to_kpara(z=z)*delays*1e-9
    k_per=factor*baseline
    theta_d=np.rad2deg(theta)
    Ae=theta_d**2
    

    steradian=(180.0/np.pi)**2 #degrees per steradian
    Omega= Ae/steradian
    V_mk=JytomK(Jy=V,f0=f1)
    constant=((Dc**2*Delta_D)/(Beff))*(1/(Omega*Beff))
    P=V_mk**2*constant*df**2
    
    # print(" at redshift "+str(z))
    return k_par,np.array(P)

def interpolate_vis(vis=None, freqs=freqs):

        data=vis   
        data_time_no_nans=data[~np.isnan(data)]
        freq_no_nans=freqs[~np.isnan(data)]
        data_fitted=[]
        if len(freq_no_nans)>150:
            coeffs = np.polyfit(freq_no_nans, data_time_no_nans, 6)
            poly = np.poly1d(coeffs)
            x_fit=freqs
            y_fit=poly(x_fit)
            
            for i in range (data.shape[0]):
                if np.isnan(data[i]):
                    data_fitted.append(poly(x_fit[i]))
                else:
                    data_fitted.append(data[i])
        data_fitted=np.array(data_fitted)
        # print(data_fitted.shape)
        # plt.semilogy(x_fit, np.abs(data_fitted),'k*')
        return data_fitted

def vis_perform_delay_trans(freqs=freqs, vis1=None, vis2=None):

    tau=np.fft.fftshift(np.fft.fftfreq(len(freqs),np.abs(freqs[1]-freqs[0])))*1e3
    index=np.where(tau>0)
    tau_pos=tau[index]
    window=windows.blackmanharris(len(freqs))
    
    vis2_here=interpolate_vis(vis=vis2, freqs=freqs)
    vis1_here=interpolate_vis(vis=vis1, freqs=freqs)
    if len(vis1_here)>0 and len(vis2_here)>0: 
        g0_delay=np.fft.fftshift(np.fft.fft(window*vis1_here))
        g1_delay=np.fft.fftshift(np.fft.fft(window*vis2_here))
        cross_power=g0_delay*np.conjugate(g1_delay)
        # print(np.imag(cross_power), np.imag(g1_delay), np.imag(g1_delay))
    else:
        cross_power=np.zeros(tau.shape)
    kpar, power_spectrum_cross_power=power_spectrum(V=cross_power, baseline=14.5, f1=62.5e6, df=122e3, delays=tau)
    return kpar, power_spectrum_cross_power

def extrac_cal_data_post_RFI_filter(filter_name=filter_name_2, mode="gleam", lst="0h", spw="low", uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs):

    k=0 
    cal_all_data=[]
    kpar=np.ones(freqs.shape)
    for bl in bls_interest:
        
        for x in range(N):
            ti=k+x
                
        
                
            cal_data1 = np.load(path_data+"red_flagged_"+mode+"_"+filter_name+"_0h_"+spw+"_"+str(ti)+".npy", allow_pickle=True).item()
            cal_data2 = np.load(path_data+"red_flagged_"+mode+"_"+filter_name+"_1h_"+spw+"_"+str(ti)+".npy", allow_pickle=True).item()
            cal_data3 = np.load(path_data+"red_flagged_"+mode+"_"+filter_name+"_2h_"+spw+"_"+str(ti)+".npy", allow_pickle=True).item()
            cal_data4 = np.load(path_data+"red_flagged_"+mode+"_"+filter_name+"_3h_"+spw+"_"+str(ti)+".npy", allow_pickle=True).item()
            if bl in cal_data1 and bl in cal_data2:
                vis1=(filter_data(cal_data1[bl]))
                vis2=(filter_data(cal_data2[bl]))
                kpar, power_spectrum_cross_power=vis_perform_delay_trans(freqs=freqs, vis1=vis1, vis2=vis2)
                cal_all_data.append(power_spectrum_cross_power)

            if bl in cal_data3 and bl in cal_data4:  
                vis3=(filter_data(cal_data3[bl]))
                vis4=(filter_data(cal_data4[bl]))
                kpar, power_spectrum_cross_power2=vis_perform_delay_trans(freqs=freqs, vis1=vis3, vis2=vis4)
                cal_all_data.append(power_spectrum_cross_power2)

                    
                    
                    
                    
    

  
    return kpar, np.array(cal_all_data)   


def filter_data(cal_data):

    index_flag=np.where((freqs>=70.2 ) & (freqs<=70.3) )
    

    if len(np.array(cal_data).shape)==2:
        cal_data_flagged=np.array(cal_data)[0,:]
    else: 
        cal_data_flagged=np.array(cal_data)    
    cal_data_flagged[index_flag]=np.nan*(np.ones(cal_data_flagged[index_flag].shape))
    
    
    return cal_data_flagged
    
def extrac_cal_data(file_name=path_data+"abs_calibrated_data_full_", lst="0h", spw="low", uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs):

    k=0 
    cal_all_data=[]
    kpar=np.ones(freqs.shape)
    for bl in bls_interest:
        
        for x in range(N):
                n=k+x
                
        
                
                  
                cal_data1 = np.load(file_name+uvh5name[4*n][43:69]+"_"+lst+"_"+spw+".npy", allow_pickle=True).item()
                cal_data2 = np.load(file_name+uvh5name[4*n+1][43:69]+"_1h_"+spw+".npy", allow_pickle=True).item()
                cal_data3 = np.load(file_name+uvh5name[4*n+2][43:69]+"_2h_"+spw+".npy", allow_pickle=True).item()
                cal_data4 = np.load(file_name+uvh5name[4*n+3][43:69]+"_3h_"+spw+".npy", allow_pickle=True).item()
                if bl in cal_data1 and bl in cal_data2:
                    vis1=(filter_data(cal_data1[bl]))
                    vis2=(filter_data(cal_data2[bl]))
                    kpar, power_spectrum_cross_power=vis_perform_delay_trans(freqs=freqs, vis1=vis1, vis2=vis2)
                    cal_all_data.append(power_spectrum_cross_power)

                if bl in cal_data3 and bl in cal_data4:  
                    vis3=(filter_data(cal_data3[bl]))
                    vis4=(filter_data(cal_data4[bl]))
                    kpar, power_spectrum_cross_power2=vis_perform_delay_trans(freqs=freqs, vis1=vis3, vis2=vis4)
                    cal_all_data.append(power_spectrum_cross_power2)

    

  
    return kpar, np.array(cal_all_data)        



# kpar, cal_data_gleam_filter_notch=extrac_cal_data_post_RFI_filter(filter_name=filter_name_2, mode="gleam", lst="0h", spw="low", uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs)
# kpar, cal_data_full_filter_notch=extrac_cal_data_post_RFI_filter(filter_name=filter_name_2, mode="full", lst="0h", spw="low", uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs)
kpar, cal_data_gleam_filter_notch=extrac_cal_data(file_name=path_data+"abs_calibrated_data_"+filter_name+"_gleam_",lst=lst , spw=spw, uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs) 
kpar, cal_data_full_filter_notch=extrac_cal_data(file_name=path_data+"abs_calibrated_data_"+filter_name+"_full_",lst=lst , spw=spw, uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs) 
np.save("all_power_spectrum_cal_data_gleam_filter_notch_"+str(baseline)+".npy", cal_data_gleam_filter_notch)
np.save("all_power_spectrum_cal_data_gleam_filter_notch_"+str(baseline)+".npy", cal_data_full_filter_notch)
print("done notch filter")

kpar, cal_data_gleam=extrac_cal_data(file_name=path_data+"abs_calibrated_data_gleam_",lst=lst , spw=spw, uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs)
kpar, cal_data_full=extrac_cal_data(file_name=path_data+"abs_calibrated_data_full_",lst=lst , spw=spw, uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs)
np.save("all_power_spectrum_cal_data_full_"+str(baseline)+".npy", cal_data_full)
np.save("all_power_spectrum_cal_data_gleam_"+str(baseline)+".npy", cal_data_gleam)
np.save("kpar_values.npy",kpar)
print("done no filter case")


kpar, cal_data_gleam_filter=extrac_cal_data(file_name=path_data+"abs_calibrated_data_"+filter_name+"_gleam_",lst=lst , spw=spw, uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs) 
kpar, cal_data_full_filter=extrac_cal_data(file_name=path_data+"abs_calibrated_data_"+filter_name+"_full_",lst=lst , spw=spw, uvh5name=uvh5name, N=N, bls_interest=bls, freqs=freqs) 

np.save("all_power_spectrum_cal_data_gleam_filter_"+str(baseline)+".npy", cal_data_gleam_filter)
np.save("all_power_spectrum_cal_data_full_filter_"+str(baseline)+".npy", cal_data_full_filter)
print("done main lobe filter")

