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
import warnings
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run RFI flag.")
    parser.add_argument("--filter_name", type=str, default="Notch_filter_40_mHz", help="Name of the filter to apply.")
    parser.add_argument( "--N", type=int, default=275, help="Number of files to calibrate.")
    parser.add_argument( "--lst", type=str, default="2h", help="lst center to be filtered")
    parser.add_argument( "--spw", type=str, default="low", help="spectral window to be filtered")
    parser.add_argument( "--mode", type=str, default="gleam", help="sky model used")
    
    args = parser.parse_args()
    lst=args.lst
    spw=args.spw
    mode=args.mode
    print(f"Filter name: {args.filter_name}")
    print(f"Number of files: {args.N}")
    filter_name_2=args.filter_name
    print(mode)
    ## load file names to calibrate
    N=args.N

    path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
    path_data2="/net/sinatra/vault2/ntsikelelo/Data_H6C/"
    def numerical_sort(value):
        # Find all digit groups and convert them to integers
        numbers = re.findall(r'\d+', value)
        return list(map(int, numbers)) if numbers else [0]
    
    uvh5name = sorted(glob.glob(path_data2+"*sum.uvh5"), key=numerical_sort)

    #make redundant groups model        
   
    Model_complete_file = uvh5name[0]
    uvd1 = UVData()
    uvd1.read("/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Full_LST_gleam_model_final_"+lst+"_"+spw+".uvh5", read_data=False)
    freqs = uvd1.freq_array*1e-6
    
    Model_complete_file = uvh5name[0]
    uvd12 = UVData()
    uvd12.read(Model_complete_file, read_data=False)  
    antpos=uvd12.antenna_positions
    ants=uvd12.antenna_numbers
    antpos_d = dict(zip(ants, antpos))
    
    hd = io.HERAData(uvh5name[0])
    reds = redcal.get_reds(hd.data_antpos, pols=['nn'], pol_mode='1pol')
    
    i=0
    red_num=[]
    for red in reds:
        
        for bl in red:
            blvec = (antpos[np.where(ants==bl[0])] - antpos[np.where(ants==bl[1])])
            
            if np.abs(np.linalg.norm(blvec)-14)<2:
                red_num.append(i)
        i+=1
    red_num_array=np.unique(red_num)
    print("max red group index ["+str(np.min(red_num_array))+","+str(np.max(red_num_array))+"]")
    
    
    def filter_data_reds(cal_data, reds_14_m):
        
        cal_data_filtered={}
        index_flag=np.where((freqs>=70.2 ) & (freqs<=70.3) )
        for bls in reds_14_m: 
            if  bls in cal_data:
            
                cal_data_flagged=np.array(cal_data[bls][0,:])
                cal_data_flagged[index_flag]=np.nan*(np.ones(cal_data_flagged[index_flag].shape))
    
                cal_data_filtered[bls]=cal_data_flagged
          
        
        return cal_data_filtered
        
    def extrac_cal_data_reds(file_name=path_data+"abs_calibrated_data_full_", uvh5name=uvh5name, i=2, reds_here=reds):
            
            
            cal_data = np.load(file_name+uvh5name[4*i][43:69]+"_0h_"+spw+".npy", allow_pickle=True).item()
            if lst=="1h":
                cal_data = np.load(file_name+uvh5name[4*i+1][43:69]+"_"+lst+"_"+spw+".npy", allow_pickle=True).item()
            if lst=="2h":    
                cal_data = np.load(file_name+uvh5name[4*i+2][43:69]+"_"+lst+"_"+spw+".npy", allow_pickle=True).item()
            if lst=="3h":
                cal_data = np.load(file_name+uvh5name[4*i+3][43:69]+"_"+lst+"_"+spw+".npy", allow_pickle=True).item()
            cal_data_filtered=filter_data_reds(cal_data, reds_here)
                                  
        
    
            
            return cal_data_filtered   

    
    def redundant_abs_individual(file_name=path_data+"abs_calibrated_data_full_", filter_name=filter_name_2, mode=mode, N=N):
      
            for ti in range (N):
                print("time number "+str(ti))
                data_red_all=[]
                data_red_all_bls=[]
                
                for red in red_num_array:
                    reds_14_m=[]
                    for k in reds[red]:          
                        reds_14_m.append(k)
                        
                    cal_data_time=extrac_cal_data_reds(file_name=file_name, uvh5name=uvh5name, i=ti, reds_here=reds_14_m)
                    print(len(cal_data_time.keys()))
                              
                    cal_data_dic=cal_data_time
                    cal_data=[]
                    for key in cal_data_dic:
                        cal_data.append(cal_data_dic[key])
                    data_mean=np.abs(np.mean(np.array(cal_data), axis=0))
                    data_no_nans=data_mean[~np.isnan(data_mean)]
                    freq_no_nans=freqs[~np.isnan(data_mean)]
    
                    
                    if len(freq_no_nans)>0:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", np.exceptions.RankWarning)
                        
                            coeffs = np.polyfit(freq_no_nans, data_no_nans, 20)
                        poly = np.poly1d(coeffs)
                        x_fit=freqs
                        y_fit=poly(x_fit)
                        sig=3
                        all_flags=[]
                        for bls in cal_data_dic: 
                            up=y_fit+sig*y_fit
                            down=y_fit/sig
                            index_flags=np.where((np.abs(cal_data_dic[bls])>up))
                            index_2 =np.where(np.abs(cal_data_dic[bls])<down)
                            data_flagged=np.array(cal_data_dic[bls])
                            data_flagged[index_flags]=np.nan
                            data_flagged[index_2]=np.nan
                            data_red_all.append(data_flagged)
                            data_red_all_bls.append(bls)
                            
    
                print(len(data_red_all_bls))
                data_red_all_dic={}
                bls_int=0
                for k in data_red_all_bls:
                    data_red_all_dic[k]=data_red_all[bls_int]
                    bls_int+=1
                np.save(path_data+"red_flagged_"+mode+"_"+filter_name_2+"_"+lst+"_"+spw+"_"+str(ti)+".npy",data_red_all_dic)
                # np.save(path_data+"test_bls_red_flagged_data_"+str(ti)+".npy",data_red_all_bls)
    
    
    redundant_abs_individual(file_name=path_data+"abs_calibrated_data_"+filter_name_2+"_"+mode+"_", filter_name=filter_name_2, mode=mode, N=N)

if __name__ == "__main__":
    main()        
        