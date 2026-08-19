from pyuvdata import UVData
import numpy as np
import hdf5plugin
import copy
import hera_cal as hc
import uvtools as uvt
import matplotlib.pyplot as plt
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import argparse
from scipy import optimize


def gauss(x, amp, loc, scale):
    return amp * np.exp(-0.5 * (x-loc)**2 / scale**2)

def chisq(x0, x, y):
    yfit = gauss(x, *x0)
    return np.sum(np.abs(yfit - y)**2)

def gauss_fit(x0, x, y, method='powell'):
    fit = optimize.minimize(chisq, x0, args=(x, y), method=method)
    ypred = gauss(x, *fit.x)
    return fit, ypred
    
m=np.load("/net/jake/home/ntsikelelo/Simulated_data_files/m_slope_filter.npy")
path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
for lst in ["0h","1h","2h","3h"]:
    for spw in ["low"]:
        mode_array=np.array(["reduced", "full"])    
        for mode in mode_array:
        
            uvd = UVData()
            Model_file=""
            print(mode, lst, spw)
            if mode=="reduced":
                Model_file=path_data+"Full_LST_gleam_model_final_"+lst+"_"+spw+".uvh5"
            if mode=="full":
                Model_file=path_data+"Model_gleam_vis_expanded_"+lst+"_"+spw+".uvh5"
            
            print("determining sigma values")
            uvd.read(Model_file, read_data=False)
            Model = hc.frf.FRFilter(Model_file)
            pol=uvd.polarization_array[0]
            uvd.read(Model_file, times=Model.times, polarizations=pol)
            F_model = hc.frf.FRFilter(uvd)
            F_model.fft_data(ax='time', window='blackman', overwrite=True, ifft=True)
            antpos = F_model.antpos
            fr_select = (0< F_model.frates) & (F_model.frates < 5)
            fr_select_negative=(-5 < F_model.frates) & (F_model.frates<0)
            x0 = np.array([1e-3, 2.0, 0.3])
            
            # filter the data!
            Model_filt_data = copy.deepcopy(F_model.data)
            filter_center_all={}
            filter_half_width_all={}
            ## determine width using only gleam model
            
            for k in Model_filt_data:
            
                blvec = (antpos[k[1]] - antpos[k[0]])
                bl_len_EW = blvec[0]
                fringe_value=m*bl_len_EW
                if bl_len_EW>30:
                
                    y = np.abs(F_model.dfft[k]).mean(1)
                    x0[1]=fringe_value
                    x = F_model.frates[fr_select]
                    fit, ypred = gauss_fit(x0, x, y[fr_select],  method='powell')
                     # make the filter
                    gmean, gsigma = fit.x[1:]   
                    filter_center = gmean * 1e-3
                    filter_half_width_here = np.abs(gsigma) * 2 * 1e-3
                   
                    if filter_half_width_here< 0.008:
                        filter_center_all[k]=filter_center
                        filter_half_width_all[k]=filter_half_width_here
                    
                    
            
                if bl_len_EW<-30:
                
                    y = np.abs(F_model.dfft[k]).mean(1)
                    x0[1]=fringe_value
                    x = F_model.frates[fr_select_negative]
                    fit, ypred = gauss_fit(x0, x, y[fr_select_negative],  method='powell')
                     # make the filter
                    gmean, gsigma = fit.x[1:]   
                    filter_center=gmean* 1e-3
                    filter_half_width_here = np.abs(gsigma) * 2 * 1e-3
                    if filter_half_width_here< 0.008:
                        filter_center_all[k]=filter_center
                        filter_half_width_all[k]=filter_half_width_here
                        
              
                    
            if mode=="reduced":  
                np.save("filter_center_mainlobe_"+lst+"_"+spw+".npy",filter_center_all)
                np.save("filter_halfwidth_mainlobe_"+lst+"_"+spw+".npy", filter_half_width_all)
                print("filter_halfwidth_mainlobe_reduced_baseline saved")
            if mode=="full":
                np.save("filter_center_mainlobe_full_baseline_"+lst+"_"+spw+".npy",filter_center_all)
                np.save("filter_halfwidth_mainlobe_full_baseline_"+lst+"_"+spw+".npy", filter_half_width_all)
                print("filter_halfwidth_mainlobe_full_baseline saved")
    
     