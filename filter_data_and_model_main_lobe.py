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

def main():
    parser = argparse.ArgumentParser(description="Filter data with custom filter.")
    parser.add_argument("--filter_type", type=str, default="Notch", help="type of the filter to apply.")
    parser.add_argument( "--mode", type=str, default="data", help="filter data or model.")
    parser.add_argument( "--filter_half_width", type=float, default=0.25e-3, help="filter width")
    parser.add_argument( "--lst", type=str, default="2h", help="lst center to be filtered")
    parser.add_argument( "--spw", type=str, default="low", help="spectral window to be filtered")
    
    
    args = parser.parse_args()
    lst=args.lst
    spw=args.spw
    filter_half_width=args.filter_half_width
    mode=args.mode
    filter_type=args.filter_type

    print(f"Filter type: {filter_type}")
    print(f"Mode: {mode}")
    print(f"Filter half width: {filter_half_width}")


    # def filter_data_notch(filter_half_width=0.25e-3, mode="data"):
    
    filter_center=0.0
    m=np.load("/net/jake/home/ntsikelelo/Simulated_data_files/m_slope_filter.npy")
    path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
    filter_factor= [1e-3]
    print("filter factor "+str(filter_factor))
    Model_complete_file =""
    if mode=="data": 
        Model_complete_file = path_data+"single_pol_"+mode+"_combined_"+lst+"_"+spw+".uvh5"
    if mode=="gleam":
        Model_complete_file="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Full_LST_gleam_model_final_"+lst+"_"+spw+".uvh5"
    if mode=="diffuse":
        Model_complete_file="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Full_LST_diffuse_model_final_"+lst+"_"+spw+".uvh5"    
    
    # only run once explicitly to determine sigma values  
    generate_sigma_value=False
    redundant_model=False
    if generate_sigma_value and filter_type=="Main_lobe_baseline_dependent": 
        uvd = UVData()
        Model_file=""
        if redundant_model:
            Model_file=path_data+"Full_LST_gleam_model_final_"+lst+"_"+spw+".uvh5"
        else:
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
   
        print("Gaussian filtering data")
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
                    print( k, gsigma, fringe_value, filter_center, filter_half_width) 
                

            if fringe_value<-30:
            
                y = np.abs(F_model.dfft[k]).mean(1)
                x0[1]=fringe_value
                x = F_model.frates[fr_select_negative]
                fit, ypred = gauss_fit(x0, x, y[fr_select_negative],  method='powell')
                 # make the filter
                gmean, gsigma = fit.x[1:]   
                filter_center=gmean* 1e-3
                filter_half_width = np.abs(gsigma) * 2 * 1e-3
                if filter_half_width_here< 0.008:
                    filter_center_all[k]=filter_center
                    filter_half_width_all[k]=filter_half_width 
                    
                    print(fringe_value, k, filter_center,filter_half_width)
                
        if redundant_model:  
            np.save("filter_center_mainlobe.npy",filter_center_all)
            np.save("filter_halfwidth_mainlobe.npy", filter_half_width_all)
        else: 
            np.save("filter_center_mainlobe_full_baseline.npy",filter_center_all)
            np.save("filter_halfwidth_mainlobe_full_baseline.npy", filter_half_width_all)
            print("filter_halfwidth_mainlobe_full_baseline saved")
    
    uvd = UVData()
    uvd.read(Model_complete_file, read_data=False)
    Model_complete = hc.frf.FRFilter(Model_complete_file)
    pol=uvd.polarization_array[0]
    uvd.read(Model_complete_file, times=Model_complete.times, polarizations=pol)
    F_model = hc.frf.FRFilter(uvd)
    Model_complete_filt_data = copy.deepcopy(F_model.data)

        
    #notch filter
    times = (F_model.times-F_model.times.min()) * 24 * 3600  # seconds
    print("data loaded")
    
    # make covariance
    C = uvt.dspec.dayenu_mat_inv(times, filter_center, filter_half_width, filter_factor, no_regularization=False)
    
    antpos = F_model.antpos
    # take inverse to get filter matrix
    R = np.linalg.pinv(C, rcond=1e-10)
    if filter_type=="Notch_filter_40_mHz" or filter_type=="Notch_filter_80_mHz" or filter_type=="Notch_filter_100_mHz":
        print("Notch filtering data")
        for k in Model_complete_filt_data:
            blvec = (antpos[k[1]] - antpos[k[0]])
            bl_len_EW = blvec[0]
            fringe_value=m*bl_len_EW
            
        
            if np.abs(bl_len_EW) > 30:
                Model_complete_filt_data[k] = R @ Model_complete_filt_data[k]

                
    if filter_type=="Main_lobe_20_mHz":
        print("Gaussian filtering data")
        
        filter_center_all=np.load("filter_center_mainlobe.npy",allow_pickle=True).item()
        if mode=="data":
            filter_center_all=np.load("filter_center_mainlobe_full_baseline.npy",allow_pickle=True).item()
            print("change filter center full baseline")
        for k in Model_complete_filt_data:
            if k in filter_center_all:
                filter_center=filter_center_all[k]
                C = uvt.dspec.dayenu_mat_inv(times, filter_center, 2*filter_half_width, filter_factor, no_regularization=False)
                R = np.linalg.pinv(C, rcond=1e-10)
                Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]
   

    
    if filter_type=="Main_lobe_baseline_dependent" and mode=="data":
        filter_center_all=np.load("filter_center_mainlobe_full_baseline.npy",allow_pickle=True).item()
        filter_half_width_all=np.load("filter_halfwidth_mainlobe_full_baseline.npy",allow_pickle=True).item()
        print("Gaussian filtering data")
        for k in Model_complete_filt_data:

            if k in filter_center_all:
                filter_center=filter_center_all[k]
                filter_half_width = filter_half_width_all[k]
                C = uvt.dspec.dayenu_mat_inv(times, filter_center, 2*filter_half_width, filter_factor, no_regularization=False)
                R = np.linalg.pinv(C, rcond=1e-10)
        
                Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]
                
    if filter_type=="Main_lobe_baseline_dependent" and mode=="gleam":
        filter_center_all=np.load("filter_center_mainlobe.npy",allow_pickle=True).item()
        filter_half_width_all=np.load("filter_halfwidth_mainlobe.npy",allow_pickle=True).item()
        print("Gaussian filtering data")
        for k in Model_complete_filt_data:

            if k in filter_center_all:
    
                filter_center=filter_center_all[k]
                filter_half_width = filter_half_width_all[k]
                C = uvt.dspec.dayenu_mat_inv(times, filter_center, 2*filter_half_width, filter_factor, no_regularization=False)
                R = np.linalg.pinv(C, rcond=1e-10)
        
                Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]


    if filter_type=="Main_lobe_baseline_dependent" and mode=="diffuse":
        filter_center_all=np.load("filter_center_mainlobe.npy",allow_pickle=True).item()
        filter_half_width_all=np.load("filter_halfwidth_mainlobe.npy",allow_pickle=True).item()
        print("Gaussian filtering data")
        for k in Model_complete_filt_data:

            if k in filter_center_all:
                filter_center=filter_center_all[k]
                filter_half_width = filter_half_width_all[k]
                C = uvt.dspec.dayenu_mat_inv(times, filter_center, 2*filter_half_width, filter_factor, no_regularization=False)
                R = np.linalg.pinv(C, rcond=1e-10)
        
                Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]  

                
    print("Model complete data filter with half_width = "+str(filter_half_width)+" and filter center = "+str(filter_center))    

            
            
    F_model.write_data(Model_complete_filt_data,path_data+"filtered_single_pol_combined_"+filter_type+"_"+mode+"_"+lst+"_"+spw+".uvh5",overwrite=True)


if __name__ == "__main__":
    main()  


    
#filter_data_notch(filter_half_width=0.25e-3, mode="data", filter_type=="Notch")

# filter_data_notch(filter_half_width=0.25e-3, mode="model_gleam")
# filter_data_notch(filter_half_width=0.25e-3, mode="model_diffuse")