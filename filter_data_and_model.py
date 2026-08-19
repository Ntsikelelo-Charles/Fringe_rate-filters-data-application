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
from hera_cal import io, utils, redcal, apply_cal, datacontainer, abscal


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
    filter_factor= [1e-8]
    print("filter factor "+str(filter_factor))
    Model_complete_file =""
    if mode=="data": 
        Model_complete_file = path_data+"single_pol_"+mode+"_combined_"+lst+"_"+spw+".uvh5"
    if mode=="gleam":
        Model_complete_file="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Full_LST_gleam_model_final_"+lst+"_"+spw+".uvh5"
    if mode=="diffuse":
        Model_complete_file="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/Full_LST_diffuse_model_final_"+lst+"_"+spw+".uvh5"    
    
    
    uvd = UVData()
    uvd.read(Model_complete_file, read_data=False)
    freqs=uvd.freq_array
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
    # rfi_flags=np.load(path_data+"Reduced_combined_RFI_flags_0h.npy")[26:231]
    # all_

    # 200, 2, 1536
    hd_unfil = io.HERAData("/net/sinatra/vault2/ntsikelelo/Data_H6C/zen.2459861.25297.sum.uvh5")
    freqs_unfil=hd_unfil.freq_array
    index=np.where((freqs_unfil>=freqs[0]) & (freqs_unfil<=freqs[-1])) 
    flags=np.load(path_data+"All_RFI_flags_"+lst+".npy")
    All_RFI_flags=[]
    for ti in range (flags.shape[0]):
        flags_here=flags[ti,0,:]
        All_RFI_flags.append(flags_here[index])   
    All_RFI_flags=np.array(All_RFI_flags)               
    if filter_type=="Notch_filter_40_mHz" or filter_type=="Notch_filter_20_mHz" or filter_type=="Notch_filter_80_mHz":
        print("Notch filtering data")
        for k in Model_complete_filt_data:
            blvec = (antpos[k[1]] - antpos[k[0]])
            bl_len_EW = blvec[0]
            fringe_value=m*bl_len_EW
            if np.abs(bl_len_EW) > 30:
                Model_complete_filt_data[k] = R @ Model_complete_filt_data[k]
    if filter_type=="Main_lobe_20_mHz":            
        if mode=="diffuse" or mode=="gleam": 
            print("Gaussian filtering model data")
            filter_center_all=np.load("filter_center_mainlobe.npy",allow_pickle=True).item()
            for k in Model_complete_filt_data:
                if k in filter_center_all:
                    filter_center=filter_center_all[k]
                    C = uvt.dspec.dayenu_mat_inv(times, filter_center, 2*filter_half_width, filter_factor, no_regularization=False)
                    R = np.linalg.pinv(C, rcond=1e-10)
                    Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]
    if filter_type=="Main_lobe_baseline_dependent":            
        if  mode=="diffuse" or mode=="gleam": 
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
                
    if filter_type=="Main_lobe_20_mHz" and mode=="data":
        print("Gaussian filtering data")
        filter_center_all=np.load("filter_center_mainlobe_full_baseline.npy",allow_pickle=True).item()
        print("change filter center full baseline")
        for k in Model_complete_filt_data:
            if k in filter_center_all:
                Model_data_fitted=np.ones(Model_complete_filt_data[k].shape, dtype=complex)
                data= Model_complete_filt_data[k]
                for t in range(data.shape[0]-1):
                    rfi_flags=All_RFI_flags[t,:]
                    data[t,:][np.where(rfi_flags==True)]=np.nan  ##flag data with rfi    
                    data_time_no_nans=data[t,:][np.where(rfi_flags==False)]
                    freq_no_nans=freqs[np.where(rfi_flags==False)]
                    coeffs = np.polyfit(freq_no_nans, data_time_no_nans, 6)
                    poly = np.poly1d(coeffs)
                    x_fit=freqs
                    y_fit=poly(x_fit)
                    data_fitted=[]
                    for i in range (data[t,:].shape[0]):
                        if np.isnan(data[t,:][i]):
                            data_fitted.append(poly(x_fit[i]))
                        else:
                            data_fitted.append(data[t,:][i])
                data_fitted=np.array(data_fitted) 
                Model_data_fitted[t,:]=data_fitted                   
                filter_center=filter_center_all[k]
                C = uvt.dspec.dayenu_mat_inv(times, filter_center, 2*filter_half_width, filter_factor, no_regularization=False)
                R = np.linalg.pinv(C, rcond=1e-10)
                Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]
       

    
    if filter_type=="Main_lobe_baseline_dependent" and mode=="data":
        filter_center_all=np.load("filter_center_mainlobe_full_baseline.npy",allow_pickle=True).item()
        filter_half_width_all=np.load("filter_halfwidth_mainlobe_full_baseline.npy",allow_pickle=True).item()
        print("Gaussian filtering data")
        for k in Model_complete_filt_data:
            blvec = (antpos[k[1]] - antpos[k[0]])
            bl_len_EW = blvec[0]
            fringe_value=m*bl_len_EW
            if k in filter_center_all:
                Model_data_fitted=np.ones(Model_complete_filt_data[k].shape, dtype=complex)
                data= Model_complete_filt_data[k]
                for t in range(data.shape[0]-1):
                    rfi_flags=All_RFI_flags[t,:]
                    data[t,:][np.where(rfi_flags==True)]=np.nan  ##flag data with rfi    
                    data_time_no_nans=data[t,:][np.where(rfi_flags==False)]
                    freq_no_nans=freqs[np.where(rfi_flags==False)]
                    coeffs = np.polyfit(freq_no_nans, data_time_no_nans, 6)
                    poly = np.poly1d(coeffs)
                    x_fit=freqs
                    y_fit=poly(x_fit)
                    data_fitted=[]
                    for i in range (data[t,:].shape[0]):
                        if np.isnan(data[t,:][i]):
                            data_fitted.append(poly(x_fit[i]))
                        else:
                            data_fitted.append(data[t,:][i])
                    data_fitted=np.array(data_fitted) 
                    Model_data_fitted[t,:]=data_fitted    
            
                filter_center=filter_center_all[k]
                filter_half_width = filter_half_width_all[k]
                C = uvt.dspec.dayenu_mat_inv(times, filter_center, 2*filter_half_width, filter_factor, no_regularization=False)
                R = np.linalg.pinv(C, rcond=1e-10)
                Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]
                

                
   
    if filter_type=="Notch_filter_40_mHz_interpolated" or filter_type=="Notch_filter_20_mHz_interpolated" or filter_type=="Notch_filter_80_mHz_interpolated" and mode =="data":
        print("Notch filtering data interpolated")
        for k in Model_complete_filt_data:
            blvec = (antpos[k[1]] - antpos[k[0]])
            bl_len_EW = blvec[0]
            fringe_value=m*bl_len_EW
            if np.abs(bl_len_EW) > 30:
                Model_data_fitted=np.ones(Model_complete_filt_data[k].shape, dtype=complex)
                data= Model_complete_filt_data[k]
                for t in range(data.shape[0]-1):
                    rfi_flags=All_RFI_flags[t,:]
                    data[t,:][np.where(rfi_flags==True)]=np.nan  ##flag data with rfi    
                    data_time_no_nans=data[t,:][np.where(rfi_flags==False)]
                    freq_no_nans=freqs[np.where(rfi_flags==False)]
                    coeffs = np.polyfit(freq_no_nans, data_time_no_nans, 6)
                    poly = np.poly1d(coeffs)
                    x_fit=freqs
                    y_fit=poly(x_fit)
                    data_fitted=[]
                    for i in range (data[t,:].shape[0]):
                        if np.isnan(data[t,:][i]):
                            data_fitted.append(poly(x_fit[i]))
                  
                        else:
                            data_fitted.append(data[t,:][i])
                    data_fitted=np.array(data_fitted)
    
                    Model_data_fitted[t,:]=data_fitted       
                Model_complete_filt_data[k] = R @ Model_data_fitted
        
    print("Model complete data filter with half_width = "+str(filter_half_width)+" and filter center = "+str(filter_center))    
    F_model.write_data(Model_complete_filt_data,path_data+"filtered_single_pol_combined_"+filter_type+"_"+mode+"_"+lst+"_"+spw+".uvh5",overwrite=True, fix_autos=True)


if __name__ == "__main__":
    main()  


