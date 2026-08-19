from pyuvdata import UVData, UVBeam, utils as uvutils
from hera_cal import io, utils, redcal, apply_cal, datacontainer, abscal
import hera_cal as hc
import time
import numpy as np
import glob
import re
import argparse

def main():
    parser = argparse.ArgumentParser(description="expand data.")
    parser.add_argument("--data_type", type=str, default="no_filter", help="type of data to expand.")
    parser.add_argument( "--spw", type=str, default="low", help="spectral window")
    parser.add_argument( "--lst", type=str, default="2h", help="lst center to be filtered")
    parser.add_argument( "--filter_name", type=str, default="Notch", help="filter name")
    
   

    
    
    args = parser.parse_args()
    lsts_center=args.lst
    filter_name=args.filter_name
    tstart = time.time()
    data_type=args.data_type
    spw=args.spw
    
    
    ## load file names to calibrate
    
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
    
    antpos=uvd1.antenna_positions
    ants=uvd1.antenna_numbers
    antpos_d = dict(zip(ants, antpos))
    
    
    #make redundant groups model        
    reds = hc.redcal.get_pos_reds(antpos_d, include_autos=True)
    
    #redundant group from data
    hd = io.HERAData(uvh5name[0])
    reds_here = redcal.get_reds(hd.data_antpos, pols=['nn'], pol_mode='1pol')
    
    all_bls=[]
    for rg in reds:
        for bl in rg:
            all_bls.append(bl)
    
    all_bls_model=[]
    for rg in reds_here:
                for bl in rg:
                    all_bls_model.append((bl[0],bl[1]))        
    
    match_baseline=[]
    for bl in all_bls_model:
        if bl in all_bls:
            bls=(bl[0],bl[1], 'nn')
            match_baseline.append(bls)
    
    if data_type=="single_time":
    
            mode_array=["diffuse","gleam"]
            for mode in mode_array:
                print(mode)
                hd_model= io.HERAData(path_data+"Single_time_gleam_model_final.uvh5")
                if mode=="diffuse":
                    hd_model= io.HERAData(path_data+"Single_time_diffuse_model_final.uvh5")
                hd_model.read()
                model_data,_,_=hd_model.build_datacontainers()
                ##expand model vis by redundancy 
                model_data_expanded={}
                for rg in reds:
                    for bl in rg:
                        bls=(bl[0],bl[1],'nn')
                        bls_model=(rg[0][0],rg[0][1], 'yy')
                        model_data_expanded[bls]=model_data[bls_model]
                
                hd_cal_model= io.HERAData(path_data+uvh5name[0])
                hd_cal_model.read(bls=match_baseline, times=hd_model.times, frequencies=hd_model.freqs)
                model_data_expanded_data,_,_=hd_cal_model.build_datacontainers()
                
                for bls in model_data_expanded_data:
                
                        if bls[0]<=bls[1]:
                    
                            model_data_expanded_data[bls]=model_data_expanded[bls]
                
                        else:   
                            model_data_expanded_data[bls]=model_data_expanded[(bls[1],bls[0],'nn')]
                ## save model vis
                hd_cal_model.update(data=model_data_expanded_data)
                
                hd_cal_model.write_uvh5(path_data+"Single_time_Model_"+mode+"_vis_expanded.uvh5", clobber=True) 
                print("this model "+mode+" has taken "+str((time.time()-tstart)/60.0)+" minutes to run")
        
    if data_type=="no_filter":            
        mode_array=["diffuse","gleam"]
        hd_cal_model= io.HERAData(path_data+"single_pol_data_combined_"+lsts_center+"_"+spw+".uvh5")
        # times_select=hd_cal_model.times
        print("do_no_filtered")
        # mode_array=["gleam"]
        for mode in mode_array:
            print(mode)
            hd_model= io.HERAData(path_data+"Full_LST_gleam_model_final_"+lsts_center+"_"+spw+".uvh5")
            if mode=="diffuse":
                hd_model= io.HERAData(path_data+"Full_LST_diffuse_model_final_"+lsts_center+"_"+spw+".uvh5")
         
            hd_model.read(times=hd_cal_model.times[0:250])
            model_data,_,_=hd_model.build_datacontainers()
            ##expand model vis by redundancy 
            model_data_expanded={}
            for rg in reds:
                for bl in rg:
                    bls=(bl[0],bl[1],'nn')
                    bls_model=(rg[0][0],rg[0][1], 'yy')
                    model_data_expanded[bls]=model_data[bls_model]
            
            hd_cal_model= io.HERAData(path_data+"single_pol_data_combined_"+lsts_center+"_"+spw+".uvh5")
            hd_cal_model.read(bls=match_baseline, times=hd_cal_model.times[0:250])
            model_data_expanded_data,_,_=hd_cal_model.build_datacontainers()
            
            for bls in model_data_expanded_data:
            
                    if bls[0]<=bls[1]:
                
                        model_data_expanded_data[bls]=model_data_expanded[bls]
       
                    else:   
                        model_data_expanded_data[bls]=model_data_expanded[(bls[1],bls[0],'nn')]
            ## save model vis
            hd_cal_model.update(data=model_data_expanded_data)
            
            hd_cal_model.write_uvh5(path_data+"Model_"+mode+"_vis_expanded_"+lsts_center+"_"+spw+".uvh5", clobber=True) 
            print("this model "+mode+" has taken "+str((time.time()-tstart)/60.0)+" minutes to run")
            
    if data_type=="filtered":
        mode_array=["gleam","diffuse"]
        # mode_array=["diffuse"]
        print("do_filtered")
        for mode in mode_array:
            print(mode)
            hd_cal_model= io.HERAData(path_data+"single_pol_data_combined_"+lsts_center+"_"+spw+".uvh5")
            hd_model= io.HERAData(path_data+"filtered_single_pol_combined_"+filter_name+"_"+mode+"_"+lsts_center+"_"+spw+".uvh5")
            hd_model.read(times=hd_cal_model.times[0:250])
            model_data,_,_=hd_model.build_datacontainers()
            ##expand model vis by redundancy 
            model_data_expanded={}
            for rg in reds:
                for bl in rg:
                    bls=(bl[0],bl[1],'nn')
                    bls_model=(rg[0][0],rg[0][1], 'yy')
                    model_data_expanded[bls]=model_data[bls_model]
           
            hd_cal_model= io.HERAData(path_data+"single_pol_data_combined_"+lsts_center+"_"+spw+".uvh5")
            hd_cal_model.read(bls=match_baseline, times=hd_cal_model.times[0:250])
            model_data_expanded_data,_,_=hd_cal_model.build_datacontainers()
            
            for bls in model_data_expanded_data:
            
                    if bls[0]<=bls[1]:
                
                        model_data_expanded_data[bls]=model_data_expanded[bls]
               
                    else:   
                        model_data_expanded_data[bls]=model_data_expanded[(bls[1],bls[0],'nn')]
            ## save model vis
            hd_cal_model.update(data=model_data_expanded_data)
            
            hd_cal_model.write_uvh5(path_data+"Filtered_model_"+filter_name+"_"+mode+"_vis_expanded_"+lsts_center+"_"+spw+".uvh5", clobber=True) 
            print("this model "+mode+" has taken "+str((time.time()-tstart)/60.0)+" minutes to run")   


if __name__ == "__main__":
    main()              