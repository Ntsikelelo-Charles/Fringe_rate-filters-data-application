from pyuvdata import UVData, UVBeam, utils as uvutils
import numpy as np
import hdf5plugin
import copy
import glob
import re




def Combine_data(input_data_path, output_data_path, uvh5name, N, mode, spw, lst):

    filename = uvh5name[0]

    n=0
    lsts_center=lsts_center_array[ls]
    if lsts_center=="0h":
            n=4
            filename = uvh5name[0]

    if lsts_center=="1h":
            n=5 
            filename = uvh5name[1]

    if lsts_center=="2h":
            n=6
            filename = uvh5name[2]

    if lsts_center=="3h":
            n=7
            filename = uvh5name[3]  
        
    uv1=UVData()
    uv1.read(filename, read_data=False)
    freqs = uv1.freq_array[np.where((uv1.freq_array>150e6) & (uv1.freq_array<175e6))] ## load 50MHz and above
    if spw=="low":
        freqs = uv1.freq_array[np.where((uv1.freq_array>50e6) & (uv1.freq_array<75e6))] ## load 50MHz and above
    uv1.read(filename)
    uv1.select(frequencies=freqs, times=uv1.time_array[0])
    uv3=copy.deepcopy(uv1)
    print("now doing "+spw+"_"+lst)
    
    
 
    for i in range (N):
        # print(4*i+n)  
        uv2=UVData()
        filename2 = uvh5name[4*i+n]
            
       
        uv2.read(filename2)
        uv2.select(frequencies=freqs, times=uv2.time_array[0])
        uv3 = uv3 + uv2

    
    uv3.write_uvh5(output_data_path+"single_pol_"+mode+"_combined_"+lst+"_"+spw+".uvh5",clobber=True)

output_data_path="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
input_data_path="/net/sinatra/vault2/ntsikelelo/Data_H6C/"

N=250


def numerical_sort(value):
    # Find all digit groups and convert them to integers
    numbers = re.findall(r'\d+', value)
    return list(map(int, numbers)) if numbers else [0]

uvh5name = sorted(glob.glob(input_data_path+"*sum.uvh5"), key=numerical_sort)
 
    
spw_array=np.array(["low","high"])
lsts_center_array=np.array([ "0h","1h","2h", "3h"])
# lsts_center_array=np.array(["3h"])
print("Combining data")
for fr in range(len(spw_array)):
    spw=spw_array[fr]
    for ls in range (len(lsts_center_array)):
        
        lsts_center=lsts_center_array[ls]
        Combine_data(input_data_path, output_data_path, uvh5name, N, mode="data", spw=spw, lst= lsts_center)



