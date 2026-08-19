from pyuvdata import UVData, UVBeam, utils as uvutils
import numpy as np
import hdf5plugin
import copy
import uvtools as uvt
import matplotlib.pyplot as plt
import hera_cal as hc
import healvis
from scipy.signal import windows
import astropy.io.fits as fits


def Combine_data(hd_full, path_data, uvh5name, N, freqs, mode):

    bls_data=[]
    data_temp=[]
    print("combining data")
    for n in range (N):
        i=2*n
        print(i)
        data_file = path_data+uvh5name[i]
        if mode=="model":
            data_file=path_data+"abscal_model_"+uvh5name[i]

        hd = hc.io.HERAData(data_file)

        pol=hd.polarization_array
        hd.read(polarizations=pol[0], frequencies=freqs)
        
        if mode=="model":
            antpos, ants = hd.get_ENU_antpos()
            antpos_d = dict(zip(ants, antpos))
            hd.inflate_by_redundancy(tol=1e-3) 
            bl_lens = {bl: np.linalg.norm(antpos_d[bl[1]]-antpos_d[bl[0]]) for bl in hd.get_antpairs()}   
            # downselect to maximum baseline length (to reduce computational load in calibration)
            bls = []
            for bl in list(hd.get_antpairs()):
            # max 40 m baseline cut (also cut autos)
                bl_len = np.linalg.norm(antpos_d[bl[1]] - antpos_d[bl[0]])
                if bl_len > 0:
                    bls.append(bl[:2])

            hd.select(bls=bls)  

        data, _, _=hd.build_datacontainers()  

        data_temp.append(data)

    data_all={}

    for key in data_temp[0].keys():

        data_per_bls=[]
        for t in range (N):
            data_per_bls.append(data_temp[t][key][0,:])

        data_all[key]=np.array(data_per_bls)      



    for key in hd_full.get_antpairs():
        
        if mode=="model":
            if key in hd.get_antpairs():
                bls=(key[0],key[1],'nn')
                ind = hd_full.antpair2ind(*key)
                data_all_index=np.ones(shape=(N,1,len(freqs),1),dtype=complex)
                data_all_index[:,0,:,0]=data_all[bls]
                hd_full.data_array[ind]=data_all_index  
        else:
            bls=(key[0],key[1],'nn')
            ind = hd_full.antpair2ind(*key)
            data_all_index=np.ones(shape=(N,1,len(freqs),1),dtype=complex)
            data_all_index[:,0,:,0]=data_all[bls]
            hd_full.data_array[ind]=data_all_index  


        
    if n==(N-1):
        for key in hd_full.get_antpairs():
            if key in hd.get_antpairs():
                bls=(key[0],key[1],'nn')
                bls_data.append(bls)
                

    print(len(bls_data))
    hd_full.select(bls_data)
    hd_full.write_uvh5(path_data+"single_pol_"+mode+"_combined_2h_test.uvh5", clobber=True)
    
    
    
    
    
    

N=200

path_data_model="/net/sinatra/vault-ike/ntsikelelo/Model_Data_H6C/"
path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"


f=open('H6C_data_name_files.txt',"r")
lines=f.readlines()
result=[]
for x in lines:
	result.append(x)
f.close()

uvh5name=[]
uvh5name_diff=[]
for i in range (len(result)):
	uvh5name.append(result[i][0:26])
	uvh5name_diff.append(result[i][0:18]+'diff.uvh5')
		
uvh5name=np.array(uvh5name)
print(len(uvh5name),uvh5name[0],uvh5name_diff[0],N)	


time_array=[]
for n in range (N):
    i=2*n
    data_file = path_data+uvh5name[i]
    uvd= UVData()
    uvd.read(data_file,read_data=False) 


    time_array.append(np.unique(uvd.time_array[0])) # seconds
    
    
time_array=np.array(time_array)[:,0] 


Model_complete_file = path_data+uvh5name[0]
uvd1 = UVData()
uvd1.read(Model_complete_file, read_data=False)  
times=time_array
freqs = uvd1.freq_array[np.where(uvd1.freq_array>50e6)] ## load 50MHz and above

# antpos=uvd1.antenna_positions
# ants=uvd1.antenna_numbers

# #create antenna file
# with open('HERA_test_layout.csv', 'w') as f:
#     f.write("Name\t Number\t BeamID\t E\t N\t U\n\n")
#     for i, ant in enumerate(ants):
#         f.write("HH{:d}\t {:d}\t 0\t {:8.4f}\t {:8.4f}\t {:8.4f}\n".format(ant, ant, *antpos[i]))

# #make redundant groups        

# uniq_bls = bls_unique=uvd1.get_antpairs()

# uvd = healvis.simulator.setup_uvdata(array_layout="HERA_test_layout.csv",
#                                 telescope_location=(-30.72152777777791, 21.428305555555557, 1073.0000000093132),
#                                 telescope_name='HERA', freq_array=freqs, time_array=times, no_autos=False,
#                                 pols=['xx'], make_full=True, bls=uniq_bls)

# uvd.write_uvh5(path_data+"data_combined_2h.uvh5", clobber=True)

# print("empty data file created")

# np.save(path_data+"time_array_JD.npy", time_array)

# Model_complete_file = path_data_model+"abscal_model_"+uvh5name[0]
# uvd1 = UVData()
# uvd1.read(Model_complete_file, read_data=False)  
# times=time_array
# freqs = uvd1.freq_array[np.where(uvd1.freq_array>50e6)] ## load 50MHz and above

# antpos=uvd1.antenna_positions
# ants=uvd1.antenna_numbers

# #create antenna file
# with open('HERA_test_layout_model.csv', 'w') as f:
#     f.write("Name\t Number\t BeamID\t E\t N\t U\n\n")
#     for i, ant in enumerate(ants):
#         f.write("HH{:d}\t {:d}\t 0\t {:8.4f}\t {:8.4f}\t {:8.4f}\n".format(ant, ant, *antpos[i]))

# #make redundant groups        

# uniq_bls = bls_unique=uvd1.get_antpairs()

# uvd = healvis.simulator.setup_uvdata(array_layout="HERA_test_layout_model.csv",
#                                 telescope_location=(-30.72152777777791, 21.428305555555557, 1073.0000000093132),
#                                 telescope_name='HERA', freq_array=freqs, time_array=times, no_autos=False,
#                                 pols=['xx'], make_full=True, bls=uniq_bls)

# uvd.write_uvh5(path_data+"model_combined_2h.uvh5", clobber=True)

# hd_f = hc.io.HERAData(path_data+"data_combined_2h.uvh5")
# pol=hd_f.polarization_array
# hd_f.read(polarizations=pol[0])


# hd_full = hc.io.HERAData(path_data+"model_combined_2h.uvh5")
# pol=hd_full.polarization_array
# hd_full.read(polarizations=pol[0])

hd_diff = hc.io.HERAData(path_data+"data_combined_2h.uvh5")
pol=hd_diff.polarization_array
hd_diff.read(polarizations=pol[0])
    
# Combine_data(hd_full=hd_f, path_data=path_data, uvh5name=uvh5name, N=N, freqs=freqs, mode="data")    
    
# Combine_data(hd_full=hd_full, path_data=path_data_model, uvh5name=uvh5name, N=N, freqs=freqs, mode="model")

Combine_data(hd_full=hd_diff, path_data=path_data, uvh5name=uvh5name_diff, N=N, freqs=freqs, mode="diff")