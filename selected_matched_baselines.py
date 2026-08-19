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

t_total=200
path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
path_data_model="/net/sinatra/vault-ike/ntsikelelo/Model_Data_H6C/"
f=open('H6C_data_name_files.txt',"r")
lines=f.readlines()
result=[]
for x in lines:
	result.append(x)
f.close()

uvh5name=[]
for i in range (len(result)):
	uvh5name.append(result[i][0:26])
	
		
uvh5name=np.array(uvh5name)
print(len(uvh5name),uvh5name[0])		


for i in range (int(t_total/2)):
    print(i)
    data_file = path_data+uvh5name[i]

    model_file = path_data_model+"abscal_model_"+uvh5name[i]

    hd = hc.io.HERAData(data_file)
    hd_m = hc.io.HERAData(model_file)

    pol=hd.polarization_array
    hd.read(polarizations=pol[0])
    hd_m.read(polarizations=pol[0])

    antpos, ants = hd_m.get_ENU_antpos()
    antpos_d = dict(zip(ants, antpos))


    hd_m.inflate_by_redundancy(tol=1e-3) 
    bl_lens = {bl: np.linalg.norm(antpos_d[bl[1]]-antpos_d[bl[0]]) for bl in hd_m.get_antpairs()}   
    # downselect to maximum baseline length (to reduce computational load in calibration)
    bls = []
    for bl in list(hd_m.get_antpairs()):
        # max 40 m baseline cut (also cut autos)
        bl_len = np.linalg.norm(antpos_d[bl[1]] - antpos_d[bl[0]])
        if bl_len > 0:
            bls.append(bl[:2])

    hd_m.select(bls=bls)  

    model, _, _=hd_m.build_datacontainers()
    data, _, _=hd.build_datacontainers()  
    baseline_pair=[]

    for k in data:
        if k in model:
            baseline_pair.append((k[0],k[1]))

    hd_m.select(bls=baseline_pair)
    hd.select(bls=baseline_pair)        
    hd.write_uvh5(path_data+"selected_data"+uvh5name[i], clobber=True)
    hd_m.write_uvh5(path_data_model+"selected_abscal_model_"+uvh5name[i], clobber=True)         