from pyuvdata import UVData
import numpy as np
import hdf5plugin
import copy
import hera_cal as hc
import uvtools as uvt
import matplotlib.pyplot as plt

channel=1
total_times=20

filter_half_width=0.40e-3

filter_center=0.0

filter_factor= [1e-8]

path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
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



Model_complete_file = path_data+uvh5name[0]
Model_complete = hc.frf.FRFilter(Model_complete_file)
uvd = UVData()
freqs_all = Model_complete.freqs/1e6
freqs=freqs_all[np.where((np.min(freqs_all)<=freqs_all) & (freqs_all<=69))] *1e6
uvd.read(filename=Model_complete_file, frequencies=freqs)
total_chan=len(freqs)
times = (Model_complete.times-Model_complete .times.min()) * 24 * 3600  # seconds
uvd.conjugate_bls()
F = hc.frf.FRFilter(uvd,frequencies=freqs)
antpos = F.antpos


print("filter factor = "+str(filter_factor)+" and n_times="+str(total_times))
m=np.load("/net/jake/home/ntsikelelo/Simulated_data_files/m_slope_filter.npy")
c=np.load("/net/jake/home/ntsikelelo/Simulated_data_files/c_intercept_filter.npy")



baselines_files={}
antpos = F.antpos
test_data = F.data
i=0
for k in test_data:
    baselines_files[k]=i
    i+=1
    
np.save(path_data+"baseline_dic.npy",baselines_files)
baselines_dic=np.load(path_data+"baseline_dic.npy", allow_pickle=True).item()
total_bls=len(baselines_dic)

data=np.ones(shape=(total_times,total_bls,total_chan),dtype=complex)
time_array=np.ones(shape=(total_times)) 
    
for i in range (int(total_times/2)):
    n=2*i
    Model_complete_file = path_data+uvh5name[i]
    print(i)
    Model_complete = hc.frf.FRFilter(Model_complete_file)
    uvd = UVData()
    uvd.read(Model_complete_file, frequencies=freqs) 
    time_array[n:n+2]=Model_complete.times * 24 * 3600 # seconds
    uvd.conjugate_bls()
    F = hc.frf.FRFilter(uvd, frequencies=freqs)
    data_per_uvh5 = F.data 
    
    for k in data_per_uvh5:
        bl=baselines_dic[k]
        data[n:n+2,bl,:]=data_per_uvh5[k]
        bl+=1     

time_array=np.array(time_array) 
time=time_array-np.min(time_array)
    
    
C = uvt.dspec.dayenu_mat_inv(time, filter_center, filter_half_width, filter_factor, no_regularization=False)
R = np.linalg.pinv(C, rcond=1e-10)
data_filter=np.ones(shape=(total_times,total_bls,total_chan),dtype=complex)
for bl in range (total_bls):
    data_filter[:,bl,:] = R @ data[:,bl,:]
    

fringe_vis=np.fft.ifftshift(np.abs(np.fft.fft(np.abs(data[:,0,0]))))
frates=np.fft.ifftshift(np.fft.fftfreq(len(time),time[2]-time[1]))
fringe_vis_filter=np.fft.ifftshift(np.abs(np.fft.fft(np.abs(data_filter[:,0,0]))))


    
for i in range (int(total_times/2)):
    n=2*i
    Model_complete_file = path_data+uvh5name[i]
    Model_complete = hc.frf.FRFilter(Model_complete_file)
    uvd = UVData()
    uvd.read(Model_complete_file, frequencies=freqs)
    uvd.conjugate_bls()
    F = hc.frf.FRFilter(uvd, frequencies=freqs)
    data_per_uvh5 = F.data 
    
    for k in data_per_uvh5:
        blvec = (antpos[k[1]] - antpos[k[0]])
        bl_len_EW = blvec[0]
        if np.abs(bl_len_EW)>30:
            bl=baselines_dic[k]
            data_per_uvh5[k]=data_filter[n:n+2,bl,:]
            bl+=1     
        else:
           data_per_uvh5[k]=data[n:n+2,bl,:] 
    
    F.write_data(data_per_uvh5,path_data+"Filtered_notch_channel_"+str(channel)+"_"+uvh5name[i],overwrite=True)
    

np.save(path_data+"time_array_JD.npy",time_array)
np.save(path_data+"times.npy",time)
np.save(path_data+"Data_time.npy",data)
np.save(path_data+"Data_filtered_time.npy",data_filter) 

plt.figure(figsize=(10,15))
plt.plot(frates,fringe_vis_filter, label="filtered")
plt.plot(frates,fringe_vis, label="unfiltered")
plt.legend()
plt.xlabel("Frate [mHz]")
plt.ylabel("|V|")   
# plt.show()
plt.savefig("Fringe_rate_plot_Data")
