from pyuvdata import UVData
import numpy as np
import hdf5plugin
import copy
import hera_cal as hc
import uvtools as uvt
import matplotlib.pyplot as plt
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

filter_factor     = [1e-8]
m=np.load("/net/jake/home/ntsikelelo/Simulated_data_files/m_slope_filter.npy")

mode="model"
path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
if mode=="model":
    path_data="/net/sinatra/vault-ike/ntsikelelo/Model_Data_H6C/"
Model_complete_file = path_data+"single_pol_"+mode+"_combined_2h_test.uvh5"
Model_complete = hc.frf.FRFilter(Model_complete_file)
freqs = Model_complete.freqs/1e6
times = (Model_complete.times-Model_complete .times.min()) * 24 * 3600  # seconds

uvd = UVData()
uvd.read(Model_complete_file)
F_model = hc.frf.FRFilter(uvd)

F_model.fft_data(ax='time', window='blackman', overwrite=True, ifft=True)
antpos=F_model.antpos

fr_select = (0< F_model.frates) & (F_model.frates < 5)

fr_select_negative=(-5 < F_model.frates) & (F_model.frates<0)



x0 = np.array([1e-3, 2.0, 0.3])


x_negative = F_model.frates[fr_select_negative]
x = F_model.frates[fr_select]

x0 = np.array([1e-3, 2.0, 0.3])

# filter the data!
Model_complete_filt_data = copy.deepcopy(F_model.data)
for k in Model_complete_filt_data:
    blvec = (antpos[k[1]] - antpos[k[0]])
    bl_len_EW = blvec[0]
    fringe_value=m*bl_len_EW
    

    if fringe_value > 0.5:


            y = np.abs(F_model.dfft[k]).mean(1)
            x0[1]=fringe_value
            fit, ypred = gauss_fit(x0, x, y[fr_select],  method='powell')
             # make the filter
            gmean, gsigma = fit.x[1:]    
            filter_center = -gmean * 1e-3
            filter_half_width = np.abs(gsigma) * 2 * 1e-3


            C = uvt.dspec.dayenu_mat_inv(times, filter_center, filter_half_width, filter_factor, no_regularization=False)
            R = np.linalg.pinv(C, rcond=1e-10)

            Model_complete_filt_data[k] = Model_complete_filt_data[k] - R @ Model_complete_filt_data[k]
            


print("Model complete data filter with half_width = "+str(filter_half_width)+" and filter center = "+str(filter_center))    

F_model.write_data(Model_complete_filt_data,path_data+"Mainlobe_filtered_single_pol_"+mode+"_combined_2h_test.uvh5",overwrite=True)