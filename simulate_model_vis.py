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
import pygdsm
import healpy

import glob
import re
path_data="/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
def numerical_sort(value):
    # Find all digit groups and convert them to integers
    numbers = re.findall(r'\d+', value)
    return list(map(int, numbers)) if numbers else [0]

uvh5name = sorted(glob.glob(path_data+"*sum.uvh5"), key=numerical_sort)

N=275
nside=64
gleam_flux_min=.5


Model_complete_file =uvh5name[0]
uvd1 = UVData()
uvd1.read(Model_complete_file, read_data=False)  
# freqs = uvd1.freq_array[np.where((uvd1.freq_array>150e6) & (uvd1.freq_array<175e6))] ## load 50MHz and above
freqs = uvd1.freq_array[np.where((uvd1.freq_array>50e6) & (uvd1.freq_array<75e6))] ## load 50MHz and above
print(len(freqs), freqs[0])


# freqs = uvd1.freq_array[np.where((uvd1.freq_array>50e6))] ## load 50MHz and above

#simulate diffuse emission
print("simulate diffuse emission")

# get GSM2016 model
GSM = pygdsm.GlobalSkyModel16()
theta, phi = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)))

gsm = np.array([healpy.get_interp_val(GSM.generate(f), theta, phi) for f in freqs/1e6])

gsm_sky = healvis.sky_model.SkyModel()
gsm_sky.set_data(gsm.T[None, :, :])
gsm_sky.Nfreqs = len(freqs)
gsm_sky.Nside = nside
gsm_sky.Nskies = 1
gsm_sky.Npix = gsm.shape[-1]
gsm_sky.freqs = freqs

gsm_sky.write_hdf5(path_data+'Low_gsm2016_nside'+str(nside)+'.hdf5', clobber=True)



# Compact sources GLEAM
print("simulate compact sources")
hdu = fits.open('GLEAM_EGC_v2.fits')
data = hdu[1].data

RA = data['RAJ2000']
RA[RA > 240] -= 360
DE = data['DEJ2000']
alpha = data['alpha']
Fint = data['int_flux_084']
Fp = data['peak_flux_084']
cut = ((Fp) > gleam_flux_min) & (RA < 120) & (RA > -90)

data = data[cut]
RA = RA[cut]
DE = DE[cut]
alpha = alpha[cut]
Fint = Fint[cut]
Fp = Fp[cut]


# fix nan spix
for i in range(len(RA)):
    if np.isfinite(alpha[i]): continue
    fstr = "Fp{:03d}"
    frq = np.array([122., 130., 143., 151., 158., 166., 174., 181., 189., 197., 204.])
    x = []
    xstr = []
    for f in frq:
        xst = fstr.format(int(f))
        if xst in data.dtype.fields:
            x.append(f)
            xstr.append(xst)
    x = np.asarray(x, dtype=np.float)
    y = np.log10([data[i][xs] for xs in xstr])
    if sum(~np.isnan(y)) < 2:
        # skip this source b/c all but 1 bins are negative or nan...
        continue
    spix = np.polyfit(np.log10(x)[~np.isnan(y)], y[~np.isnan(y)], deg=1)[0]

    # if this is unreasonable, set to 0
    if spix < -3 or spix > 1:
        spix = 0
    alpha[i] = spix
    
# load gleam files
select = (Fp > 0.1) & ~np.isnan(alpha)
RA, DE, Fp, alpha = RA[select], DE[select], Fp[select], alpha[select]
alpha[(alpha < -3) & (alpha > 1)] = 0    
    
# make dummy sky model object
gleam = healvis.sky_model.construct_skymodel('flat_spec', freqs=freqs, Nside=nside, sigma=1)
gleam.data[:] = 0    

# get indices for map
inds = healpy.ang2pix(nside, RA, DE, lonlat=True)

# insert sources
for i, ind in enumerate(inds):
    gleam.data[0, ind] = Fp[i] * (freqs / 184e6)**alpha[i]
    
# load A Team
anames, ara, adec, aF, afreq, aspix, _ = np.loadtxt('bright_sources.txt', dtype=str, skiprows=3).T
ara, adec, aF = ara.astype(float), adec.astype(float), aF.astype(float)
afreq, aspix = afreq.astype(float), aspix.astype(float)    

# get indices for map
inds = healpy.ang2pix(nside, ara, adec, lonlat=True)

# insert sources
for i, ind in enumerate(inds):
    gleam.data[0, ind, :] = aF[i] * (freqs / 1e6 / afreq[i])**aspix[i]
    
gleam.data *= healvis.utils.jy2Tsr(freqs, bm=healpy.nside2pixarea(nside))[None, None, :]

gleam.write_hdf5(path_data+'Low_gleam_nside'+str(nside)+'.hdf5', clobber=True)



# simulate visibilities 
print("simulating vis")
	
Model_complete_file = uvh5name[0]
uvd1 = UVData()
uvd1.read(Model_complete_file, read_data=False)  

antpos=uvd1.antenna_positions
ants=uvd1.antenna_numbers
antpos_d = dict(zip(ants, antpos))
times_total=[]
#get metadata 
for i in range (N):
    print("file number "+str(4*i))
    Model_complete_file =uvh5name[4*i]
    uvd1 = UVData()
    uvd1.read(Model_complete_file, read_data=False)  
    times=np.unique(uvd1.time_array)
    times_total.append(times[0])


times_total=np.array(times_total)

print(times_total.shape)

    
#create antenna file
with open('HERA_test_layout.csv', 'w') as f:
    f.write("Name\t Number\t BeamID\t E\t N\t U\n\n")
    for i, ant in enumerate(ants):
        f.write("HH{:d}\t {:d}\t 0\t {:8.4f}\t {:8.4f}\t {:8.4f}\n".format(ant, ant, *antpos[i]))

# make redundant groups        
reds = hc.redcal.get_pos_reds(antpos_d, include_autos=True)
uniq_bls = [red[0] for red in reds]

print("here now")
# uniq_bls=np.load("baselines_to_simulate.npy")[0]

uvd = healvis.simulator.setup_uvdata(array_layout="HERA_test_layout.csv",
                                telescope_location=(-30.72152777777791, 21.428305555555557, 1073.0000000093132),
                                telescope_name='HERA', freq_array=freqs, time_array=times_total, no_autos=False,
                                pols=['yy'], make_full=True, bls=uniq_bls)

uvd.write_uvh5(path_data+"Low_diffuse_model_final.uvh5", clobber=True)
uvd.write_uvh5(path_data+"Low_gleam_model_final.uvh5", clobber=True)
    
    
    
    
    
    
#apply primary beam to vis
# get data files
uvd_file_gleam = path_data+"Low_gleam_model_final.uvh5"
uvd_file_diffuse = path_data+"Low_diffuse_model_final.uvh5"
# get sky files
sky_file_diffuse = path_data+'Low_gsm2016_nside'+str(nside)+'.hdf5'
sky_file_gleam = path_data+'Low_gleam_nside'+str(nside)+'.hdf5'

# get beam files
# bfile1 = 'NF_HERA_Dipole_CCBeam_Port21_PowerBeam.fits'
uvb_file = path_data+"HERA-Beams/NicolasFagnoniBeams/NF_HERA_Dipole_efield_beam.fits"


interp_freqs = freqs
fchans = np.arange(interp_freqs.shape[0])

# load beam and interpolate to data frequnciesfreqs
beam = healvis.beam_model.PowerBeam(uvb_file)
beam.interpolation_function = 'az_za_simple'
beam.interp_freq(interp_freqs, inplace=True, kind='cubic')
# make jobs
jobs = []
for j in range(len(    for i in range(len(dfiles)):
        if j < 2: continue
        if i % 2 == 1: continue
        jobs.append([sfiles[i], dfiles[i], bfiles[i], fchans[j]])


healvis.simulator.run_simulation_partial_freq(fchans, uvd_file_gleam, sky_file_gleam, beam=beam,
                                                       fov=180, smooth_beam=False)
print("now simulating diffuse")

healvis.simulator.run_simulation_partial_freq(fchans, uvd_file_diffuse, sky_file_diffuse, beam=beam,
                                                       fov=180, smooth_beam=False)

print("all done!")
