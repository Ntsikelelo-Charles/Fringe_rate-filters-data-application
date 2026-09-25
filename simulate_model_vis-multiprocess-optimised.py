import os
import copy
import glob
import re
import hdf5plugin
import numpy as np
import healpy
import pygdsm
import healvis
import hera_cal as hc
import multiprocess as mp
from astropy.io import fits
from pyuvdata import UVData
from functools import lru_cache

path_data = "/net/sinatra/vault-ike/ntsikelelo/Data_H6C/"
path_data2 = "/net/sinatra/vault2/ntsikelelo/Data_H6C/"

N = 250
nside = 64
gleam_flux_min = 0.5
spw_array = np.array(["mid", "high"])
lsts_center_array = np.array(["0h", "1h", "2h", "3h"])

BEAM_FILE = path_data + "HERA-Beams/NicolasFagnoniBeams/NF_HERA_Dipole_efield_beam.fits"
GLEAM_FITS = "/home/ntsikelelo/Projects/Fringe_rate_filtering/GLEAM_EGC_v2.fits"
A_TEAM_FILE = "/home/ntsikelelo/Projects/Fringe_rate_filtering/bright_sources.txt"
TELESCOPE_LOCATION = (-30.72152777777791, 21.428305555555557, 1073.0000000093132)

def numerical_sort(value):
    numbers = re.findall(r"\d+", value)
    return list(map(int, numbers)) if numbers else [0]

uvh5name = sorted(glob.glob(path_data2 + "*sum.uvh5"), key=numerical_sort)

def freq_selection(freqs, spw):
    freqs = np.asarray(freqs).ravel()
    if spw == "low":
        m = (freqs > 50e6) & (freqs < 75e6)
    elif spw == "mid":
        m = (freqs > 110e6) & (freqs < 135e6)
    else:
        m = (freqs > 150e6) & (freqs < 175e6)
    return freqs[m]

def load_static_catalogs():
    # Load once, outside loops
    hdu = fits.open(GLEAM_FITS, memmap=True)
    data = hdu[1].data

    RA = np.array(data["RAJ2000"], copy=True)
    RA[RA > 240] -= 360
    DE = np.array(data["DEJ2000"], copy=True)
    alpha = np.array(data["alpha"], copy=True)
    Fint = np.array(data["int_flux_084"], copy=True)
    Fp = np.array(data["peak_flux_084"], copy=True)

    a_names, ara, adec, aF, afreq, aspix, _ = np.loadtxt(
        A_TEAM_FILE, dtype=str, skiprows=3
    ).T
    ara = ara.astype(float)
    adec = adec.astype(float)
    aF = aF.astype(float)
    afreq = afreq.astype(float)
    aspix = aspix.astype(float)

    return data, RA, DE, alpha, Fint, Fp, ara, adec, aF, afreq, aspix

def build_gleam_sky(freqs, nside, gleam_data, RA, DE, alpha, Fint, Fp, ara, adec, aF, afreq, aspix):
    # Filter catalog once
    cut = (Fp > gleam_flux_min) & (RA < 120) & (RA > -90)
    RA = RA[cut]
    DE = DE[cut]
    alpha = alpha[cut]
    Fp = Fp[cut]

    # Fix NaN spectral indices only where needed
    nan_mask = ~np.isfinite(alpha)
    if np.any(nan_mask):
        frq = np.array([122., 130., 143., 151., 158., 166., 174., 181., 189., 197., 204.], dtype=float)
        valid_cols = []
        valid_freqs = []
        fstr = "Fp{:03d}"
        for f in frq:
            xst = fstr.format(int(f))
            if xst in gleam_data.dtype.fields:
                valid_cols.append(xst)
                valid_freqs.append(f)
        valid_freqs = np.asarray(valid_freqs, dtype=float)

        if len(valid_cols) >= 2:
            for i in np.where(nan_mask)[0]:
                y = np.log10([gleam_data[i][xs] for xs in valid_cols])
                ok = np.isfinite(y)
                if ok.sum() < 2:
                    continue
                spix = np.polyfit(np.log10(valid_freqs[ok]), y[ok], deg=1)[0]
                if spix < -3 or spix > 1:
                    spix = 0.0
                alpha[i] = spix

    select = (Fp > 0.1) & np.isfinite(alpha)
    RA, DE, Fp, alpha = RA[select], DE[select], Fp[select], alpha[select]

    # Vectorized map insertion
    gleam = healvis.sky_model.construct_skymodel("flat_spec", freqs=freqs, Nside=nside, sigma=1)
    gleam.data.fill(0)

    inds = healpy.ang2pix(nside, RA, DE, lonlat=True)
    flux = Fp[:, None] * (freqs[None, :] / 184e6) ** alpha[:, None]
    np.add.at(gleam.data[0], inds, flux)

    a_inds = healpy.ang2pix(nside, ara, adec, lonlat=True)
    a_flux = aF[:, None] * (freqs[None, :] / 1e6 / afreq[:, None]) ** aspix[:, None]
    np.add.at(gleam.data[0], a_inds, a_flux)

    gleam.data *= healvis.utils.jy2Tsr(freqs, bm=healpy.nside2pixarea(nside))[None, None, :]
    gleam.data = gleam.data.astype(np.float32, copy=False)
    return gleam

def build_gsm_sky(freqs, nside):
    GSM = pygdsm.GlobalSkyModel16()
    theta, phi = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)))

    gsm = np.empty((len(freqs), theta.size), dtype=np.float32)
    for i, f in enumerate(freqs / 1e6):
        gsm[i] = healpy.get_interp_val(GSM.generate(f), theta, phi).astype(np.float32, copy=False)

    gsm_sky = healvis.sky_model.SkyModel()
    gsm_sky.set_data(gsm.T[None, :, :])
    gsm_sky.Nfreqs = len(freqs)
    gsm_sky.Nside = nside
    gsm_sky.Nskies = 1
    gsm_sky.Npix = gsm.shape[-1]
    gsm_sky.freqs = freqs
    return gsm_sky

def write_antenna_layout_csv(ants, antpos, outname="HERA_test_layout.csv"):
    with open(outname, "w") as f:
        f.write("Name\t Number\t BeamID\t E\t N\t U\n\n")
        for ant, pos in zip(ants, antpos):
            f.write("HH{:d}\t {:d}\t 0\t {:8.4f}\t {:8.4f}\t {:8.4f}\n".format(ant, ant, *pos))

def get_times_total(uvh5_files):
    times_total = []
    for fname in uvh5_files:
        uvd = UVData()
        uvd.read(fname, read_data=False)
        times_total.append(np.unique(uvd.time_array)[0])
    return np.asarray(times_total)

def make_uvdata_template(freqs, times_total, ants, antpos, outname):
    write_antenna_layout_csv(ants, antpos)
    antpos_d = dict(zip(ants, antpos))
    reds = hc.redcal.get_pos_reds(antpos_d, include_autos=True)
    uniq_bls = [red[0] for red in reds]

    uvd = healvis.simulator.setup_uvdata(
        array_layout="HERA_test_layout.csv",
        telescope_location=TELESCOPE_LOCATION,
        telescope_name="HERA",
        freq_array=freqs,
        time_array=times_total,
        no_autos=False,
        pols=["yy"],
        make_full=True,
        bls=uniq_bls,
    )
    uvd.write_uvh5(outname, clobber=True)

# Load static catalog data once
gleam_data, RA, DE, alpha, Fint, Fp, ara, adec, aF, afreq, aspix = load_static_catalogs()

# Preload one metadata file for antenna positions
uvd0 = UVData()
uvd0.read(uvh5name[0], read_data=False)
antpos = uvd0.antenna_positions
ants = uvd0.antenna_numbers

# Main sky generation
for spw in spw_array:
    base_uv = UVData()
    base_uv.read(uvh5name[0], read_data=False)
    freqs = freq_selection(base_uv.freq_array, spw)

    print("simulate diffuse emission")
    gsm_sky = build_gsm_sky(freqs, nside)
    gsm_sky.write_hdf5(path_data + f"gsm2016_nside{nside}_{spw}.hdf5", clobber=True)

    print("simulate compact sources")
    gleam_sky = build_gleam_sky(freqs, nside, gleam_data, RA, DE, alpha.copy(), Fint, Fp, ara, adec, aF, afreq, aspix)
    gleam_sky.write_hdf5(path_data + f"gleam_nside{nside}_{spw}.hdf5", clobber=True)

    for lsts_center in lsts_center_array:
        if lsts_center == "0h":
            n0 = 0
        elif lsts_center == "1h":
            n0 = 1
        elif lsts_center == "2h":
            n0 = 2
        else:
            n0 = 3

        times_total = get_times_total([uvh5name[4 * i + n0] for i in range(N)])
        print("times_total shape:", times_total.shape)

        template_name = path_data + f"Full_LST_template_{lsts_center}_{spw}.uvh5"
        make_uvdata_template(freqs, times_total, ants, antpos, template_name)

        # If the two downstream runs really need separate files, keep both.
        # Otherwise, you can use one template and reference it twice.
        import shutil
        diff_name = path_data + f"Full_LST_diffuse_model_final_{lsts_center}_{spw}.uvh5"
        gleam_name = path_data + f"Full_LST_gleam_model_final_{lsts_center}_{spw}.uvh5"
        shutil.copyfile(template_name, diff_name)
        shutil.copyfile(template_name, gleam_name)


_BEAM = None

def init_worker(beam_file):
    global _BEAM
    _BEAM = healvis.beam_model.PowerBeam(beam_file)
    _BEAM.interpolation_function = "az_za_simple"

def simulate_data(job):
    sky_file, uvd_file, freqs = job
    global _BEAM

    beam = copy.deepcopy(_BEAM)
    beam.interp_freq(freqs, inplace=True, kind="cubic")

    fchans = np.arange(len(freqs))
    print("now simulating", uvd_file)

    healvis.simulator.run_simulation_partial_freq(
        fchans,
        uvd_file,
        sky_file,
        beam=beam,
        fov=180,
        smooth_beam=False,
    )
    print("all done!")
    return 0

jobs = []
for spw in spw_array:
    sky_file_diffuse = path_data + f"gsm2016_nside{nside}_{spw}.hdf5"
    sky_file_gleam = path_data + f"gleam_nside{nside}_{spw}.hdf5"

    for lsts_center in lsts_center_array:
        uvd_file_gleam = path_data + f"Full_LST_gleam_model_final_{lsts_center}_{spw}.uvh5"
        uvd_file_diffuse = path_data + f"Full_LST_diffuse_model_final_{lsts_center}_{spw}.uvh5"

        uvd1 = UVData()
        uvd1.read(uvd_file_gleam, read_data=False)
        freqs = uvd1.freq_array.ravel()

        jobs.append([sky_file_diffuse, uvd_file_diffuse, freqs])
        jobs.append([sky_file_gleam, uvd_file_gleam, freqs])

nproc = min(mp.cpu_count(), len(jobs))
with mp.Pool(processes=nproc, initializer=init_worker, initargs=(BEAM_FILE,)) as pool:
    list(pool.imap_unordered(simulate_data, jobs, chunksize=1))        