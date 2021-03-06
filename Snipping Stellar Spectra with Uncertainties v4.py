# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 14:11:18 2020

@author: David
"""

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import optimize

snip_width_13CO = 0.0006
snip_width_C18O = 0.00015

#set base directory: Laptop
#base_dir = 'C:\\Users\\there\\Desktop\\ExoStuff_GitHub\\'

#set base directory: Desktop
base_dir = 'C:\\Users\\d338c921\\GitHub\\ExoStuff'

base_dir_spect = base_dir + '\\ST_spectra\\'
base_dir_models = base_dir + '\\Solar_Models\\'


obj_list = []

for root, dirs, files in os.walk(base_dir_spect):
    for file in files:
        if file.endswith('final.fits'):
            obj_list.append(file)
            
beta_shifts = [-1.9263267897412533e-04, -1.1545041431507331e-04, 1.9379144076145112e-04, 6.888621612746786e-05, 2.4225025924415754e-05, -1.5905756391913187e-04]

from hapi import *
db_begin("hapi_data")

from astropy import units as u
nu_min = (4.7e-6 * u.m).to(u.cm)**(-1)
nu_max = (4.6e-6 * u.m).to(u.cm)**(-1)

#12CO lines
fetch('CO', 5, 1, nu_min.value, nu_max.value)
nu_12_CO, sw_12_CO = getColumns('CO', ['nu', 'sw'])
#13CO lines
fetch('CO', 5, 2, nu_min.value, nu_max.value)
nu_13_CO, sw_13_CO = getColumns('CO', ['nu', 'sw'])
#C18O lines
fetch('CO', 5, 3, nu_min.value, nu_max.value)
nu_C_18_O, sw_C_18_O = getColumns('CO', ['nu', 'sw'])

#convert lines from wavenumber to wavelength in microns
wl_12_CO = nu_12_CO**(-1) * 10000
wl_12_CO = np.sort(wl_12_CO)
wl_13_CO = nu_13_CO**(-1) * 10000
wl_13_CO = np.sort(wl_13_CO)
wl_C_18_O = nu_C_18_O**(-1) * 10000
wl_C_18_O = np.sort(wl_C_18_O)

####################################################################################################################################################################
#write a function to load C & O isotopic abundance ratios from Kobayashi+ 2011
def kobayashi_GCE(species):
    if(species == "C"):
        filename = base_dir + "\\gcemodels\\k11_feh-crat_solarneighborhood.csv"
        feh_crat = pd.read_csv(filename, sep = ',')
        feh_crat = (pd.DataFrame.to_numpy(feh_crat)).T
        feh = feh_crat[0]
        C_13 = feh_crat[1]
        
        return feh, C_13
    
    elif(species == "O"):
        filename = base_dir + "\\gcemodels\\k11_feh-orat_solarneighborhood.csv"
        feh_crat = pd.read_csv(filename, sep = ',')
        feh_crat = (pd.DataFrame.to_numpy(feh_crat)).T
        feh = feh_crat[0]
        O_18 = feh_crat[1]
        O_17 = feh_crat[2]
        
        return feh, O_17, O_18
    
####################################################################################################################################################################
#load star data: wl, flux, err ONLY
def star_data(star):
#select a star
#     star:
#         1 = HIP 102040
#         2 = HIP 29432
#         3 = HIP 42333
#         4 = HIP 77052
#         5 = HIP 79672
#         6 = HIP 85042

    if star in [1, 2, 3, 4, 5, 6]:
    #star = 1, 2, 3, 4, 5, 6
        choose_file = int(star - 1)
        obj_path = base_dir_spect + obj_list[choose_file]
        spect = pyfits.getdata(obj_path)
        obj_name = pyfits.getval(obj_path, 'OBJECT')
        
        #create variables for wavelength, flux, and uncertainty; fits files 1, 2, and 3 respectively
        wl = spect[0,:]
        flux = spect[1,:]
        err = spect[2,:]
        
        #normalize err and flux
        wl = wl / (beta_shifts[choose_file] + 1)
        err = err/np.nanmean(flux)
        flux = flux/np.nanmean(flux) #nanmean or nanmedian?
        
    return wl, flux, err

####################################################################################################################################################################
#load solar model data
def solar_data(CO_species, temp, model_num):
    ###13CO models
    if(CO_species == "13CO" and temp == 5690):
        if(model_num == 1):
            sun = pyfits.getdata(base_dir_models + "nlte5690_4.39.0x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 2):
            sun = pyfits.getdata(base_dir_models + "nlte5690_4.39.0.3x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 3):
            sun = pyfits.getdata(base_dir_models + "nlte5690_4.39.1x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 4):
            sun = pyfits.getdata(base_dir_models + "nlte5690_4.39.3x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
        
    elif(CO_species == "13CO" and temp == 5780):
        if(model_num == 1):
            sun = pyfits.getdata(base_dir_models + "nlte5780_4.44.0x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 2):
            sun = pyfits.getdata(base_dir_models + "nlte5780_4.44.0.3x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 3):
            sun = pyfits.getdata(base_dir_models + "nlte5780_4.44.1x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 4):
            sun = pyfits.getdata(base_dir_models + "nlte5780_4.44.3x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?

    elif(CO_species == "13CO" and temp == 5870):
        if(model_num == 1):
            sun = pyfits.getdata(base_dir_models + "nlte5870_4.49.0x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 2):
            sun = pyfits.getdata(base_dir_models + "nlte5870_4.49.0.3x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 3):
            sun = pyfits.getdata(base_dir_models + "nlte5870_4.49.1x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
    
        elif(model_num == 4):
            sun = pyfits.getdata(base_dir_models + "nlte5870_4.49.3x.fromsun.hires.spec.fits")
            wl = sun[0,:] / 10000
            sun_flux_log = sun[1,:] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?

    
    ###C18O models
    elif(CO_species == "C18O" and temp == 5690):
        if(model_num == 1): 
            filename = base_dir_models + "nlte5.69e+03.4.39e+00.0.0.O18=0.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 2): 
            filename = base_dir_models + "nlte5.69e+03.4.39e+00.0.0.O18=3.00e-01x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 3): 
            filename = base_dir_models + "nlte5.69e+03.4.39e+00.0.0.O18=5.60e-01x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 4): 
            filename = base_dir_models + "nlte5.69e+03.4.39e+00.0.0.O18=1.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 5): 
            filename = base_dir_models + "nlte5.69e+03.4.39e+00.0.0.O18=1.78e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 6): 
            filename = base_dir_models + "nlte5.69e+03.4.39e+00.0.0.O18=3.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
    elif(CO_species == "C18O" and temp == 5780):
        if(model_num == 1): 
            filename = base_dir_models + "nlte5.78e+03.4.44e+00.0.0.O18=0.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 2): 
            filename = base_dir_models + "nlte5.78e+03.4.44e+00.0.0.O18=3.00e-01x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 3): 
            filename = base_dir_models + "nlte5.78e+03.4.44e+00.0.0.O18=5.60e-01x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 4): 
            filename = base_dir_models + "nlte5.78e+03.4.44e+00.0.0.O18=1.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 5): 
            filename = base_dir_models + "nlte5.78e+03.4.44e+00.0.0.O18=1.78e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 6): 
            filename = base_dir_models + "nlte5.78e+03.4.44e+00.0.0.O18=3.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
    elif(CO_species == "C18O" and temp == 5870):
        if(model_num == 1): 
            filename = base_dir_models + "nlte5.87e+03.4.49e+00.0.0.O18=0.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 2): 
            filename = base_dir_models + "nlte5.87e+03.4.49e+00.0.0.O18=3.00e-01x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 3): 
            filename = base_dir_models + "nlte5.87e+03.4.49e+00.0.0.O18=5.60e-01x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 4): 
            filename = base_dir_models + "nlte5.87e+03.4.49e+00.0.0.O18=1.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 5): 
            filename = base_dir_models + "nlte5.87e+03.4.49e+00.0.0.O18=1.78e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
            
        elif(model_num == 6): 
            filename = base_dir_models + "nlte5.87e+03.4.49e+00.0.0.O18=3.00e+00x.fromsun.hires.7.csv"
            sun = pd.read_csv(filename, sep = ' ')
            sun = (pd.DataFrame.to_numpy(sun)).T
            wl = sun[0] / 10000
            sun_flux_log = sun[1] #log base 10 of flux
            sun_flux = 10**sun_flux_log
            flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?

    return wl, flux

#####################################################################################################################################################################
#snipping functions
#no errors
def snip_simp(wavegrid, spectrum, linelocation, width):
    index = np.logical_and(wavegrid > (linelocation - width/2), wavegrid < (linelocation + width/2))
    return wavegrid[index], spectrum[index]
#with errors
def snip(wavegrid, spectrum, error, linelocation, width):
    index = np.logical_and(wavegrid > (linelocation - width/2), wavegrid < (linelocation + width/2))
    return wavegrid[index], spectrum[index], error[index]

#####################################################################################################################################################################
def useable_lines(line_list, skip_list):
    new_list = []
    for i in range(len(line_list)):
        if i not in skip_list:
            new_list = np.append(new_list, line_list[i])
    return new_list

#####################################################################################################################################################################
def snip_spectrum(line_list, spectrum, error, wavegrid, snip_width, interp_sample_size): #indices_per_width uniform in solar spectrum but NOT in solar twin spectra
    snips_wl = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_flux = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_flux_norm = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_err = np.ndarray(shape = (len(line_list), interp_sample_size))
    for i in range(len(line_list)):
        wavegrid_mini = np.linspace(-snip_width/2, snip_width/2, interp_sample_size)
        wl_mini , spec_mini, err_mini = snip(wavegrid, spectrum, error, line_list[i], snip_width)
        spec_mini_interp = np.interp(wavegrid_mini, wl_mini - line_list[i], spec_mini)
        err_mini_interp = np.interp(wavegrid_mini, wl_mini - line_list[i], err_mini) #Do I need to normalize this as well?
        spec_norm = spec_mini_interp / np.nanmean(spec_mini_interp) #normalize each snip to flux ~1
        snips_wl[i, :] = np.linspace(line_list[i] - (snip_width/2), line_list[i] + (snip_width/2), interp_sample_size)
        snips_flux[i, :] = spec_mini_interp
        snips_flux_norm[i, :] = spec_norm
        snips_err[i, :] = err_mini_interp / np.nanmean(spec_mini_interp)
    return snips_wl.T, snips_flux.T, snips_flux_norm.T, snips_err.T

#no errors
def snip_spectrum_simp(line_list, spectrum, wavegrid, snip_width, interp_sample_size): #indices_per_width uniform in solar spectrum but NOT in solar twin spectra
    snips_wl = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_flux = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_flux_norm = np.ndarray(shape = (len(line_list), interp_sample_size))
    for i in range(len(line_list)):
        wavegrid_mini = np.linspace(-snip_width/2, snip_width/2, interp_sample_size)
        wl_mini , spec_mini = snip_simp(wavegrid, spectrum, line_list[i], snip_width)
        spec_mini_interp = np.interp(wavegrid_mini, wl_mini - line_list[i], spec_mini)
        spec_norm = spec_mini_interp / np.nanmean(spec_mini_interp) #normalize each snip to flux ~1
        snips_wl[i, :] = np.linspace(line_list[i] - (snip_width/2), line_list[i] + (snip_width/2), interp_sample_size)
        snips_flux[i, :] = spec_mini_interp
        snips_flux_norm[i, :] = spec_norm
    return snips_wl.T, snips_flux.T, snips_flux_norm.T
###WILL RETURN AN ERROR IF LINES FROM DISCONTINOUS REGIONS IN THE SPECTRUM ARE USED

#####################################################################################################################################################################
#Use this function to identify bad lines
def id_bad_lines(star_num, solar_skiplist, star_skiplist, CO_species, line_by_line, save, label = "ST X"):
    
    #load star data
    wl, flux, err = star_data(star_num)
    
    if(CO_species == "13CO"):
        full_line_list = wl_13_CO
    elif(CO_species == "C18O"):
        full_line_list = wl_C_18_O
    
    lines = useable_lines(full_line_list, np.append(solar_skiplist, star_skiplist))
    snips_wl, snips_flux, snips_flux_norm, snips_err = snip_spectrum(lines, flux, err, wl, 0.0006, 50)
    
    if line_by_line == True:
        plt.figure(figsize=(15,15))
        plt.suptitle(label + ": Selected " +  CO_species + " Lines")
        for i in range(len(lines)):
            plt.subplot(int(np.ceil(len(lines) / 2)), 2, i+1) ###adjust size
            plt.title("Line " + str(i + 1))
            plt.xlabel("Wavelength (um)")
            plt.ylabel("Normalized Flux Intensity")
            plt.plot(snips_wl[:, i], snips_flux[:, i])
            plt.errorbar(snips_wl[:, i], snips_flux[:, i], yerr = snips_err[:, i], fmt = 'none', linewidth = 0.5, color = "black")
            plt.axvline(x = lines[i], color = "k", linestyle = ':')
            # if save == True:
            #     fig_path = base_dir + "\\5_28_2021 C18O Snips\\" + label
            #     fig_filename = "\\line" + str(i + 1) + "snip.png"
            #     fig.savefig(fig_path + fig_filename)
        plt.subplots_adjust(top = 0.935, wspace = 0.20, hspace = 0.55)
        #plt.subplots_adjust()
    if(line_by_line == True and save == False):
        print(lines)
            
    # #plot snips individually in original wavlength grid
    # plt.figure()
    # plt.title(label + ": Spectral Snips for Selected Lines: " + CO_species)
    # plt.xlabel("Wavelength (um)")
    # plt.ylabel("Normalized Flux Intensity")
    # plt.plot(snips_wl, snips_flux)
    # for i in range(len(wl_C_18_O)):
    #     plt.axvline(x = wl_C_18_O[i], color = "k", linestyle = ':')
    
    # #plot stacked snips centered around v = 0
    # plt.figure()
    # plt.title(label + ": Stacked Spectral Snips: " + CO_species)
    # plt.xlabel("Velocity (km/s)")
    # plt.ylabel("Normalized Flux Intensity")
    # plt.ylim(0.8, 1.2)
    # plt.plot((snips_wl - lines) / lines * 3e5, snips_flux)
    
    return

#####################################################################################################################################################################
#Use this function to identify bad lines
def id_bad_lines_SM(solar_skiplist, star_skiplist, CO_species, temp, line_by_line, save, label = " ST X"):
    
    if(CO_species == "13CO"):
        full_line_list = wl_13_CO
    elif(CO_species == "C18O"):
        full_line_list = wl_C_18_O
    
    lines = useable_lines(full_line_list, np.append(solar_skiplist, star_skiplist))
    
    if line_by_line == True:
        plt.figure(figsize=(15,45))
        plt.suptitle("T = " + str(temp) + label + ": Selected " +  CO_species + " Lines")
        plt.xlabel("Wavelength (um)")
        plt.ylabel("Normalized Flux Intensity")
        for i in range(len(lines)):
            if (len(lines) <= 12):
                plt.subplot(int(np.ceil(len(lines) / 2)), 2, i + 1) ###adjust size
                plt.subplots_adjust(top = 0.935, wspace = 0.20, hspace = 0.55)
            elif (len(lines) >= 12):
                plt.subplot(int(np.ceil(len(lines) / 4)), 4, i + 1) ###adjust size
                plt.subplots_adjust(top = 0.935, wspace = 0.20, hspace = 0.55)
            plt.title("Line " + str(i + 1))
            plt.xlim(lines[i] - 0.0003, lines[i] + 0.0003)
            plt.axvline(x = lines[i], color = "k", linestyle = ':')
            for j in range(5):
                wl, flux = solar_data(CO_species, temp, j + 1)
                wl, flux = snip_simp(wl, flux, lines[i], 0.0006)
                plt.ylim(np.mean(flux) - 0.005, np.mean(flux) + 0.005)
                plt.plot(wl, flux)
            # if save == True:
            #     fig_path = base_dir + "\\5_28_2021 C18O Snips\\" + label
            #     fig_filename = "\\line" + str(i + 1) + "snip.png"
            #     fig.savefig(fig_path + fig_filename)
        #plt.subplot_tool()
    if(line_by_line == True and save == False):
        print(lines)
    
    return

####################################################################################################################################################################
#Write a function to output stack_vel and stack_flux for a given star/model i.e. same result as snip_plots with output_stack = True but without the plots
def stack_data(star, useable_line_list, snip_width, CO_species, temp, model):
            
    if (model == True):
        #Solar Models
        wl, flux = solar_data(CO_species, temp, star)
        snips_wl, snips_flux, snips_flux_norm = snip_spectrum_simp(useable_line_list, flux, wl, snip_width, 50)
        stack_vel = np.mean((snips_wl - useable_line_list) / useable_line_list * 3e5, axis = 1)
        stack_flux = np.mean(snips_flux_norm, axis = 1)
        stack_flux = stack_flux - (np.max(stack_flux) - 1) #remember to lower models so they max at ~1
        
        return stack_vel, stack_flux
    
    elif (model == False):
        wl, flux, err = star_data(star)
        snips_wl, snips_flux, snips_flux_norm, snips_err = snip_spectrum(useable_line_list, flux, err, wl, snip_width, 50)
        #stars w/ velocity shift corrections
        if CO_species == "13CO":
            stack_vel = np.mean((snips_wl - useable_line_list) / useable_line_list * 3e5, axis = 1) + velocity_shifts[star - 1]
        elif CO_species == "C18O":
            stack_vel = np.mean((snips_wl - useable_line_list) / useable_line_list * 3e5, axis = 1) + velocity_shifts_C18O[star - 1]
        stack_flux = np.average(snips_flux_norm, axis = 1, weights = snips_err**-2) ###weighted average using weights 1/variance
        stack_err = np.sqrt(1./ (snips_err**-2).sum(axis = 1))
        
        return stack_vel, stack_flux, stack_err

#####################################################################################################################################################################
#Write a fuction to find this line for a given stacked absorption line; outputs np.poly1d(fit_coefficients) to be divided from the spectrum later
def de_slant(stack_vel, stack_flux, L_shoulder, R_shoulder, shoulder_width):
    #returns a np.poly1d polynomial (a line) to be divided out from the spectrum in order to correct the slant
    
    vel_snipL, flux_snipL = snip_simp(stack_vel, stack_flux, L_shoulder, shoulder_width) #snips the left "shoulder" of the stacked absorption line i.e. [-10, -5] km/s
    vel_snipR, flux_snipR = snip_simp(stack_vel, stack_flux, R_shoulder, shoulder_width) #snips the right "shoulder" of the stacked absorption line i.e. [-10, -5] km/s
    
    avgL = np.nanmean(flux_snipL) #average flux of the left shoulder
    avgR = np.nanmean(flux_snipR) #average flux of the right shoulder
    
    l_coeff = np.polyfit([L_shoulder, R_shoulder], [avgL, avgR], deg = 1) #fits a line through the midpoints of each shoulder i.e. (- 7.5 km/s, avgL) and ( 7.5 km/s, avgR); returns coeffecients 
    line = np.poly1d(l_coeff) #puts the equation of the line into a "useable form"; divide out of spectrum as line(stack_vel)
    
    return line

#####################################################################################################################################################################
#Write a function to calculate chi-squared values for the line profiles plotted in abundance_plot()
def chi_sqr_lp(star_num, stack_vel_obs, stack_vel_model, stack_flux_obs, stack_flux_model, err, CO_species):
    #custom chi_sqr snip widths and line centers; lp = line profiles
    if CO_species == "13CO":
        snip_bounds_L = [-5.864, -3.405, -4.429, -4.016, -5.227, -6.35]
        snip_bounds_R = [6.118, 3.627, 6.235, 4.444, 4.273, 4.83]
    elif CO_species == "C18O":
        snip_bounds_L = [-4.84, -5.2, -3.59, -3.55, -3.7, -4.64]
        snip_bounds_R = [4.4, 4.8, 3.91, 3.2, 3.67, 4.64]
        
    chi_sqr_snipwidth = snip_bounds_R[star_num - 1] - snip_bounds_L[star_num - 1]
    snip_center = snip_bounds_R[star_num - 1] - (chi_sqr_snipwidth / 2)
    
    #resize arrays to exclude the wings
    stack_vel_obs, stack_flux_obs, err = snip(stack_vel_obs, stack_flux_obs, err, snip_center, chi_sqr_snipwidth)
    stack_vel_model, stack_flux_model = snip_simp(stack_vel_model, stack_flux_model, snip_center, chi_sqr_snipwidth)
    
    chisqr = 0
    for i in range(len(stack_flux_obs) - 1):
        chisqr = chisqr + ((stack_flux_obs[i] - stack_flux_model[i])**2)/(err[i]**2)
    return chisqr

#####################################################################################################################################################################
#Write a function to calculate chi-squared values for regions of spectrum
def chi_sqr(flux_obs, flux_model, err):
    
    chisqr = 0
    for i in range(len(flux_obs) - 1):
        chisqr = chisqr + ((flux_obs[i] - flux_model[i])**2)/(err[i]**2)
    return chisqr

#####################################################################################################################################################################
def test_stellar_parameters(star_num, temp_guess, plot, label = "HIP ###"):
    #this function resembles the calc_abundance() function; however, here we will use the original stellar spectrum and compare it to 1x Solar models of 
    #different temperatures and log g values to see which set of models we should use for line profile comparison
    
    #load data
    wl, flux, err = star_data(star_num)
    sun_wl1, sun_flux1 = solar_data(CO_species = "13CO", temp = 5690, model_num = 3) #1x Solar with T = 5690, log g = 4.39
    sun_wl2, sun_flux2 = solar_data(CO_species = "13CO", temp = 5780, model_num = 3) #1x Solar with T = 5780, log g = 4.44
    sun_wl3, sun_flux3 = solar_data(CO_species = "13CO", temp = 5870, model_num = 3) #1x Solar with T = 5870, log g = 4.49
    
    #Use snip functions to select a wavelength region on which to perform a chi-squared test

    sun_wl1, sun_flux1 = snip_simp(sun_wl1, sun_flux1, 4.618, width = 0.016) #snip_center = center of the region of interest; width = width of that region
    sun_wl2, sun_flux2 = snip_simp(sun_wl2, sun_flux2, 4.618, width = 0.016)
    sun_wl3, sun_flux3 = snip_simp(sun_wl3, sun_flux3, 4.618, width = 0.016)
    
    wl, flux, err = snip(wl, flux, err, 4.618, width = 0.016) #snip_center = center of the region of interest; width = width of that region
    
    #interpolate models onto ST wavegrid
    sun_flux1 = np.interp(wl, sun_wl1, sun_flux1)
    sun_flux2 = np.interp(wl, sun_wl2, sun_flux2)
    sun_flux3 = np.interp(wl, sun_wl3, sun_flux3)
    
    #normalize models ~1
    sun_flux1 = sun_flux1 / np.mean(sun_flux1)
    sun_flux2 = sun_flux2 / np.mean(sun_flux2)
    sun_flux3 = sun_flux3 / np.mean(sun_flux3)
    
    # #deslant sun_fluxes???
    # slant1 = de_slant(wl, sun_flux1, 4.603, 4.633, 0.006)
    # slant2 = de_slant(wl, sun_flux2, 4.603, 4.633, 0.006)
    # slant3 = de_slant(wl, sun_flux3, 4.603, 4.633, 0.006)
    
    # sun_flux1 = sun_flux1 / slant1
    # sun_flux2 = sun_flux2 / slant2
    # sun_flux3 = sun_flux3 / slant3
    
    #calculate 3 chi_sqr values and put them in an array
    chi_sqr_list = [chi_sqr(flux, sun_flux1, err), chi_sqr(flux, sun_flux2, err), chi_sqr(flux, sun_flux3, err)]
        
    #print(chi_sqr_list)
    #create a matching "x-array" of xSolar abundances i.e. [0, 1/3, 1, 3]
    x_list = [5690, 5780, 5870]
    
    #fit them with a degree 2 polynomial
    p_coeff = np.polyfit(x_list, chi_sqr_list, deg = 2)
    parabola = np.poly1d(p_coeff)
    
    #find the minimum of the parabola
    fit_min, fit_min_val, num_its, f_evals, warn_flag = optimize.fmin(parabola, temp_guess, full_output = True)
    #find where polynomial fit = fit_min_val + 1
    err_bounds = np.roots(parabola - fit_min_val - 1)
    temp_err = np.abs(err_bounds[1] - err_bounds[0]) / 2
    
    #plot if option plot= True
    if (plot == True):
        plt.figure()
        plt.title(label + ": Temperature Fit")
        plt.xlabel("Temperature")
        plt.ylabel("$\chi^2$")
        plt.scatter(x_list, chi_sqr_list, marker = 'o', label = "$\chi^2$ Interp")
        plt.scatter(err_bounds, parabola(err_bounds), color = "r", marker = "x")
        plt.plot(np.linspace(5000, 6000, 1000), parabola(np.linspace(5000, 6000, 1000)), label = "Data Fit")
        plt.plot(fit_min, parabola(fit_min), marker = "X", label = ("min @ x = %s" % np.round(fit_min, 0))) #different color than previous
        plt.legend()
        
    else:
        pass
    
    return np.round(fit_min, 0), np.round(temp_err, 0)

#####################################################################################################################################################################
def calc_abundance(star_num, full_line_list, solar_skiplist, star_skiplist, snip_width, abundance_guess, plot, CO_species, temp, label = "HIP ###"):
    ###requires stack_data() and chi_sqr_lp() functions
    #requires scipy.optimize
    
    #star_num = a number 1-6; used to pick which star we're looking at
    #line_list = list of lines to skip in stacked absorption line analysis
    #how wide to make the snips; default = 0.0006
    #abundance_guess = number provided for finding the chi-sqr minimum i.e. your best guess at what the xSolar abundance is
    #if plot = true, show the plot of chi-sqr values, data fit, and min
    
    #stack_data for chosen star and all four solar models
    line_list = useable_lines(full_line_list, np.append(star_skiplist, solar_skiplist))
    stack_vel, stack_flux, stack_err = stack_data(star_num, line_list, snip_width, CO_species, temp, model = False)
    #SLANT CORRECTION!!!
    
    if (star_num == 6 and CO_species == "13CO"):
        slant = de_slant(stack_vel, stack_flux, -6.5, 4.8, 2)
    else:
        slant = de_slant(stack_vel, stack_flux, -7.5, 7.5, 5)
    
    stack_flux = stack_flux / slant(stack_vel)
    stack_err = stack_err / slant(stack_vel)
    
    if CO_species == "13CO":
        stack_vel1, stack_flux1 = stack_data(1, line_list, snip_width, CO_species, temp, model = True)
        stack_vel2, stack_flux2 = stack_data(2, line_list, snip_width, CO_species, temp, model = True)
        stack_vel3, stack_flux3 = stack_data(3, line_list, snip_width, CO_species, temp, model = True)
        stack_vel4, stack_flux4 = stack_data(4, line_list, snip_width, CO_species, temp, model = True)
        plt_title = " Abundance Calculation: $\chi^2$ Test Interp: 13CO"
            
        #calculate 4 chi_sqr values and put them in an array
        chi_sqr_list = [chi_sqr_lp(star_num, stack_vel, stack_vel1, stack_flux, stack_flux1, stack_err, CO_species = "13CO"), chi_sqr_lp(star_num, stack_vel, stack_vel2, stack_flux, stack_flux2, stack_err, CO_species = "13CO"), chi_sqr_lp(star_num, stack_vel, stack_vel3, stack_flux, stack_flux3, stack_err, CO_species = "13CO"), chi_sqr_lp(star_num, stack_vel, stack_vel4, stack_flux, stack_flux4, stack_err, CO_species = "13CO")]
        
        #print(chi_sqr_list)
        #create a matching "x-array" of xSolar abundances i.e. [0, 1/3, 1, 3]
        x_list = [0.0, 1/3, 1.0, 3.0]
        
    elif CO_species == "C18O":
        stack_vel1, stack_flux1 = stack_data(1, line_list, snip_width, CO_species, temp, model = True)
        stack_vel2, stack_flux2 = stack_data(2, line_list, snip_width, CO_species, temp, model = True)
        stack_vel3, stack_flux3 = stack_data(3, line_list, snip_width, CO_species, temp, model = True)
        stack_vel4, stack_flux4 = stack_data(4, line_list, snip_width, CO_species, temp, model = True)
        stack_vel5, stack_flux5 = stack_data(5, line_list, snip_width, CO_species, temp, model = True)
        stack_vel6, stack_flux6 = stack_data(6, line_list, snip_width, CO_species, temp, model = True)
        plt_title = " Abundance Calculation: $\chi^2$ Test Interp: C18O"
            
        #calculate 6 chi_sqr values and put them in an array
        chi_sqr_list = [chi_sqr_lp(star_num, stack_vel, stack_vel1, stack_flux, stack_flux1, stack_err, CO_species = "C18O"), chi_sqr_lp(star_num, stack_vel, stack_vel2, stack_flux, stack_flux2, stack_err, CO_species = "C18O"), chi_sqr_lp(star_num, stack_vel, stack_vel3, stack_flux, stack_flux3, stack_err, CO_species = "C18O"), chi_sqr_lp(star_num, stack_vel, stack_vel4, stack_flux, stack_flux4, stack_err, CO_species = "C18O"), chi_sqr_lp(star_num, stack_vel, stack_vel5, stack_flux, stack_flux5, stack_err, CO_species = "C18O"), chi_sqr_lp(star_num, stack_vel, stack_vel6, stack_flux, stack_flux6, stack_err, CO_species = "C18O")]
        
        #print(chi_sqr_list)
        #create a matching "x-array" of xSolar abundances i.e. [0, 1/3, 1, 3]
        x_list = [0.0, 0.3, 0.56, 1.0, 1.78, 3.0]
    
    #fit them with a degree 2 polynomial
    p_coeff = np.polyfit(x_list, chi_sqr_list, deg = 2)
    parabola = np.poly1d(p_coeff)
    
    #find the minimum of the parabola
    fit_min, fit_min_val, num_its, f_evals, warn_flag = optimize.fmin(parabola, abundance_guess, full_output = True)
    #find where polynomial fit = fit_min_val + 1
    err_bounds = np.roots(parabola - fit_min_val - 1)
    abundance_err = np.abs(err_bounds[1] - err_bounds[0]) / 2
    fit_min_dex = np.log10(fit_min)
    abundance_err_dex = abundance_err / fit_min #err through a log: log(x+dx) = log(x) + dx/x
    
    #plot if option plot= True
    if (plot == True):
        plt.figure()
        plt.title(label + plt_title)
        plt.xlabel("xSolar Abundance")
        plt.ylabel("$\chi^2$")
        plt.scatter(x_list, chi_sqr_list, marker = 'o', label = "$\chi^2$ Interp")
        plt.scatter(err_bounds, parabola(err_bounds), color = "r", marker = "x")
        plt.plot(np.linspace(0, 3, 1000), parabola(np.linspace(0, 3, 1000)), label = "Data Fit")
        plt.plot(fit_min, parabola(fit_min), marker = "X", label = ("min @ x = %s" % np.round(fit_min, 2))) #different color than previous
        plt.legend()
        
    else:
        pass
    
    return np.round(fit_min, 2), np.round(fit_min_dex, 2), np.round(abundance_err, 2), np.round(abundance_err_dex[0], 2)

#####################################################################################################################################################################
#Write a function to automate the abundance plots i.e. stellar stacked absorption lines + 4 model stacked lines
def abundance_plot(star_num, full_line_list, solar_skiplist, star_skiplist, snip_width, CO_species, temp, label = "HIP ###"):
    #make suree to use same line lists for solar models and stars
    line_list = useable_lines(full_line_list, np.append(star_skiplist, solar_skiplist))
    
    stack_vel, stack_flux, stack_err = stack_data(star_num, line_list, snip_width, CO_species, temp, model = False)
    
    if (star_num == 6 and CO_species == "13CO"):
        slant = de_slant(stack_vel, stack_flux, -6.5, 4.8, 2)
    else:
        slant = de_slant(stack_vel, stack_flux, -7.5, 7.5, 5)
    
    if CO_species == "13CO":
        
        stack_vel_mod1, stack_flux_mod1 = stack_data(1, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod2, stack_flux_mod2 = stack_data(2, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod3, stack_flux_mod3 = stack_data(3, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod4, stack_flux_mod4 = stack_data(4, line_list, snip_width, CO_species, temp, model = True)
        
        label1 = "0x Solar $^{13}$CO"
        label2 = "1/3x Solar $^{13}$CO"
        label3 = "Solar $^{13}$CO"
        label4 = "3x Solar $^{13}$CO"
        
        plt.figure()
        plt.title("Abundance Plot: " + label + ": 13CO")
        plt.ylim(0.91, 1.02)
        plt.xlim(-15, 15)
        plt.plot(stack_vel, stack_flux / slant(stack_vel), label = label, linewidth = 3)
        plt.plot(stack_vel_mod1, stack_flux_mod1, linestyle='--', label = label1)
        plt.plot(stack_vel_mod2, stack_flux_mod2, linestyle='--', label = label2)
        plt.plot(stack_vel_mod3, stack_flux_mod3, linestyle='--', label = label3)
        plt.plot(stack_vel_mod4, stack_flux_mod4, linestyle='--', label = label4)
        plt.errorbar(stack_vel, stack_flux / slant(stack_vel), yerr = stack_err / slant(stack_vel), fmt = 'none', linewidth = 0.5, color = "black")
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("Normalized Flux Intensity")
        plt.axvline(x = 0, color= 'k', linestyle='--')
        plt.legend()
    
    elif CO_species == "C18O":
        
        stack_vel_mod1, stack_flux_mod1 = stack_data(1, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod2, stack_flux_mod2 = stack_data(2, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod3, stack_flux_mod3 = stack_data(3, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod4, stack_flux_mod4 = stack_data(4, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod5, stack_flux_mod5 = stack_data(5, line_list, snip_width, CO_species, temp, model = True)
        stack_vel_mod6, stack_flux_mod6 = stack_data(6, line_list, snip_width, CO_species, temp, model = True)
        
        label1 = "0x Solar C$^{18}$O"
        label2 = "0.3x Solar C$^{18}$O"
        label3 = "0.56x Solar C$^{18}$O"
        label4 = "Solar C$^{18}$O"
        label5 = "1.78x Solar C$^{18}$O"
        label6 = "3x Solar C$^{18}$O"
        
        plt.figure()
        plt.title("Abundance Plot: " + label + ": C18O")
        plt.ylim(0.94, 1.01)
        plt.xlim(-10, 10)
        plt.plot(stack_vel, stack_flux / slant(stack_vel), label = label, linewidth = 3)
        plt.plot(stack_vel_mod1, stack_flux_mod1, linestyle='--', label = label1)
        plt.plot(stack_vel_mod2, stack_flux_mod2, linestyle='--', label = label2)
        plt.plot(stack_vel_mod3, stack_flux_mod3, linestyle='--', label = label3)
        plt.plot(stack_vel_mod4, stack_flux_mod4, linestyle='--', label = label4)
        plt.plot(stack_vel_mod5, stack_flux_mod5, linestyle='--', label = label5)
        plt.plot(stack_vel_mod6, stack_flux_mod6, linestyle='--', label = label6)
        plt.errorbar(stack_vel, stack_flux / slant(stack_vel), yerr = stack_err / slant(stack_vel), fmt = 'none', linewidth = 0.5, color = "black")
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("Normalized Flux Intensity")
        plt.axvline(x = 0, color= 'k', linestyle='--')
        plt.legend()
    
    return

#####################################################################################################################################################################
velocity_shifts = [1.97337 - 1.1885, 0.944558, 0.600705 - 0.425, 0.94817 - 0.455, 0.772493 - 0.88, -0.845]###0.845] #add to velocity grid of stacked absorption lines; + moves spectrum right, - moves left

#plot for sun using solar skip_list ###no flux 13CO lines = 11, 12, 24, 25
solar_skiplist = [0, 1, 2, 6, 12, 14, 15, 16, 19, 21, 23, 25, 26, 27, 29, 32, 33]
sun_13_CO_lines = useable_lines(wl_13_CO, solar_skiplist)


#star 1: HIP 102040
star1_skiplist = [8, 11, 24]
abundance_plot(1, wl_13_CO, solar_skiplist, star1_skiplist, snip_width_13CO, CO_species = "13CO", temp = 5870, label = "Solar Twin 1")
fit_min1, fit_min_dex1, ab_err1, ab_err_dex1 = calc_abundance(1, wl_13_CO, solar_skiplist, star1_skiplist, snip_width_13CO, 1.0, plot=True, CO_species = "13CO", temp = 5870, label = "Solar Twin 1")
print("HIP 102040 has a 13CO abundance of " + str(fit_min1) + " +/- " + str(ab_err1) + " xSolar or " + str(fit_min_dex1) + " dex \n")


#star 2: HIP 29432
star2_skiplist = [8, 10, 11, 18, 24]
abundance_plot(2, wl_13_CO, solar_skiplist, star2_skiplist, snip_width_13CO, CO_species = "13CO", temp = 5780, label = "Solar Twin 2")
fit_min2, fit_min_dex2, ab_err2, ab_err_dex2 = calc_abundance(2, wl_13_CO, solar_skiplist, star2_skiplist, snip_width_13CO, 1.55, plot=True, CO_species = "13CO", temp = 5780, label = "Solar Twin 2")
print("HIP 29432 has a 13CO abundance of " + str(fit_min2) + " +/- " + str(ab_err2) + " xSolar or " + str(fit_min_dex2) + " dex \n")


#star 3: HIP 42333
star3_skiplist = [11, 24]
abundance_plot(3, wl_13_CO, solar_skiplist, star3_skiplist, snip_width_13CO, CO_species = "13CO", temp = 5870, label = "Solar Twin 3")
fit_min3, fit_min_dex3, ab_err3, ab_err_dex3 = calc_abundance(3, wl_13_CO, solar_skiplist, star3_skiplist, snip_width_13CO, 1.5, plot=True, CO_species = "13CO", temp = 5870, label = "Solar Twin 3")
print("HIP 42333 has a 13CO abundance of " + str(fit_min3) + " +/- " + str(ab_err3) + " xSolar or " + str(fit_min_dex3) + " dex \n")


#star 4: HIP 77052
star4_skiplist = [11, 24]
star4_useable_sun_lines = useable_lines(sun_13_CO_lines, star4_skiplist)
abundance_plot(4, wl_13_CO, solar_skiplist, star4_skiplist, snip_width_13CO, CO_species = "13CO", temp = 5690, label = "Solar Twin 4")
fit_min4, fit_min_dex4, ab_err4, ab_err_dex4 = calc_abundance(4, wl_13_CO, solar_skiplist, star4_skiplist, snip_width_13CO, 1.5, plot=True, CO_species = "13CO", temp = 5690, label = "Solar Twin 4")
print("HIP 77052 has a 13CO abundance of " + str(fit_min4) + " +/- " + str(ab_err4) + " xSolar or " + str(fit_min_dex4) + " dex \n")

#star 5: HIP 79672
star5_skiplist = [11, 24]
star5_useable_sun_lines = useable_lines(sun_13_CO_lines, star5_skiplist)
abundance_plot(5, wl_13_CO, solar_skiplist, star5_skiplist, snip_width_13CO, CO_species = "13CO", temp = 5780, label = "Solar Twin 5")
fit_min5, fit_min_dex5, ab_err5, ab_err_dex5 = calc_abundance(5, wl_13_CO, solar_skiplist, star5_skiplist, snip_width_13CO, 0.8, plot=True, CO_species = "13CO", temp = 5780, label = "Solar Twin 5")
print("HIP 79672 has a 13CO abundance of " + str(fit_min5) + " +/- " + str(ab_err5) + " xSolar or " + str(fit_min_dex5) + " dex \n")


#star 6: HIP 85042
star6_skiplist = [8, 11, 24]
star6_useable_sun_lines = useable_lines(sun_13_CO_lines, star6_skiplist)
abundance_plot(6, wl_13_CO, solar_skiplist, star6_skiplist, snip_width_13CO, CO_species = "13CO", temp = 5690, label = "Solar Twin 6")
fit_min6, fit_min_dex6, ab_err6, ab_err_dex6 = calc_abundance(6, wl_13_CO, solar_skiplist, star6_skiplist, snip_width_13CO, 1.7, plot=True, CO_species = "13CO", temp = 5690, label = "Solar Twin 6")
print("HIP 85042 has a 13CO abundance of " + str(fit_min6) + " +/- " + str(ab_err6) + " xSolar or " + str(fit_min_dex6) + " dex \n")

###
print("xSolar abundances: ")
print(fit_min1, fit_min2, fit_min3, fit_min4, fit_min5, fit_min6)
print("Uncertainties: ")
print(ab_err1, ab_err2, ab_err3, ab_err4, ab_err5, ab_err6)

print("Abundances (dex): ")
print(fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6)
print("Uncertainties: ")
print(ab_err_dex1, ab_err_dex2, ab_err_dex3, ab_err_dex4, ab_err_dex5, ab_err_dex6)
print("\n")


###################################################################################################################################################################
# #Identify unuseable C18O lines in the solar model; add them to the solar_skiplist for C18O
# sun_wl, sun_flux = solar_data(8)
# sun_flux = sun_flux + 0.8
# wl, flux, err = star_data(1)
# plt.figure()
# plt.xlim(4.6, 4.7)
# plt.ylim(0.9, 1.1)
# plt.plot(sun_wl, sun_flux)
# for i in range(len(wl_C_18_O)):
#     plt.axvline(x = wl_C_18_O[i], color = "k", linestyle = ':')
# for i in range(len(wl_13_CO)):
#     plt.axvline(x = wl_13_CO[i], color = "k", linestyle = '--')
# for i in range(len(wl_C_18_O)):
#     plt.axvline(x = wl_12_CO[i], color = "k", linestyle = '-')
# plt.figure()
# plt.xlim(4.6, 4.7)
# plt.ylim(0.8, 1.2)
# plt.plot(wl, flux)
# for i in range(len(wl_C_18_O)):
#     plt.axvline(x = wl_C_18_O[i], color = "k", linestyle = '--')

solar_skiplist_C18O = [7, 9, 11, 12, 13, 16, 17, 18, 19, 21, 22, 25, 26, 27, 30, 31] #[12, 13, 26, 27] = no flux
# #NOTE: lines 15 or 24 are the strongest

#####################################################################################################################################################################
snip_width_C18O = 0.0006 #make sure snip width is wide enough for deslant to do its thing!!!
velocity_shifts_C18O = [0, -1.19, 0.405, 0, 0.79, 0.585]

#star 1: HIP 102040
star1_C18O_skiplist = [0, 2, 3, 6, 10, 15, 24] #[0, 1, 3, 5, 10, 14, 15, 23, 28, 29] #[6] #[0, 2, 5, 6, 15, 24] #strongest lines i.e. 15 and 24 are messy; so is line 6
abundance_plot(1, wl_C_18_O, solar_skiplist_C18O, star1_C18O_skiplist, snip_width_C18O, CO_species = "C18O", temp = 5870, label = "Solar Twin 1")
fit_min1_C18O, fit_min_dex1_C18O, ab_err1_C18O, ab_err_dex1_C18O = calc_abundance(1, wl_C_18_O, solar_skiplist_C18O, star1_C18O_skiplist, snip_width_C18O, 0.0, plot=True, CO_species = "C18O", temp = 5870, label = "Solar Twin 1")
print("HIP 102040 has a C18O abundance of " + str(fit_min1_C18O) + " +/- " + str(ab_err1_C18O) + " xSolar or " + str(fit_min_dex1_C18O) + " dex \n")


#star 2: HIP 29432
star2_C18O_skiplist = [0, 1, 2, 3, 15, 24] #[0, 1, 2, 3, 5, 10, 14, 15, 23, 24, 28, 29]
abundance_plot(2, wl_C_18_O, solar_skiplist_C18O, star2_C18O_skiplist, snip_width_C18O, CO_species = "C18O", temp = 5780, label = "Solar Twin 2")
fit_min2_C18O, fit_min_dex2_C18O, ab_err2_C18O, ab_err_dex2_C18O = calc_abundance(2, wl_C_18_O, solar_skiplist_C18O, star2_C18O_skiplist, snip_width_C18O, 3.0, plot=True, CO_species = "C18O", temp = 5780, label = "Solar Twin 2")
print("HIP 29432 has a C18O abundance of " + str(fit_min2_C18O) + " +/- " + str(ab_err2_C18O) + " xSolar or " + str(fit_min_dex2_C18O) + " dex \n")


#star 3: HIP 42333
star3_C18O_skiplist = [0, 2, 3, 15, 24] #[24] #24 optional; profile looked fine before it was removed
abundance_plot(3, wl_C_18_O, solar_skiplist_C18O, star3_C18O_skiplist, snip_width_C18O, CO_species = "C18O", temp = 5870, label = "Solar Twin 3")
fit_min3_C18O, fit_min_dex3_C18O, ab_err3_C18O, ab_err_dex3_C18O = calc_abundance(3, wl_C_18_O, solar_skiplist_C18O, star3_C18O_skiplist, snip_width_C18O, 0.5, plot=True, CO_species = "C18O", temp = 5870, label = "Solar Twin 3")
print("HIP 42333 has a C18O abundance of " + str(fit_min3_C18O) + " +/- " + str(ab_err3_C18O) + " xSolar or " + str(fit_min_dex3_C18O) + " dex \n")


#star 4: HIP 77052
star4_C18O_skiplist = [0, 2, 3, 15, 24] #[14] #14 optional; profile looked fine before it was removed
abundance_plot(4, wl_C_18_O, solar_skiplist_C18O, star4_C18O_skiplist, snip_width_C18O, CO_species = "C18O", temp = 5690, label = "Solar Twin 4")
fit_min4_C18O, fit_min_dex4_C18O, ab_err4_C18O, ab_err_dex4_C18O = calc_abundance(4, wl_C_18_O, solar_skiplist_C18O, star4_C18O_skiplist, snip_width_C18O, 1.5, plot=True, CO_species = "C18O", temp = 5690, label = "Solar Twin 4")
print("HIP 77052 has a C18O abundance of " + str(fit_min4_C18O) + " +/- " + str(ab_err4_C18O) + " xSolar or " + str(fit_min_dex4_C18O) + " dex \n")


#star 5: HIP 79672
star5_C18O_skiplist = [0, 2, 3, 15, 24] #stellar & model profiles look GREAT with these lines removed
abundance_plot(5, wl_C_18_O, solar_skiplist_C18O, star5_C18O_skiplist, snip_width_C18O, CO_species = "C18O", temp = 5780, label = "Solar Twin 5")
fit_min5_C18O, fit_min_dex5_C18O, ab_err5_C18O, ab_err_dex5_C18O = calc_abundance(5, wl_C_18_O, solar_skiplist_C18O, star5_C18O_skiplist, snip_width_C18O, 0.75, plot=True, CO_species = "C18O", temp = 5780, label = "Solar Twin 5")
print("HIP 79672 has a C18O abundance of " + str(fit_min5_C18O) + " +/- " + str(ab_err5_C18O) + " xSolar or " + str(fit_min_dex5_C18O) + " dex \n")


#star 6: HIP 85042
star6_C18O_skiplist = [0, 2, 3, 15, 24] #[1, 6]
abundance_plot(6, wl_C_18_O, solar_skiplist_C18O, star6_C18O_skiplist, snip_width_C18O, CO_species = "C18O", temp = 5690, label = "Solar Twin 6")
fit_min6_C18O, fit_min_dex6_C18O, ab_err6_C18O, ab_err_dex6_C18O = calc_abundance(6, wl_C_18_O, solar_skiplist_C18O, star6_C18O_skiplist, snip_width_C18O, 1.78, plot=True, CO_species = "C18O", temp = 5690, label = "Solar Twin 6")
print("HIP 85042 has a C18O abundance of " + str(fit_min6_C18O) + " +/- " + str(ab_err6_C18O) + " xSolar or " + str(fit_min_dex6_C18O) + " dex \n")

###
print("xSolar abundances: ")
print(fit_min1_C18O, fit_min2_C18O, fit_min3_C18O, fit_min4_C18O, fit_min5_C18O, fit_min6_C18O)
print("Uncertainties: ")
print(ab_err1_C18O, ab_err2_C18O, ab_err3_C18O, ab_err4_C18O, ab_err5_C18O, ab_err6_C18O)

print("Abundances (dex): ")
print(fit_min_dex1_C18O, fit_min_dex2_C18O, fit_min_dex3_C18O, fit_min_dex4_C18O, fit_min_dex5_C18O, fit_min_dex6_C18O)
print("Uncertainties: ")
print(ab_err_dex1_C18O, ab_err_dex2_C18O, ab_err_dex3_C18O, ab_err_dex4_C18O, ab_err_dex5_C18O, ab_err_dex6_C18O)


####################################################################################################################################################################
plt.figure()
plt.title("$^{13}$CO Abundance vs. Stellar Age")
plt.xlabel("Stellar Age (Gyr)")
plt.ylabel("Calculated $^{13}$CO Abundance (dex)")
plt.scatter(4.603, 0.0,  color = '#FF6700', marker = "o", label = "Sun")
plt.scatter([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], marker = "o", label = "Solar Twins")
plt.errorbar([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], yerr = [ab_err_dex1, ab_err_dex2, ab_err_dex3, ab_err_dex4, ab_err_dex5, ab_err_dex6], fmt = 'none', linewidth = 0.5, color = "black")
plt.errorbar([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], xerr = [0.91, 0.71, 0.52, 0.91, 0.39, 0.62], fmt = 'none', linewidth = 0.5, color = "black")
age_abundance_fit = np.polyfit([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], deg = 1)
#age_abundance_line = np.poly1d(age_abundance_fit)
#plt.plot(np.linspace(0, 7, 1000), age_abundance_fit[0]*(np.linspace(0, 7, 1000)) + age_abundance_fit[1])
plt.legend()

metallicity = [-0.093, -0.096, 0.138, 0.036, 0.056, 0.015] #Fe/H from spreadsheet
err_m = [0.006, 0.005, 0.008, 0.006, 0.003, 0.004]

plt.figure()
plt.title("$^{13}$CO Abundance vs. Metallicity")
plt.xlabel("Metallicity = Fe/H (dex)")
plt.ylabel("Calculated $^{13}$CO Abundance (dex)")
plt.scatter(0.0, 0.0,  color = '#FF6700', marker = "o", label = "Sun")
plt.scatter(metallicity, [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], marker = "o", label = "Solar Twins")
plt.errorbar(metallicity, [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], yerr = [ab_err_dex1, ab_err_dex2, ab_err_dex3, ab_err_dex4, ab_err_dex5, ab_err_dex6], fmt = 'none', linewidth = 0.5, color = "black")
plt.errorbar(metallicity, [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], xerr = err_m, fmt = 'none', linewidth = 0.5, color = "black")
#age_abundance_fit = np.polyfit(metallicity, [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], deg = 1)
plt.legend()

Teff = [5838, 5758, 5848, 5683, 5814, 5694]
err_Teff = [6, 5, 8, 5, 3, 5]

radius = [0.96, 0.95, 1.0, 0.96, 1.03, 1.04] #Rstar_DR2

plt.figure()
plt.title("T$_{eff}$ vs. Radius")
plt.xlabel("Radius (xSolar)")
plt.ylabel("Effective Temperature (K)")
plt.scatter(1.0, 5780,  color = '#FF6700', marker = "o", label = "Sun")
plt.scatter(radius, Teff, marker = "o", label = "Solar Twins")
plt.errorbar(radius, Teff, yerr = err_Teff, fmt = 'none', linewidth = 0.5, color = "black")
#age_abundance_fit = np.polyfit(radius, Teff, deg = 1)
plt.legend()

C_ab = [-0.05, -0.11, -0.02, 0.02, 0.03] #elemental C (C/H) abundance from Brewer 2016; HIP42333 = HD 73350 not in Brewer? IN SUN: log(8.39/12)

plt.figure()
plt.title("$^{13}$CO vs. Elemental C Abundance")
plt.xlabel("C Abundance = C/H (dex)")
plt.ylabel("Calculated $^{13}$CO Abundance (dex)")
plt.scatter(0.0, 0.0,  color = '#FF6700', marker = "o", label = "Sun")
plt.scatter(C_ab, [fit_min_dex1, fit_min_dex2, fit_min_dex4, fit_min_dex5, fit_min_dex6], marker = "o", label = "Solar Twins")
plt.errorbar(C_ab, [fit_min_dex1, fit_min_dex2, fit_min_dex4, fit_min_dex5, fit_min_dex6], yerr = [ab_err_dex1, ab_err_dex2, ab_err_dex4, ab_err_dex5, ab_err_dex6], fmt = 'none', linewidth = 0.5, color = "black")
#age_abundance_fit = np.polyfit(C_ab, [fit_min_dex1, fit_min_dex2, fit_min_dex4, fit_min_dex5, fit_min_dex6], deg = 1)
plt.legend()

O_ab = [-0.01, -0.05, -0.02, 0.03, -0.03] #elemental C (C/H) abundance from Brewer 2016, HIP42333 = HD 73350 not in Brewer? IN SUN: log(8.66/12)

plt.figure()
plt.title("$^{13}$CO vs. Elemental O Abundance")
plt.xlabel("O Abundance = O/H (dex)")
plt.ylabel("Calculated $^{13}$CO Abundance (dex)")
plt.scatter(0.0, 0.0,  color = '#FF6700', marker = "o", label = "Sun")
plt.scatter(O_ab, [fit_min_dex1, fit_min_dex2, fit_min_dex4, fit_min_dex5, fit_min_dex6], marker = "o", label = "Solar Twins")
plt.errorbar(O_ab, [fit_min_dex1, fit_min_dex2, fit_min_dex4, fit_min_dex5, fit_min_dex6], yerr = [ab_err_dex1, ab_err_dex2, ab_err_dex4, ab_err_dex5, ab_err_dex6], fmt = 'none', linewidth = 0.5, color = "black")
#age_abundance_fit = np.polyfit(O_ab, [fit_min_dex1, fit_min_dex2, fit_min_dex4, fit_min_dex5, fit_min_dex6], deg = 1)
plt.legend()