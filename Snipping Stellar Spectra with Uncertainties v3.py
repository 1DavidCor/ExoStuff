# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 14:11:18 2020

@author: David
"""

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import optimize

snip_width = 0.0006
chi_sqr_snipwidth = 10
CO_species = 13

#base directory
base_dir = 'C:\\Users\\d338c921\\Desktop\\M-band_Spectra_of_Solar_Twins\\'
#base_dir_solar_models = 'C:\\Users\\there\\Desktop\\M-band_Spectra_of_Solar_Twins\\Solar_Models\\'

obj_list = []

for root, dirs, files in os.walk(base_dir):
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

#snipping functions
def snip_simp(wavegrid, spectrum, linelocation, width):
    index = np.logical_and(wavegrid > (linelocation - width/2), wavegrid < (linelocation + width/2))
    return wavegrid[index], spectrum[index]
def snip(wavegrid, spectrum, error, linelocation, width):
    index = np.logical_and(wavegrid > (linelocation - width/2), wavegrid < (linelocation + width/2))
    return wavegrid[index], spectrum[index], error[index]

def useable_lines(line_list, skip_list):
    new_list = []
    for i in range(len(line_list)):
        if i not in skip_list:
            new_list = np.append(new_list, line_list[i])
    return new_list

def snip_spectrum(line_list, spectrum, error, wavegrid, snip_width, interp_sample_size): #indices_per_width uniform in solar spectrum but NOT in solar twin spectra
    snips_wl = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_flux = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_flux_norm = np.ndarray(shape = (len(line_list), interp_sample_size))
    snips_err = np.ndarray(shape = (len(line_list), interp_sample_size))
    for i in range(len(line_list)):
        wavegrid_mini = np.linspace(-snip_width/2, snip_width/2, interp_sample_size)
        wl_mini , spec_mini, err_mini = snip(wavegrid, spectrum, error, line_list[i], snip_width)
        spec_mini_interp = np.interp(wavegrid_mini, wl_mini - line_list[i], spec_mini)
        err_mini_interp = np.interp(wavegrid_mini, wl_mini - line_list[i], err_mini)
        spec_norm = spec_mini_interp / np.nanmean(spec_mini_interp) #normalize each snip to flux ~1
        snips_wl[i, :] = np.linspace(line_list[i] - (snip_width/2), line_list[i] + (snip_width/2), interp_sample_size)
        snips_flux[i, :] = spec_mini_interp
        snips_flux_norm[i, :] = spec_norm
        snips_err[i, :] = err_mini_interp / np.nanmean(spec_mini_interp)
    return snips_wl.T, snips_flux.T, snips_flux_norm.T, snips_err.T

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

#Write a function to output stack_vel and stack_flux for a given star/model i.e. same result as snip_plots with output_stack = True but without the plots
def stack_data(star, useable_line_list, snip_width):
#select a star
#     star:
#         1 = HIP 102040
#         2 = HIP 29432
#         3 = HIP 42333
#         4 = HIP 77052
#         5 = HIP 79672
#         6 = HIP 85042
#         96 = Solar Model (T = 5780 K; 0x)
#         97 = Solar Model (T = 5780 K; 0.3x)
#         98 = Solar Model (T = 5780 K; 1x)
#         99 = Solar Model (T = 5780 K; 3x)
        
#     useable_line_list: list of line locations for desired CO species; function will return an error if lines in NAN region are input
        
#     snip_width: width, in microns, of the spectrum we'd like to cut out at each line center
    
#     CO_species: integer input used for plot labels
#         12 = 12CO
#         13 = 13CO
#         18 = C18O


    if(star == 96):
        sun = pyfits.getdata(base_dir + "Solar_Models\\" + "nlte5780_4.44.0x.fromsun.hires.spec.fits")
        wl = sun[0,:] / 10000
        sun_flux_log = sun[1,:] #log base 10 of flux
        sun_flux = 10**sun_flux_log
        flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?

    elif(star == 97):
        sun = pyfits.getdata(base_dir + "Solar_Models\\" + "nlte5780_4.44.0.3x.fromsun.hires.spec.fits")
        wl = sun[0,:] / 10000
        sun_flux_log = sun[1,:] #log base 10 of flux
        sun_flux = 10**sun_flux_log
        flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?

    elif(star == 98):
        sun = pyfits.getdata(base_dir + "Solar_Models\\" + "nlte5780_4.44.1x.fromsun.hires.spec.fits")
        wl = sun[0,:] / 10000
        sun_flux_log = sun[1,:] #log base 10 of flux
        sun_flux = 10**sun_flux_log
        flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?

    elif(star == 99):
        sun = pyfits.getdata(base_dir + "Solar_Models\\" + "nlte5780_4.44.3x.fromsun.hires.spec.fits")
        wl = sun[0,:] / 10000
        sun_flux_log = sun[1,:] #log base 10 of flux
        sun_flux = 10**sun_flux_log
        flux = sun_flux / np.nanmean(sun_flux) #nanmean or nanmedian?
        
    else:
    #star = 1, 2, 3, 4, 5, 6
        choose_file = int(star - 1)
        obj_path = base_dir + obj_list[choose_file]
        spect = pyfits.getdata(obj_path)
        obj_name = pyfits.getval(obj_path, 'OBJECT')
        
        #create variables for wavelength, flux, and uncertainty; fits files 1, 2, and 3 respectively
        wl = spect[0,:]
        flux = spect[1,:]
        err = spect[2,:]
        
        #normalize flux and wl
        wl = wl / (beta_shifts[choose_file] + 1)
        err = err/np.nanmean(flux)
        flux = flux/np.nanmean(flux) #nanmean or nanmedian?
        
        
    
    if (star in [96, 97, 98, 99]):
        #Solar Models
        snips_wl, snips_flux, snips_flux_norm = snip_spectrum_simp(useable_line_list, flux, wl, snip_width, 50)
        stack_vel = np.mean((snips_wl - useable_line_list) / useable_line_list * 3e5, axis = 1)
        stack_flux = np.mean(snips_flux_norm, axis = 1)
        stack_flux = stack_flux - (np.max(stack_flux) - 1) #remember to lower models so they max at ~1
        
        return stack_vel, stack_flux
    
    else:
        snips_wl, snips_flux, snips_flux_norm, snips_err = snip_spectrum(useable_line_list, flux, err, wl, snip_width, 50)
        #stars w/ velocity shift corrections
        stack_vel = np.mean((snips_wl - useable_line_list) / useable_line_list * 3e5, axis = 1) + velocity_shifts[star - 1]
        stack_flux = np.average(snips_flux_norm, axis = 1, weights = snips_err**-2) ###weighted average using weights 1/variance
        stack_err = np.sqrt(1./ (snips_err**-2).sum(axis = 1))
        
        return stack_vel, stack_flux, stack_err
    
#Write a fuction to find this line for a given stacked absorption line; outputs np.poly1d(fit_coefficients) to be divided from the spectrum later
def de_slant(stack_vel, stack_flux):
    #returns a np.poly1d polynomial (a line) to be divided out from the spectrum in order to correct the slant
    
    vel_snipL, flux_snipL = snip_simp(stack_vel, stack_flux, -7.5, 5) #snips the left "shoulder" of the stacked absorption line i.e. [-10, -5] km/s
    vel_snipR, flux_snipR = snip_simp(stack_vel, stack_flux, 7.5, 5) #snips the right "shoulder" of the stacked absorption line i.e. [-10, -5] km/s
    
    avgL = np.mean(flux_snipL) #average flux of the left shoulder
    avgR = np.mean(flux_snipR) #average flux of the right shoulder
    
    l_coeff = np.polyfit([-7.5, 7.5], [avgL, avgR], deg = 1) #fits a line through the midpoints of each shoulder i.e. (- 7.5 km/s, avgL) and ( 7.5 km/s, avgR); returns coeffecients 
    line = np.poly1d(l_coeff) #puts the equation of the line into a "useable form"; divide out of spectrum as line(stack_vel)
    
    return line

#Write a function to calculate chi-squared values
def chi_sqr(stack_vel_obs, stack_vel_model, stack_flux_obs, stack_flux_model, err):
    #resize arrays to exclude the wings
    stack_vel_obs, stack_flux_obs, err = snip(stack_vel_obs, stack_flux_obs, err, 0, chi_sqr_snipwidth)
    stack_vel_model, stack_flux_model = snip_simp(stack_vel_model, stack_flux_model, 0, chi_sqr_snipwidth)
    chisqr = 0
    for i in range(len(stack_flux_obs) - 1):
        chisqr = chisqr + ((stack_flux_obs[i] - stack_flux_model[i])**2)/(err[i]**2)
    return chisqr

def calc_abundance(star_num, line_list, solar_line_list, star_skiplist, snip_width, abundance_guess, plot, label = "HIP ###"):
    ###requires stack_data() and chi_sqr() functions
    #requires scipy.optimize
    
    #star_num = a number 1-6; used to pick which star we're looking at
    #line_list = list of lines to skip in stacked absorption line analysis
    #how wide to make the snips; default = 0.0006
    #abundance_guess = number provided for finding the chi-sqr minimum i.e. your best guess at what the xSolar abundance is
    #if plot = true, show the plot of chi-sqr values, data fit, and min
    
    #stack_data for chosen star and all four solar models
    star_lines = useable_lines(line_list, star_skiplist)
    sun_lines = useable_lines(solar_line_list, star_skiplist)
    stack_vel, stack_flux, stack_err = stack_data(star_num, star_lines, snip_width)
    
    stack_vel1, stack_flux1 = stack_data(96, sun_lines, snip_width)
    stack_vel2, stack_flux2 = stack_data(97, sun_lines, snip_width)
    stack_vel3, stack_flux3 = stack_data(98, sun_lines, snip_width)
    stack_vel4, stack_flux4 = stack_data(99, sun_lines, snip_width)
        
    #calculate 4 chi_sqr values and put them in an array
    chi_sqr_list = [chi_sqr(stack_vel, stack_vel1, stack_flux, stack_flux1, stack_err), chi_sqr(stack_vel, stack_vel2, stack_flux, stack_flux2, stack_err), chi_sqr(stack_vel, stack_vel3, stack_flux, stack_flux3, stack_err), chi_sqr(stack_vel, stack_vel4, stack_flux, stack_flux4, stack_err)]
    
    #print(chi_sqr_list)
    #create a matching "x-array" of xSolar abundances i.e. [0, 1/3, 1, 3]
    x_list = [0.0, 1/3, 1.0, 3.0]
    
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
        plt.title("$\chi^2$ Fit: " + label)
        plt.xlabel("xSolar Abundance")
        plt.ylabel("$\chi^2$")
        plt.scatter(x_list, chi_sqr_list, marker = 'o', label = "$\chi^2$ Values")
        plt.scatter(err_bounds, parabola(err_bounds), color = "r", marker = "x")
        plt.plot(np.linspace(0, 3, 1000), parabola(np.linspace(0, 3, 1000)), label = "Data Fit")
        plt.plot(fit_min, parabola(fit_min), marker = "X", label = ("min @ x = %s" % fit_min)) #different color than previous
        plt.legend()
        
    else:
        pass
    
    return np.round(fit_min, 2), np.round(fit_min_dex, 2), np.round(abundance_err, 2), np.round(abundance_err_dex[0], 2)


#Write a function to automate the abundance plots i.e. stellar stacked absorption lines + 4 model stacked lines
def abundance_plot(star_num, line_list, solar_line_list, star_skiplist, label = "HIP ###"):
    star_lines = useable_lines(line_list, star_skiplist)
    sun_lines = useable_lines(solar_line_list, star_skiplist)
    
    stack_vel, stack_flux, stack_err = stack_data(star_num, star_lines, snip_width)
    slant = de_slant(stack_vel, stack_flux)
    
    stack_vel_mod1, stack_flux_mod1 = stack_data(96, sun_lines, snip_width)
    stack_vel_mod2, stack_flux_mod2 = stack_data(97, sun_lines, snip_width)
    stack_vel_mod3, stack_flux_mod3 = stack_data(98, sun_lines, snip_width)
    stack_vel_mod4, stack_flux_mod4 = stack_data(99, sun_lines, snip_width)
    
    plt.figure()
    plt.ylim(0.91, 1.02)
    plt.xlim(-15, 15)
    plt.plot(stack_vel, stack_flux / slant(stack_vel), label = label, linewidth = 3)
    plt.plot(stack_vel_mod1, stack_flux_mod1, linestyle='--', label = "0x Solar $^{13}$CO")
    plt.plot(stack_vel_mod2, stack_flux_mod2, linestyle='--', label = "1/3x Solar $^{13}$CO")
    plt.plot(stack_vel_mod3, stack_flux_mod3, linestyle='--', label = "Solar $^{13}$CO")
    plt.plot(stack_vel_mod4, stack_flux_mod4, linestyle='--', label = "3x Solar $^{13}$CO")
    plt.errorbar(stack_vel, stack_flux / slant(stack_vel), yerr = stack_err / slant(stack_vel), fmt = 'none', linewidth = 0.5, color = "black")
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Normalized Flux Intensity")
    plt.axvline(x = 0, color= 'k', linestyle='--')
    plt.legend()
    return

velocity_shifts = [1.97337, 0.944558, 0.600705 - 0.425, 0.94817 - 0.455, 0.772493 - 0.88, 0.943784 - 0.845] #add to velocity grid of stacked absorption lines; + moves spectrum right, - moves left

#plot for sun using solar skip_list
solar_skip_list = [0, 1, 2, 6, 12, 14, 15, 16, 19, 21, 23, 25, 26, 27, 29, 32, 33]
sun_13_CO_lines = useable_lines(wl_13_CO, solar_skip_list)


#star 1: HIP 102040
star1_skiplist = [8, 11, 12, 24, 25]
abundance_plot(1, wl_13_CO, sun_13_CO_lines, star1_skiplist, label = "HIP 102040")

fit_min1, fit_min_dex1, ab_err1, ab_err_dex1 = calc_abundance(1, wl_13_CO, sun_13_CO_lines, star1_skiplist, snip_width, 1.0, plot=True, label = "HIP 102040")
print("HIP 102040 has a 13CO abundance of " + str(fit_min1) + " xSolar or " + str(fit_min_dex1) + " dex \n")


#star 2: HIP 29432
star2_skiplist = [8, 10, 11, 12, 15, 18, 19, 23, 24, 25, 27, 32]
abundance_plot(2, wl_13_CO, sun_13_CO_lines, star2_skiplist, label = "HIP 29432")

fit_min2, fit_min_dex2, ab_err2, ab_err_dex2 = calc_abundance(2, wl_13_CO, sun_13_CO_lines, star2_skiplist, snip_width, 1.55, plot=True, label = "HIP 29432")
print("HIP 29432 has a 13CO abundance of " + str(fit_min2) + " xSolar or " + str(fit_min_dex2) + " dex \n")


#star 3: HIP 42333
star3_skiplist = [11, 12, 24, 25]
abundance_plot(3, wl_13_CO, sun_13_CO_lines, star3_skiplist, label = "HIP 42333")

fit_min3, fit_min_dex3, ab_err3, ab_err_dex3 = calc_abundance(3, wl_13_CO, sun_13_CO_lines, star3_skiplist, snip_width, 1.5, plot=True, label = "HIP 42333")
print("HIP 42333 has a 13CO abundance of " + str(fit_min3) + " xSolar or " + str(fit_min_dex3) + " dex \n")


#star 4: HIP 77052
star4_skiplist = [2, 11, 12, 23, 24, 25]
star4_useable_sun_lines = useable_lines(sun_13_CO_lines, star4_skiplist)
abundance_plot(4, wl_13_CO, sun_13_CO_lines, star4_skiplist, label = "HIP 77052")

fit_min4, fit_min_dex4, ab_err4, ab_err_dex4 = calc_abundance(4, wl_13_CO, sun_13_CO_lines, star4_skiplist, snip_width, 1.5, plot=True, label = "HIP 77052")
print("HIP 77052 has a 13CO abundance of " + str(fit_min4) + " xSolar or " + str(fit_min_dex4) + " dex \n")

#star 5: HIP 79672
star5_skiplist = [2, 11, 12, 15, 23, 24, 25, 33]
star5_useable_sun_lines = useable_lines(sun_13_CO_lines, star5_skiplist)
abundance_plot(5, wl_13_CO, sun_13_CO_lines, star5_skiplist, label = "HIP 79672")


fit_min5, fit_min_dex5, ab_err5, ab_err_dex5 = calc_abundance(5, wl_13_CO, sun_13_CO_lines, star5_skiplist, snip_width, 0.8, plot=True, label = "HIP 79672")
print("HIP 79672 has a 13CO abundance of " + str(fit_min5) + " xSolar or " + str(fit_min_dex5) + " dex \n")


#star 6: HIP 85042
star6_skiplist = [8, 11, 12, 24, 25]
star6_useable_sun_lines = useable_lines(sun_13_CO_lines, star6_skiplist)
abundance_plot(6, wl_13_CO, sun_13_CO_lines, star6_skiplist, label = "HIP 85042")

fit_min6, fit_min_dex6, ab_err6, ab_err_dex6 = calc_abundance(6, wl_13_CO, sun_13_CO_lines, star6_skiplist, snip_width, 1.7, plot=True, label = "HIP 85042")
print("HIP 85042 has a 13CO abundance of " + str(fit_min6) + " +/- " + str(ab_err6) + " xSolar or " + str(fit_min_dex6) + " dex \n")

plt.figure()
plt.title("$^{13}$CO Abundance vs. Stellar Age")
plt.xlabel("Stellar Age (Gyr)")
plt.ylabel("Calculated Abundance (dex)")
plt.scatter([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], marker = "o")
plt.errorbar([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], yerr = [ab_err_dex1, ab_err_dex2, ab_err_dex3, ab_err_dex4, ab_err_dex5, ab_err_dex6], fmt = 'none', linewidth = 0.5, color = "black")
plt.errorbar([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], xerr = [0.91, 0.71, 0.52, 0.91, 0.39, 0.62], fmt = 'none', linewidth = 0.5, color = "black")
age_abundance_fit = np.polyfit([2.42, 5.51, 1.01, 3.67, 3.09, 6.66], [fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6], deg = 1)
#age_abundance_line = np.poly1d(age_abundance_fit)
plt.plot(np.linspace(0, 7, 1000), age_abundance_fit[0]*(np.linspace(0, 7, 1000)) + age_abundance_fit[1])

print("xSolar abundances: ")
print(fit_min1, fit_min2, fit_min3, fit_min4, fit_min5, fit_min6)
print("Uncertainties: ")
print(ab_err1, ab_err2, ab_err3, ab_err4, ab_err5, ab_err6)

print("Abundances (dex): ")
print(fit_min_dex1, fit_min_dex2, fit_min_dex3, fit_min_dex4, fit_min_dex5, fit_min_dex6)
print("Uncertainties: ")
print(ab_err_dex1, ab_err_dex2, ab_err_dex3, ab_err_dex4, ab_err_dex5, ab_err_dex6)
