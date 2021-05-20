# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 21:55:20 2020

@author: David
"""
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import os

#base directory
base_dir = 'C:\\Users\\there\\Desktop\\M-band_Spectra_of_Solar_Twins\\'

obj_list = []

for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith('final.fits'):
            obj_list.append(file)

#enter an integer to determine what file to load
# ftl = input("Please enter an integer to decide with spectrum to load:\n 1) HIP 29432\n 2) HIP 42333\n 3) HIP 77052\n 4) HIP 79672\n 5) HIP 85042\n 6) HIP 102040\n")
# ftl = int(ftl)

# if (ftl == 1):
#     spect = pyfits.getdata(obj_path + 'hip29432_m1_final.fits')
#     obj = "HIP 29432"
# elif(ftl == 2):
#     spect = pyfits.getdata(obj_path + 'hip42333_m1_final.fits')
#     obj = "HIP 42333"
# elif(ftl == 3):
#     spect = pyfits.getdata(obj_path + 'hip77052_m1_final.fits')
#     obj = "HIP 77052"
# elif(ftl == 4):
#     spect = pyfits.getdata(obj_path + 'hip79672_m1_final.fits')
#     obj = "HIP 79672"
# elif(ftl == 5):
#     spect = pyfits.getdata(obj_path + 'hip85042_m1_final.fits')
#     obj = "HIP 85042"
# else:
#     spect = pyfits.getdata(obj_path + 'hip102040_m1_final.fits')
#     obj = "HIP 102040"

#choose file from obj_list to load
choose_file = int(input("Choose a file to load:\n 1) HIP 102040\n 2) HIP 29432\n 3) HIP 42333\n 4) HIP 77052\n 5) HIP 79672\n 6) HIP 85042\n OR\n 99) Sun\n")) - 1

##if/else to to choose solar spectrum instead of solar_twin
if (choose_file == 98):
    spect = pyfits.getdata(base_dir + "nlte5780_4.44.1x.fromsun.hires.spec.fits")
    wl = spect[0,:]
    flux_log = spect[1,:]
    flux = 10**flux_log     #log base 10 of flux
else:
    obj_path = base_dir + obj_list[choose_file]
    spect = pyfits.getdata(obj_path)
    obj_name = pyfits.getval(obj_path, 'OBJECT')
    #create variables for wavelength, flux, and uncertainty; fits files 1, 2, and 3 respectively
    wl = spect[0,:]
    flux = spect[1,:]
    err = spect[2, :]

#check spect type and dimensions
print(type(spect))
print(spect.shape)

CO_rest_lines_12 = [45803.8234711, 45876.3599396, 46024.4369507, 46099.98703, 46176.5575409, 46412.4202728, 46493.1201935, 46574.8596191, 46826.4293671, 46999.4974136, 45821.8812943, 45887.4750137, 46090.1927948, 46159.7537994, 46301.8846512, 46448.0733871, 46522.693634, 46752.7675629, 46831.536293, 46992.2494888, 46079.9837112, 46141.5290833, 46204.0519714, 46267.5666809, 46397.5763321, 46464.0855789, 46531.5961838, 46669.6548462, 46884.4127655, 46958.0602646]

CO_rest_lines_13 = [46116.476, 46178.775, 46205.730, 46241.999, 46317.954, 46404.576, 46433.864, 46556.783, 46641.064, 46710.901, 46717.315, 46926.227, 46935.010]

#shift the wl according to Doppler Shift eq.
beta_shifts = [-1.9263267897412533e-04, -1.1545041431507331e-04, 1.9379144076145112e-04, 6.888621612746786e-05, 2.4225025924415754e-05, -1.5905756391913187e-04] 

#convert list to array and convert from Angstroms to microns
CO_rest_lines_12 = np.array(CO_rest_lines_12) / 10000
CO_rest_lines_13 = np.array(CO_rest_lines_13) / 10000

#create a plot
#if else for solar twin plot vs. Sun plot
# if (choose_file == 98):
#     plt.figure()
#     plt.title("Sun Reduced/Calibrated Spectrum")
#     plt.xlabel("Wavelength (microns)")
#     plt.ylabel("Flux (???)")
#     #plt.ylim(0, 1e-13)
#     plt.xlim(4.6,4.7)
#     plt.plot(wl/10000, flux)
#     #plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
#     for i in range(13):
#         plt.axvline(x = CO_rest_lines_13[i], color= 'k', linestyle=':')
# else:
#     plt.figure()
#     plt.title(obj_name + " Reduced/Calibrated Spectrum")
#     plt.xlabel("Wavelength (microns)")
#     plt.ylabel("Flux (???)")
#     plt.ylim(0, 1e-13)
#     plt.xlim(4.6,4.7)
#     plt.plot(wl, flux)
#     plt.plot(wl / (beta_shifts[choose_file] + 1), flux)
#     #plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
#     for i in range(30):
#         plt.axvline(x = CO_rest_lines_12[i], color= 'k', linestyle=':')

#look at the plot to determine shift; calculate b
# lam_obs = 4.644725
# lam_rest = 12_CO_rest_lines[15]
# b1 = (lam_obs - lam_rest)/lam_rest - 0.0000150
# b2 = (lam_obs - lam_rest)/lam_rest - 0.00006
# b3 = (lam_obs - lam_rest)/lam_rest - 0.00015 #change this value for fine tuning: - shifts right, + shifts left 
# plt.plot(wl/(b1 + 1), flux)
# plt.plot(wl/(b2 + 1), flux)
# plt.plot(wl/(b3 + 1), flux)

# beta_shifts = [8.063560340486664e-06, 4.378369209256946e-06, 0, 6.888621612746786e-05, 2.4225025924415754e-05, -1.9847063134193604e-05]

if not (choose_file == 98):
    #SNR = flux / err ? Or is it more ecomplicated than that?
    snr = flux / err
    #calcaulate the snr from 4.6 to 4.7 microns
    ## index at 4.6 microns
    min_range_index = np.array(np.where((wl >= 4.6) & (wl <= 4.7)))[0,0]
    ##index at 4.7 microns
    max_range_index = np.array(np.where((wl >= 4.6) & (wl <= 4.7)))[0,-1]
    ##array of snr over the given range
    snr_range = flux[min_range_index:max_range_index] / err[min_range_index:max_range_index]
    #average snr over the range to be displayed in plot
    snr_range_avg = np.nansum(snr_range) / (max_range_index - min_range_index)

#calculate the average by summing the snr at each wl and dividing by # entries?
#snr_avg = np.nansum(snr) / np.size(snr) #np.nansum treats nan = 0
### np.sum = nan; the spectrum is discontinuous? 
###how would we average with this in mind? 

##Why is the spectrum centered around 0.5??
##What spectral features do we look for in the infrared?
###In optical ranges, we looked for H-alpha, Calcium triplet, and K-lines
##How do we distinguish between starlight and hot pixels/cosmic rays? 
###We shouldn't get many cosmic rays since exposure times are short, right?

if not (choose_file == 98):
    plt.figure()
    plt.title(obj_name + "  SNR (4.6 to 4.7) = " + str("%.3f" % snr_range_avg))
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("SNR")
    #plt.ylim(0, 1e-13)
    plt.xlim(4.6,4.7)
    plt.plot(wl, snr)

#make plots pretty!!
rainbow = ['#ff0000', '#ffa500', '#ffff00', '#008000','#10a5f5', '#673ab7', '#ee82ee']

#stacked plot of all 6 spectra + solar spectra
# plt.figure()
# plt.title("Reduced/Calibrated Spectra: Stacked")
# plt.xlabel("Wavelength (microns)")
# plt.ylabel("Flux (???)")
# plt.ylim(0, 1e-13)
# plt.xlim(4.6,4.7)
# for i in range(6):
#     spect_stack = pyfits.getdata(base_dir + obj_list[i])
#     wl_stack = spect_stack[0,:]
#     flux_stack = spect_stack[1,:]
#     plt.plot(wl_stack, flux_stack, color = rainbow[i])
# #plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
# for i in range(30):
#     plt.axvline(x = CO_rest_lines_12[i], color= 'k', linestyle=':')
    
#stacked plot of all 6 spectra; wl corrected
# plt.figure()
# plt.title("Reduced/Calibrated Spectra: Corrected")
# plt.xlabel("Wavelength (microns)")
# plt.ylabel("Flux (???)")
# plt.ylim(0, 1e-13)
# plt.xlim(4.6,4.7)
# for i in range(6):
#     spect_stack = pyfits.getdata(base_dir + obj_list[i])
#     wl_stack = spect_stack[0,:]
#     flux_stack = spect_stack[1,:]
#     plt.plot(wl_stack / (beta_shifts[i] + 1), flux_stack, color = rainbow[i])
# #plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
# for i in range(30):
#     plt.axvline(x = CO_rest_lines_12[i], color= 'k', linestyle=':')
    
#stacked plot of all spectra; wl corrected; normalized; waterfalled
plt.figure()
#plt.title("ALL Reduced/Calibrated Spectra: Corrected")
#plt.xlabel("Wavelength (microns)")
#plt.ylabel("Flux (???)")
plt.ylim(0,12)
plt.xlim(4.6,4.7)
for i in range(6):
    spect_stack = pyfits.getdata(base_dir + obj_list[i])
    wl_stack = spect_stack[0,:]
    flux_stack = spect_stack[1,:]
    flux_norm = spect_stack[1,:] / np.nanmedian(spect_stack[1,:]) + 6 - i
    plt.plot(wl_stack / (beta_shifts[i] + 1), flux_norm, color = rainbow[i])
#plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
for i in range(30):
    plt.axvline(x = CO_rest_lines_12[i], color= 'k', linestyle=':')
for i in range(13):
    plt.axvline(x = CO_rest_lines_13[i], color= 'k', linestyle='--')
#normalize sun
sun = pyfits.getdata(base_dir + "nlte5780_4.44.1x.fromsun.hires.spec.fits")
sun_wl = sun[0,:] / 10000
sun_flux_log = sun[1,:] #log base 10 of flux
sun_flux = 10**sun_flux_log
sun_flux_norm = sun_flux / np.nanmedian(sun_flux)
for i in range(6):
    plt.plot(sun_wl, sun_flux_norm + 7 - i, color = rainbow[-1])
    

#FIND CO LINES!!!
from hapi import *
db_begin("hapi_data")

from astropy import units as u

#M = HITRAN Molecule Number
#I = HITRAN Isotopologue Number
#12CO has ISO_ID = 26; M = 5, I = 1
#13CO has ISO_ID = 27; M = 5, I = 2
#C18O has ISO_ID = 28; M = 5, I = 3
#wavenumber = wavelength**-1 in cm**-1

nu_min = (4.7e-6 * u.m).to(u.cm)**(-1)
nu_max = (4.6e-6 * u.m).to(u.cm)**(-1)

#fetch('HOH',1,1,nu_min.value,nu_max.value)
#Get columns from table: p1,p2,p3 = getColumns('sampletab',('p1','p2','p3'))
#We want nu (wavenumber) and sw (intensity)

#12CO
fetch('CO', 5, 1, nu_min.value, nu_max.value)
nu_12_CO, sw_12_CO = getColumns('CO', ['nu', 'sw'])
print(nu_12_CO.shape)
#13CO
fetch('CO', 5, 2, nu_min.value, nu_max.value)
nu_13_CO, sw_13_CO = getColumns('CO', ['nu', 'sw'])
print(nu_13_CO.shape)
#C18O
fetch('CO', 5, 3, nu_min.value, nu_max.value)
nu_C_18_O, sw_C_18_O = getColumns('CO', ['nu', 'sw'])
print(nu_C_18_O.shape)

print("\n" + str(nu_12_CO.shape[0]) +  " 12_CO Lines from 4.6 to 4.7 microns (cm^-1):\n")
print(np.sort(nu_12_CO))
print("\n" + str(nu_12_CO.shape[0]) +  " 13_CO Lines from 4.6 to 4.7 microns (cm^-1):\n")
print(np.sort(nu_13_CO))
print("\n" + str(nu_C_18_O.shape[0]) +  " C_18_O Lines from 4.6 to 4.7 microns (cm^-1):\n")
print(np.sort(nu_C_18_O))

#convert from wavenumber to wavelength in microns
wl_12_CO = nu_12_CO**(-1) * 10000
wl_12_CO = np.sort(wl_12_CO)
wl_13_CO = nu_13_CO**(-1) * 10000
wl_13_CO = np.sort(wl_13_CO)
wl_C_18_O = nu_C_18_O**(-1) * 10000
wl_C_18_O = np.sort(wl_C_18_O)

print("\n" + str(wl_12_CO.shape[0]) +  " 12_CO Lines from 4.6 to 4.7 microns (microns):\n")
print(wl_12_CO)
print("\n" + str(wl_13_CO.shape[0]) +  " 13_CO Lines from 4.6 to 4.7 microns (microns):\n")
print(wl_13_CO)
print("\n" + str(wl_C_18_O.shape[0]) +  " C_18_O Lines from 4.6 to 4.7 microns (microns):\n")
print(wl_C_18_O)

# #FIND LINES OBSCURED BY TELLURICS
# #We want lines located where telluric transmission ~1

# #load telluric data file; Column 1: wavelength in microns; Column 2: fractional telluric transmission
# tellurics = pyfits.getdata(base_dir + "transdata_0,5-14_mic_hires.fits")
# wl_tellurics = tellurics[:, 0] #in microns
# trans_tellurics = tellurics[:, 1]

# plt.figure()
# plt.xlim(4.6, 4.7)
# plt.plot(wl_tellurics, trans_tellurics)
# plt.title("Tellurics")
# plt.xlabel("Wavelength (microns)")
# plt.ylabel("Flux (???)")
# for i in range(wl_12_CO.shape[0] - 1):
#     plt.axvline(x = wl_12_CO[i], color= 'k', linestyle=':')
# for i in range(wl_13_CO.shape[0] - 1):
#     plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
# for i in range(wl_C_18_O.shape[0] - 1):
#     plt.axvline(x = wl_C_18_O[i], color= 'k', linestyle='-')
# label1 = plt.axvline(x = wl_12_CO[-1], color= 'k', linestyle=':')
# label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--')
# label3 = plt.axvline(x = wl_C_18_O[-1], color= 'k', linestyle='-')
# plt.legend([label1,label2,label3], ["12_CO","13_CO","C_18_O"])

# obscured_lines = [wl_12_CO[2], wl_12_CO[5], wl_12_CO[8], wl_12_CO[11], wl_12_CO[14], wl_12_CO[15], wl_12_CO[18], wl_12_CO[21], wl_12_CO[24], wl_12_CO[29], wl_12_CO[30], wl_12_CO[32], wl_12_CO[33], wl_12_CO[36], wl_12_CO[37], wl_12_CO[38], wl_12_CO[39], wl_12_CO[40], wl_13_CO[5], wl_13_CO[8], wl_13_CO[15], wl_13_CO[23], wl_13_CO[32], wl_13_CO[33], wl_C_18_O[0], wl_C_18_O[11], wl_C_18_O[30], wl_C_18_O[31]]

# print("\n Obscured Lines:\n")
# print(obscured_lines)

# #FIND 12_CO and 13_CO lines that are too close together
# plt.figure()
# plt.xlim(4.6, 4.7)
# plt.title("Sun Reduced/Calibrated Spectrum")
# plt.xlabel("Wavelength (microns)")
# plt.ylabel("Flux (???)")
# plt.plot(sun_wl, sun_flux_norm)
# for i in range(wl_12_CO.shape[0] - 1):
#     plt.axvline(x = wl_12_CO[i], color= 'k', linestyle=':')
# for i in range(wl_13_CO.shape[0] - 1):
#     plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
# label1 = plt.axvline(x = wl_12_CO[-1], color= 'k', linestyle=':')
# label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--')
# plt.legend([label1,label2], ["12_CO","13_CO"])


#normalize sun
sun = pyfits.getdata(base_dir + "nlte5780_4.44.1x.fromsun.hires.spec.fits")
sun_wl = sun[0,:] / 10000
sun_flux_log = sun[1,:] #log base 10 of flux
sun_flux = 10**sun_flux_log
sun_flux_norm = sun_flux / np.nanmedian(sun_flux)

#PLOTS UPDATED WITH HITRAN LINE LISTS
if (choose_file == 98):
    plt.figure()
    plt.title("Solar Spectrum")
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("Flux (???)")
    #plt.ylim(0, 1e-13)
    plt.xlim(4.6,4.7)
    plt.plot(wl/10000, flux)
    #plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
    for i in range(wl_12_CO.shape[0] - 1):
        plt.axvline(x = wl_12_CO[i], color= 'k', linestyle=':')
    for i in range(wl_13_CO.shape[0] - 1):
        plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
    for i in range(wl_C_18_O.shape[0] - 1):
        plt.axvline(x = wl_C_18_O[i], color= 'k', linestyle='-')
    label1 = plt.axvline(x = wl_12_CO[-1], color= 'k', linestyle=':')
    label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--')
    label3 = plt.axvline(x = wl_C_18_O[-1], color= 'k', linestyle='-')
    plt.legend([label1,label2,label3], ["12_CO","13_CO","C_18_O"])
else:
    plt.figure()
    plt.title(obj_name + " Reduced/Calibrated Spectrum")
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("Flux (???)")
    plt.ylim(0.8, 1.5)
    plt.xlim(4.6,4.7)
    #plt.plot(sun_wl, sun_flux_norm + 0.9, label = "Sun")
    plt.plot(wl / (beta_shifts[choose_file] + 1), flux/np.nanmedian(flux), color = rainbow[choose_file], label = obj_name)
    plt.legend()
    #plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
    for i in range(wl_12_CO.shape[0] - 1):
        plt.axvline(x = wl_12_CO[i], color= 'k', linestyle=':')
    for i in range(wl_13_CO.shape[0] - 1):
        plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
    for i in range(wl_C_18_O.shape[0] - 1):
        plt.axvline(x = wl_C_18_O[i], color= 'k', linestyle='-')
    label1 = plt.axvline(x = wl_12_CO[-1], color= 'k', linestyle=':')
    label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--')
    label3 = plt.axvline(x = wl_C_18_O[-1], color= 'k', linestyle='-')
    plt.legend([label1,label2,label3], ["12_CO","13_CO","C_18_O"])


#updated stacked plot 
plt.figure()
plt.title("ALL Reduced/Calibrated Spectra: Corrected")
plt.xlabel("Wavelength (microns)")
plt.ylabel("Flux (???)")
plt.ylim(0,10)
plt.xlim(4.6,4.7)
for i in range(6):
    obj_name = pyfits.getval(base_dir + obj_list[i], 'OBJECT')
    spect_stack = pyfits.getdata(base_dir + obj_list[i])
    wl_stack = spect_stack[0,:]
    flux_stack = spect_stack[1,:]
    flux_norm = spect_stack[1,:] / np.nanmedian(spect_stack[1,:]) + 6 - i
    plt.plot(wl_stack / (beta_shifts[i] + 1), flux_norm, color = rainbow[i], label = obj_name)
    plt.legend()
#plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
for i in range(wl_12_CO.shape[0] - 1):
    plt.axvline(x = wl_12_CO[i], color= 'k', linestyle=':')
for i in range(wl_13_CO.shape[0] - 1):
    plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
for i in range(wl_C_18_O.shape[0] - 1):
    plt.axvline(x = wl_C_18_O[i], color= 'k', linestyle='-')
label1 = plt.axvline(x = wl_12_CO[-1], color= 'k', linestyle=':')
label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--')
label3 = plt.axvline(x = wl_C_18_O[-1], color= 'k', linestyle='-')
plt.legend([label1,label2,label3], ["12_CO","13_CO","C_18_O"])
for i in range(6):
    plt.plot(sun_wl, sun_flux_norm + 7 - i, color = rainbow[-1])
plt.legend()


#Are C18O lines visible in the solar spectrum?
plt.figure()
plt.title("Are the C_18_O lines visible in Solar Spectrum?")
plt.xlabel("Wavelength (microns)")
plt.ylabel("Flux (???)")
plt.xlim(4.6,4.7)
plt.ylim(0.2,0.5)
plt.plot(sun_wl, sun_flux_norm)
for i in range(wl_C_18_O.shape[0] - 1):
    plt.axvline(x = wl_C_18_O[i], color= 'k', linestyle='-')
label3 = plt.axvline(x = wl_C_18_O[-1], color= 'k', linestyle='-')
plt.legend([label3], ["C_18_O"])
#Yes, they're visible. NOTE: They don't all line up very well.

#Create a list of useable lines for each star:
obj0_useable_lines = []
#Rejects [20]: wl_12_CO[2] & wl_C_18_O[0] (cannot discern between the two), wl_C_18_O[6] (obscured by tellurics), wl_13_CO[8] & wl_12_CO[12] & wl_C_18_O[9] (obscured by tellurics), wl_13_CO[11] & wl_12_CO[16] & wl_13_CO[12] (spectrum discontinuous in this region), wl_12_CO[20] & wl_C_18_O[16] (obscured by tellurics), wl_12_CO[24] (obscured by tellurics), wl_C_18_O[19] (obscured by tellurics), wl_12_CO[31] (obscured by tellurics), wl_13_CO[24] & wl_C_18_O[26] & wl_13_CO[25] & wl_12_CO[32] (spectrum discontinuous in this region), wl_12_CO[40] & wl_C_18_O[33] (cannot discern between the two)

#Weak [9]: wl_C_18_O[5], wl_C_18_O[7], wl_C_18_O[13], wl_C_18_O[17], wl_C_18_O[18], wl_C_18_O[21], wl_C_18_O[25], wl_C_18_O[27], wl_C_18_O[28]

#Note: wl_C_18_O[15] is STRONG; wl_13_CO[23] & wl_12_CO[29] very close together BUT discernible in stellar spectrum, wl_C_18_O[31] discernible in stellar spectrum

obj1_useable_lines = []
#Rejects [22]: wl_12_CO[2] & wl_C_18_O[0] (cannot discern between the two), wl_C_18_O[1] (obscured by tellurics), wl_13_CO[10] (cannot discern), wl_13_CO[11] & wl_12_CO[16] & wl_13_CO[12] (spectrum discontinuous in this region), wl_C_18_O[15] (obscured by tellurics), wl_12_CO[21] & wl_13_CO[15] (cannot discern between the two), wl_13_CO[18] &  wl_12_CO[24] (obscured by tellurics), wl_13_CO[19] (obscured by tellurics), wl_13_CO[23] & wl_12_CO[29] (cannot discern between the two), wl_13_CO[24] & wl_C_18_O[26] & wl_13_CO[25] & wl_12_CO[32] (spectrum discontinuous in this region), wl_13_CO[27] (cannot discern), wl_C_18_O[30] & wl_13_CO[32] (part of a larger spectral feature; can't discern)

#Weak [8]: wl_C_18_O[3], wl_C_18_O[10], wl_C_18_O[14], wl_C_18_O[18], wl_C_18_O[25], wl_C_18_O[28], wl_C_18_O[29], wl_C_18_O[31]

#Note: Telluric region: 4.6250 - 4.670 & 4.6825-4.6835 (but lines still look good??), wl_C_18_O[11] & wl_12_CO[14] & wl_12_CO[14] very close together BUT discernible in stellar spectrum, wl_12_CO[40] & wl_13_CO[33] very close together BUT discernible in stellar spectrum



#Define a function to snip parts out of the spectrum
def snip(wavegrid, spectrum, linelocation, width):
    index = np.logical_and(wavegrid > (linelocation - width/2), wavegrid < (linelocation + width/2))
    return wavegrid[index], spectrum[index]

#snip_wl, snip_spect = snip_spectrum(sun_wl, sun_flux_norm, 4.60454, 0.001) *creates a snippet of the given spectrum (wl and flux) using a specified linelocation and line width

#"Down-the-line" average:
#1)define a single wavelength grid (pick favorite wl array and use for all spectra) 
wl_grid = pyfits.getdata(base_dir + obj_list[0])[0,:] / (beta_shifts[0] + 1) #use doppler corrected wavelegnth array from 1st object i.e. HIP 102040
#2) use new spect2 = np.interp(wl1, wl2, spec2) to bring 2nd spectrum onto wl grid #1
spect_0 = pyfits.getdata(base_dir + obj_list[0])[1,:]
spect_1 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[1])[0,:] / (beta_shifts[1] + 1), pyfits.getdata(base_dir + obj_list[1])[1,:])
spect_2 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[2])[0,:] / (beta_shifts[2] + 1), pyfits.getdata(base_dir + obj_list[2])[1,:])
spect_3 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[3])[0,:] / (beta_shifts[3] + 1), pyfits.getdata(base_dir + obj_list[3])[1,:])
spect_4 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[4])[0,:] / (beta_shifts[4] + 1), pyfits.getdata(base_dir + obj_list[4])[1,:])
spect_5 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[5])[0,:] / (beta_shifts[5] + 1), pyfits.getdata(base_dir + obj_list[5])[1,:])
#3)similarly write out the errors for each star; Do they need to be interpolated???
err_0 = pyfits.getdata(base_dir + obj_list[0])[2,:]
err_1 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[1])[0,:] / (beta_shifts[1] + 1), pyfits.getdata(base_dir + obj_list[1])[2,:])
err_2 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[2])[0,:] / (beta_shifts[2] + 1), pyfits.getdata(base_dir + obj_list[2])[2,:])
err_3 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[3])[0,:] / (beta_shifts[3] + 1), pyfits.getdata(base_dir + obj_list[3])[2,:])
err_4 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[4])[0,:] / (beta_shifts[4] + 1), pyfits.getdata(base_dir + obj_list[4])[2,:])
err_5 = np.interp(wl_grid, pyfits.getdata(base_dir + obj_list[5])[0,:] / (beta_shifts[5] + 1), pyfits.getdata(base_dir + obj_list[5])[2,:])

#4)vstack spectra and uncertainty arrays
vstacked_spect = np.vstack((spect_0, spect_1, spect_2, spect_3, spect_4, spect_5))
vstacked_err = np.stack((err_0, err_1, err_2, err_3, err_4, err_5))

#5)take the weighted average (using Ian’s an.wmean()) across the whole spectrum (pixel-by-pixel)
from analysis import wmean
weighted_avg_spect = wmean(vstacked_spect, 1 / (vstacked_err)**2, axis = 0)

plt.figure()
plt.title("Solar Twin Weighted Average Spectrum")
plt.xlabel("Wavelength (microns)")
plt.ylabel("Flux (???)")
plt.ylim(2e-14,3.5e-14)
plt.xlim(4.6,4.7)
plt.plot(wl_grid, weighted_avg_spect[0,:], color = '#fd5e53')
#plot CO rest lines (in a loop) along with spectrum from 4.6 to 4.7 microns
for i in range(wl_12_CO.shape[0] - 1):
    plt.axvline(x = wl_12_CO[i], color= 'k', linestyle=':')
for i in range(wl_13_CO.shape[0] - 1):
    plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
for i in range(wl_C_18_O.shape[0] - 1):
    plt.axvline(x = wl_C_18_O[i], color= 'k', linestyle='-')
label1 = plt.axvline(x = wl_12_CO[-1], color= 'k', linestyle=':')
label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--')
label3 = plt.axvline(x = wl_C_18_O[-1], color= 'k', linestyle='-')
plt.legend([label1,label2,label3], ["12_CO","13_CO","C_18_O"])

#define arrays for the wl and flux; each [i,:] is a snip; dimensions are [# of lines, snip_width]
snips_wl = np.ndarray(shape = (wl_13_CO.shape[0], 30))
snips_flux = np.ndarray(shape = (wl_13_CO.shape[0], 30))
#snip the lines of a single species from solar spectrum
for i in range(wl_13_CO.shape[0]):
    x , y = snip(sun_wl, sun_flux_norm, wl_13_CO[i], 0.0003)
    snips_wl[i, :] = x
    snips_flux[i, :] = y

plt.figure()
plt.title("Snipping Absorption Lines from Solar Spectrum")
plt.xlabel("Wavelength (um)")
plt.ylabel("Normalized Flux")
plt.xlim(4.6, 4.7)
plt.ylim(0.3, 0.5)
for i in range(wl_13_CO.shape[0] - 1):
    plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--', label = "13_CO")
plt.legend()
plt.plot(sun_wl, sun_flux_norm, color = rainbow[-1])
for i in range(wl_13_CO.shape[0]):
    plt.plot(snips_wl[i,:], snips_flux[i, :], color = rainbow[0])
    
#Will obviously be difficult to automate as we want the snip taken at each line center; NEXT: will want to combine fluxes using a vstack? and sum; *want on same wl_grid; can convert to velocity and thus we have a line profile

#plt.plot(snips_wl.T, snips_flux.T)
#plt.plot(snips_wl.T - wl_13_CO, snips_flux.T)
plt.figure()
plt.plot((snips_wl.T - wl_13_CO) / wl_13_CO * 3e5, snips_flux.T)
#plt.plot(np.mean((snips_wl.T - wl_13_CO) / wl_13_CO * 3e5, axis = 1), np.mean(snips_flux.T, axis = 1))
#plt.plot(np.mean((snips_wl.T - wl_13_CO) / wl_13_CO * 3e5, axis = 1), np.median(snips_flux.T, axis = 1))

#Write a function to return an updated line list, snips_wl, and snips_flux using a "skip list" of unuseable lines
def useable_lines(line_list, skip_list):
    new_list = []
    for i in range(len(line_list)):
        if i not in skip_list:
            new_list = np.append(new_list, line_list[i])
    return new_list

solar_skip_list = [0 , 1, 2, 6, 12, 14, 15, 16, 19, 21, 23, 25, 26, 27, 29, 32, 33]
useable_13_CO_lines = useable_lines(wl_13_CO, solar_skip_list)

def snip_spectrum(line_list, spectrum, wavegrid, snip_width):
    snips_wl = np.ndarray(shape = (len(line_list), indices_per_width)) # /2 for solar twins; remove for solar spectrum
    snips_flux = np.ndarray(shape = (len(line_list), indices_per_width)) # /2 for solar twins; remove for solar spectrum
    for i in range(len(line_list)):
        x , y = snip(wavegrid, spectrum, line_list[i], snip_width)
        snips_wl[i, :] = x
        snips_flux[i, :] = y
    return snips_wl.T, snips_flux.T

snips_wl, snips_flux = snip_spectrum(useable_13_CO_lines, sun_flux_norm, sun_wl, 0.0003)

plt.figure()
plt.title("Snipping Absorption Lines from Solar Spectrum")
plt.xlabel("Wavelength (um)")
plt.ylabel("Normalized Flux")
plt.xlim(4.6, 4.7)
plt.ylim(0.3, 0.5)
for i in range(wl_13_CO.shape[0] - 1):
    plt.axvline(x = wl_13_CO[i], color= 'k', linestyle='--')
label2 = plt.axvline(x = wl_13_CO[-1], color= 'k', linestyle='--', label = "13_CO")
plt.legend()
plt.plot(sun_wl, sun_flux_norm, color = rainbow[-1])
plt.plot(snips_wl, snips_flux)
    
#Will obviously be difficult to automate as we want the snip taken at each line center; NEXT: will want to combine fluxes using a vstack? and sum; *want on same wl_grid; can convert to velocity and thus we have a line profile

#plt.plot(snips_wl.T, snips_flux.T)
#plt.plot(snips_wl.T - wl_13_CO, snips_flux.T)
plt.figure()
plt.title("Useable 13CO Lines (Centered)")
plt.xlabel("Velocity (km/s)")
plt.ylabel("Normalized Flux Intensity")
plt.plot((snips_wl - useable_13_CO_lines) / useable_13_CO_lines * 3e5, snips_flux)

plt.figure()
plt.title("Stacked 13CO Lines")
plt.xlabel("Velocity (km/s)")
plt.ylabel("Normalized Flux Intensity")
plt.plot(np.mean((snips_wl - useable_13_CO_lines) / useable_13_CO_lines * 3e5, axis = 1), np.mean(snips_flux, axis = 1))


#select a star
star = 1
choose_file = int(star - 1)
obj_path = base_dir + obj_list[choose_file]
spect = pyfits.getdata(obj_path)
obj_name = pyfits.getval(obj_path, 'OBJECT')

#create variables for wavelength, flux, and uncertainty; fits files 1, 2, and 3 respectively
wl = spect[0,:]
flux = spect[1,:]
err = spect[2, :]

#normalize flux and wl
wl = wl / (beta_shifts[choose_file] + 1)
flux = flux/np.nanmedian(flux)

#determine remaining variables
useable_line_list = useable_lines(wl_13_CO, []) #choose a CO species line list
snips_wl, snips_flux = snip_spectrum(useable_line_list, flux, wl, snip_width, indices_per_width)

#create plots
plt.figure()
plt.title("Snipping Absorption Lines from " + str(obj_name))
plt.xlabel("Wavelength (um)")
plt.ylabel("Normalized Flux")
plt.xlim(4.6, 4.7)
plt.ylim(0.3, 0.5)
for i in range(line_list.shape[0] - 1):
    plt.axvline(x = line_list[i], color= 'k', linestyle='--')
label2 = plt.axvline(x = line_list[-1], color= 'k', linestyle='--', label = "13_CO")
plt.legend()
plt.plot(wl, flux)
plt.plot(snips_wl, snips_flux)

CO_species = str(CO_species)
plt.figure()
plt.title("Useable " + CO_species + " Lines (Centered)")
plt.xlabel("Velocity (km/s)")
plt.ylabel("Normalized Flux Intensity")
plt.plot((snips_wl - useable_line_list) / useable_line_list * 3e5, snips_flux)

plt.figure()
plt.title("Stacked " + CO_species + " Lines")
plt.xlabel("Velocity (km/s)")
plt.ylabel("Normalized Flux Intensity")
plt.plot(np.mean((snips_wl - useable_line_list) / useable_line_list * 3e5, axis = 1), np.mean(snips_flux, axis = 1))
    

plt.figure()


        

