# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 14:10:45 2021

@author: David
"""

#Calculate 4 chi-sqared values for each star; Plot those points and the parabola that fits the data; determine the minimum and dispaly as xSolar abundance and as a dex value

def calc_abundance(star_num, skip_list, snip_width, abundance_guess, plot):
    ###requires stack_data() and chi_sqr() functions
    #requires scipy.optimize
    
    #star_num = a number 1-6; used to pick which star we're looking at
    #skip_list = list of lines to skip in stacked absorption line analysis
    #how wide to make the snips; default = 0.0006
    #abundance_guess = number provided for finding the chi-sqr minimum i.e. your best guess at what the xSolar abundance is
    #if plot = true, show the plot of chi-sqr values, data fit, and min
    
    #load chosen star and ALL FOUR SOLAR MODELS; prep them for stack_data function
    
    
    #run chosen star through stack_data to obtain: stack_vel, stack_flux, and stack_err
    
    
    #run all four solar models through stack_data USING GIVEN SKIPLIST to produce: stack_vel_mod, stack_flux_mod
    ###must use skip_list corresponding to given star_num
    
    
    #calculate 4 chi_sqr values and put them in an array
    chi_sqr_list = [chi_sqr(stack_vel, stack_vel_mod1, stack_flux, stack_flux_mod1, stack_err), chi_sqr(stack_vel, stack_vel_mod2, stack_flux, stack_flux_mod2, stack_err), chi_sqr(stack_vel, stack_vel_mod3, stack_flux, stack_flux_mod3, stack_err), chi_sqr(stack_vel, stack_vel_mod4, stack_flux, stack_flux_mod4, stack_err)]
    
    #create a matching "x-array" of xSolar abundances i.e. [0, 1/3, 1, 3]
    x_list = [0.0, 1/3, 1.0, 3.0]
    
    #fit them with a degree 2 polynomial
    p_coeff = np.polyfit(x_list, chi_sqr_list, deg = 2)
    parabola = np.poly1d(p_coeff)
    
    #find the minimum of the parabola
    fit_min = optimize.fmin(parabola, abundance_guess)
    fit_min_dex = np.log10(fit_min)
    
    #plot if option plot= True
    if (plot == True):
        plt.figure()
        plt.title("$\chi^2$ Fit")
        plt.xtitle("xSolar Abundance")
        plt.ytitle("$\chi^2$")
        plt.scatter(x_list, chi_sqr_list, marker = 'o', label = "$\chi^2$ Values")
        plt.plot(np.linspace(0, 3, 1000), parabola(np.linspace(0, 3, 1000)), label = "Data Fit")
        plt.plot(np.linspace(0, 3, 1000), fit_min, marker = "X") #different color than previous
        
    else:
        pass
    
    return fit_min, fit_min_dex
    
    