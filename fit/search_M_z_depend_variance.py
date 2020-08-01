#!/usr/bin/env python
# coding: utf-8

import math
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm
import scipy.interpolate as interpolate
import scipy.optimize as optimization
from scipy.optimize import minimize
from scipy.special import gamma, factorial, erf
from scipy import integrate
import scipy.stats as stats
import os

argvs = sys.argv
argc = len(argvs)
if (argc != 5):
    print('Usage: # python %s input output t_or_r OPT' % argvs[0])
    print('t_or_r =0: sigma_t, 1: sigma_r')
    print('OPT = 0-2')
    quit()

path = argvs[1]

ofname_txt = argvs[2] + ".txt"
ofname     = argvs[2] + ".pdf"

t_or_r = int(argvs[3])
OPT = int(argvs[4])

sns.set()
sns.set_palette('gray')
sns.set_style("ticks", {'grid.linestyle': '--'})
fig = plt.figure(figsize=(9,7))
fs=24
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

def fit_func_1D_1(x, a1, a2):
    
    return np.log(a1*np.power(x, a2))

def fit_func_1D_2(x, a1, a2):
    
    return np.log(a1*np.exp(-(a2/x)**2))

def read_data_and_fit(file):

    df = pd.read_csv(file, delim_whitespace=True, header=None)
    df.columns=['m1', 'm2', 'b1', 'b2', 'p1', 'p2', 'p3', 'q1', 'q2', 'q3']

    if t_or_r == 0:
        tab_ylabel = [r"${\rho}_{\mathrm{200t},0}$", r"${\rho}_{\mathrm{200t},1}$", 
                      r"${\rho}_{\mathrm{200t},2}$"]
    if t_or_r == 1:
        tab_ylabel = [r"${\rho}_{\mathrm{200r},0}$", r"${\rho}_{\mathrm{200r},1}$", 
                      r"${\rho}_{\mathrm{200r},2}$"]

    ax.grid(color='gray', linestyle='--', linewidth=1, alpha=0.5)
    lab_left = True
    lab_bottom = True
    plt.tick_params(labelbottom=lab_bottom, labelleft=lab_left, labelsize=fs)
    ax.set_xlabel(r'$F(R_1, R_2)$', fontsize=0.8*fs) 
    ax.set_ylabel(tab_ylabel[OPT], fontsize=fs)

    title=r"$z=$"+"${0:.2f}$".format(z)
    ax.set_title(title, fontsize=fs)

    if OPT == 0:
        if t_or_r==0:
            ax.set_xlim(20.0, 200.0)
            ax.set_ylim(20.0, 200.0)                
        if t_or_r==1:
            ax.set_xlim(50.0, 2000.0)
            ax.set_ylim(50.0, 2000.0)                

    if OPT == 1:
        if t_or_r==0:
            ax.set_xlim(0.50, 10.0)
            ax.set_ylim(0.50, 10.0)
        if t_or_r==1:
            ax.set_xlim(14.0, 28.0)
            ax.set_ylim(14.0, 28.0)
            
    if OPT == 2:
        if t_or_r==0:
            ax.set_xlim(0.5, 5.0)
            ax.set_ylim(0.5, 5.0)                
        if t_or_r==1:
            ax.set_xlim(0.3, 0.6)
            ax.set_ylim(0.3, 0.6)

    xxx = np.logspace(-5, 5, 100)
        
    rhom = 2.7754e11 * 0.315
                            
    if t_or_r == 0:
        m1 = np.array(df[(df['p1'] > -1e10)]['m1'])
        m2 = np.array(df[(df['p1'] > -1e10)]['m2'])
        b1 = np.array(df[(df['p1'] > -1e10)]['b1'])
        b2 = np.array(df[(df['p1'] > -1e10)]['b2'])
        if OPT == 0:
            data  = np.array(df[(df['p1'] > -1e10)]['p1'])
            data2 = np.array(df[(df['p1'] > -1e10)]['p3'])
        if OPT == 1:
            data  = np.array(df[(df['p1'] > -1e10)]['p2']) 
        if OPT == 2:
            data  = np.array(df[(df['p1'] > -1e10)]['p3']) 

    if t_or_r == 1:
        m1 = np.array(df[(df['q1'] > -1e10)]['m1'])
        m2 = np.array(df[(df['q1'] > -1e10)]['m2'])
        b1 = np.array(df[(df['q1'] > -1e10)]['b1'])
        b2 = np.array(df[(df['q1'] > -1e10)]['b2'])
        if OPT == 0:
            data  = np.array(df[(df['q1'] > -1e10)]['q1']) 
        if OPT == 1:
            data  = np.array(df[(df['q1'] > -1e10)]['q2']) 
        if OPT == 2:
            data  = np.array(df[(df['q1'] > -1e10)]['q3']) 
        
    r1 = np.power(3*m1/(4*np.pi*200.0*rhom), 0.3333)
    r2 = np.power(3*m2/(4*np.pi*200.0*rhom), 0.3333)
        
    xbin = np.array([])
    ybin = np.array([])
    data_bin = np.array([])
    weight_bin = np.array([])

    for i in range(len(m1)):

        if r1[i]+r2[i] < 2:

            if OPT == 0:
                xbin = np.append(xbin, r1[i]+r2[i])
                ybin = np.append(ybin, r2[i])
                weight_bin = np.append(weight_bin, 1.0)
                if t_or_r == 0:
                    zbin = data[i]/np.sqrt(r2[i])**data2[i]
                    data_bin = np.append(data_bin, np.log(zbin))
                if t_or_r == 1:
                    data_bin = np.append(data_bin, np.log(data[i]))

            if OPT == 1:            
                xbin = np.append(xbin, r1[i]+r2[i])
                ybin = np.append(ybin, r2[i])
                weight_bin = np.append(weight_bin, 1.0)
                data_bin = np.append(data_bin, np.log(data[i]))
        
            if OPT == 2:
                xbin = np.append(xbin, r1[i]+r2[i])
                ybin = np.append(ybin, r2[i])
                weight_bin = np.append(weight_bin, 1.0)
                data_bin = np.append(data_bin, np.log(abs(data[i])))
                            
    xdata = np.vstack((xbin, ybin))

    if t_or_r == 0:
        popt, pcov = optimization.curve_fit(fit_func_1D_1, xbin, data_bin, sigma=weight_bin, method='lm', maxfev=10000)
        xplot=fit_func_1D_1(xbin, popt[0], popt[1])

    if t_or_r == 1:
        if OPT == 0:
            pinit = [1000, 1.0]
            popt, pcov = optimization.curve_fit(fit_func_1D_2, xbin, data_bin, sigma=weight_bin, p0=pinit, method='lm', maxfev=10000)
            xplot=fit_func_1D_2(xbin, popt[0], popt[1])
        else:
            popt, pcov = optimization.curve_fit(fit_func_1D_1, xbin, data_bin, sigma=weight_bin, method='lm', maxfev=10000)
            xplot=fit_func_1D_1(xbin, popt[0], popt[1])

    res = 0
    for i in range(len(xplot)):
        res += (data_bin[i]-xplot[i])*(data_bin[i]-xplot[i])
        res = res / (len(xplot))
        
    print(z, res)
    
    label=r"$(a_1, a_2)=$"+"$({0:.2f}, {1:.2f})$".format(popt[0], popt[1])
        
    ax.errorbar(np.exp(xplot), np.exp(data_bin), yerr=None, xerr=None, fmt='o', mec='black', markersize=7.0, 
                mfc=CB_color_cycle[0], color=CB_color_cycle[0], alpha=0.8, label=label)
            
    ax.plot(xxx, xxx, color='black', linestyle="--", linewidth=4)

    ax.legend(loc='best', fontsize=0.4*fs)
    
    return popt[0], np.sqrt(pcov[0][0]), popt[1], np.sqrt(pcov[1][1])

with open(ofname_txt, mode='w') as f:
    file = path
    a1, err_a1, a2, err_a2 = read_data_and_fit(file)                  
    f.write("{0:.6E} {1:.6E} {2:.6E} {3:.6E} {4:.6E}\n".format(tab_z[iz], a1, err_a1, a2, err_a2))
                               
f.close()

plt.tight_layout()
fig.savefig(ofname)
