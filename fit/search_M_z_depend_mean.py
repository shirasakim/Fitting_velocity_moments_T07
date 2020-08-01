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
if (argc != 4):
    print('Usage: # python %s input output OPT' % argvs[0])
    quit()

path = argvs[1]

ofname_txt = argvs[2] + ".txt"
ofname     = argvs[2] + ".pdf"

OPT = int(argvs[3])

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

def fit_func_1D_2(x, a1, a2, a3):
    
    sx = a2*x

    return a1 * np.power(sx, a3) / (1.0 + np.power(sx, a3))

def read_data_and_fit(file):

    df = pd.read_csv(file, delim_whitespace=True, header=None)
    df.columns=['m1', 'm2', 'b1', 'b2', 'p1', 'p2', 'p3']

    ax = fig.add_subplot(1,1,1)

    tab_ylabel = [r"${\rho}_{\mathrm{PDF},0}$", r"${\rho}_{\mathrm{PDF},1}$", r"${\rho}_{\mathrm{PDF},2}$"]

    ax.grid(color='gray', linestyle='--', linewidth=1, alpha=0.5)
    lab_left = True
    lab_bottom = True
    plt.tick_params(labelbottom=lab_bottom, labelleft=lab_left, labelsize=fs)
    if OPT <= 1: ax.set_xlabel(r'$a_1 \, (M_{1,13}+M_{2,13})^{a_2}$', fontsize=0.8*fs) 
    if OPT == 2: ax.set_xlabel(r'$a_1 \, ((b_1+b_2)/a_2)^{a_3} \, [1+((b_1+b_2)/a_2)^{a_3}]^{-1}$', fontsize=0.8*fs)
    ax.set_yscale("log")
    ax.set_ylabel(tab_ylabel[OPT], fontsize=fs)
    
    ax.set_xscale("log")

    if OPT == 0:
        ax.set_xlim(10.0, 300.0)
        ax.set_ylim(10.0, 300.0)                
    if OPT == 1:
        ax.set_xlim(1.0, 2.5)
        ax.set_ylim(1.0, 2.5)                
    if OPT == 2:
        ax.set_xlim(0.7, 1.4)
        ax.set_ylim(0.7, 1.4)                

    xxx = np.logspace(-3, 3, 100)
        
    rhom = 2.7754e11 * 0.315
                                
    m1 = np.array(df[(df['p1'] > -1e10) & (df['p1'] < 300) & (df['p2'] > -2.5)]['m1'])
    m2 = np.array(df[(df['p1'] > -1e10) & (df['p1'] < 300) & (df['p2'] > -2.5)]['m2'])
    b1 = np.array(df[(df['p1'] > -1e10) & (df['p1'] < 300) & (df['p2'] > -2.5)]['b1'])
    b2 = np.array(df[(df['p1'] > -1e10) & (df['p1'] < 300) & (df['p2'] > -2.5)]['b2'])

    if OPT == 0:
        data  = np.array(df[(df['p1'] > -1e10) & (df['p1'] < 300) & (df['p2'] > -2.5)]['p1']) 
    if OPT == 1:
        data  = np.array(df[(df['p1'] > -1e10) & (df['p1'] < 300) & (df['p2'] > -2.5)]['p2']) 
    if OPT == 2:
        data  = np.array(df[(df['p1'] > -1e10) & (df['p1'] < 300) & (df['p2'] > -2.5)]['p3']) 
        
    r1 = np.power(3*m1/(4*np.pi*200.0*rhom), 0.3333)
    r2 = np.power(3*m2/(4*np.pi*200.0*rhom), 0.3333)
        
    xbin = np.array([])
    ybin = np.array([])
    data_bin = np.array([])
    weight_bin = np.array([])

    for i in range(len(m1)):

        if OPT == 0:
            judge = 0
            if (m1[i]+m2[i])/1e13 > 50 and data[i] < 100: judge = 1
            if judge == 0:
                xbin = np.append(xbin, m1[i]/1e13+m2[i]/1e13)
                data_bin = np.append(data_bin, np.log(np.abs(data[i])))
                weight_bin = np.append(weight_bin, 1.0)

        if OPT == 1:
            judge = 0
            if (m1[i]+m2[i])/1e13 > 10 and abs(data[i]) < 1.2: judge = 1
            if judge == 0:
                xbin = np.append(xbin, m1[i]/1e13+m2[i]/1e13)
                data_bin = np.append(data_bin, np.log(np.abs(data[i])))
                weight_bin = np.append(weight_bin, 1.0)

        if OPT == 2:
            if (b1[i]+b2[i]) < 5:
                xbin = np.append(xbin, (b1[i]+b2[i]))
                data_bin = np.append(data_bin, data[i])
                weight_bin = np.append(weight_bin, 1.0)
                
    if OPT == 0:
        pinit = [30.0, 0.5]
        popt, pcov = optimization.curve_fit(fit_func_1D_1, xbin, data_bin, sigma=weight_bin, method='lm')
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

        return popt[0], np.sqrt(pcov[0][0]), popt[1], np.sqrt(pcov[1][1]) #, popt[2], np.sqrt(pcov[2][2])

    if OPT == 1:
        pinit = [0.70, 0.25]
        popt, pcov = optimization.curve_fit(fit_func_1D_1, xbin, data_bin, sigma=weight_bin, method='lm')

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

    if OPT == 2:
        pinit = [1.0, 1.0, 3.0]
        popt, pcov = optimization.curve_fit(fit_func_1D_2, xbin, data_bin, sigma=weight_bin, method='lm')

        xplot=fit_func_1D_2(xbin, popt[0], popt[1], popt[2])

        res = 0
        for i in range(len(xplot)):
            res += (data_bin[i]-xplot[i])*(data_bin[i]-xplot[i])
            res = res / (len(xplot))
        
        print(z, res)
    
        label=r"$(a_1, a_2, a_3)=$"+"$({0:.2f}, {1:.2f}, {2:.2f})$".format(popt[0], popt[1], popt[2])
        
        ax.errorbar((xplot), (data_bin), yerr=None, xerr=None, fmt='o', mec='black', markersize=7.0, 
                    mfc=CB_color_cycle[0], color=CB_color_cycle[0], alpha=0.8, label=label)
            
        ax.plot(xxx, xxx, color='black', linestyle="--", linewidth=4)

        ax.legend(loc='best', fontsize=0.4*fs)

        return popt[0], np.sqrt(pcov[0][0]), popt[1], np.sqrt(pcov[1][1]), popt[2], np.sqrt(pcov[2][2])

with open(ofname_txt, mode='w') as f:

        file = path

        if OPT == 0:
            a1, err_a1, a2, err_a2 = read_data_and_fit(file)                  
            f.write("{0:.6E} {1:.6E} {2:.6E} {3:.6E} {4:.6E}\n".format(
                 tab_z[iz], a1, err_a1, a2, err_a2))
        if OPT == 1:
            a1, err_a1, a2, err_a2 = read_data_and_fit(file)                  
            f.write("{0:.6E} {1:.6E} {2:.6E} {3:.6E} {4:.6E}\n".format(
                 tab_z[iz], a1, err_a1, a2, err_a2))
        if OPT == 2:
            a1, err_a1, a2, err_a2, a3, err_a3 = read_data_and_fit(file)                  
            f.write("{0:.6E} {1:.6E} {2:.6E} {3:.6E} {4:.6E} {5:.6E} {6:.6E}\n".format(
                 tab_z[iz], a1, err_a1, a2, err_a2, a3, err_a3))
                               
f.close()

plt.tight_layout()
fig.savefig(ofname)
