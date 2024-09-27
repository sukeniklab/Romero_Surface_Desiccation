# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:34:52 2023

@author: ssukenik
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
from scipy.stats import ttest_ind, false_discovery_control, pearsonr
from matplotlib_venn import venn2, venn3, venn3_circles
import requests

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20

#%% helper functions
def peval(pval,n=1):
    if pval < 0.00001/n:
        return('****')
    elif pval < 0.0001/n:
        return('***')
    elif pval < 0.001/n:
        return('**')
    elif pval < 0.01/n:
        return('*')
    else:
        return('ns')
#%% data imports

data=pd.read_csv('../SI/Table_S1.csv')
N_proteome=6060 # proteins in S.cer proteome UP000002311

#%% TSP coverage diagram (Fig. 1A)

fig, ax = plt.subplots(figsize=[1.75,5])
neg_change = len(data[(data['S/P_mean']<1) & (data['-logPval_S/P']>1.3)])
pos_change = len(data[(data['S/P_mean']>1) & (data['-logPval_S/P']>1.3)])
no_sig_change = len(data)-neg_change-pos_change
not_detected = N_proteome-len(data)
ax.bar(0,not_detected,color='lightgrey',width=0.5,edgecolor='k',label='not detected')
ax.bar(0,no_sig_change,bottom=not_detected,color='grey',width=0.5,edgecolor='k',label='no sig change')
ax.bar(0,neg_change,bottom=not_detected+no_sig_change,color='lightgreen',width=0.5,edgecolor='k',label='enriched in pellet')
ax.bar(0,pos_change,bottom=no_sig_change+not_detected+neg_change,color='blue',width=0.5,edgecolor='k',label='enriched in sup')
ax.legend(bbox_to_anchor=(1, 0.8),fontsize=18)
ax.set_xticks([])
ax.set_ylabel('protein count', fontsize=20)
ax.spines[['right', 'top']].set_visible(False)
plt.savefig('Fig.2A.svg')

#%% reproducibility plots (Fig. S1)
fig,ax = plt.subplots(2,3,figsize=[15,10],sharex=False,sharey='row')
for i, fxn in enumerate(['T','S','P']):
    sliced = data[[fxn+'1_cf',fxn+'2_cf',fxn+'3_cf']].dropna().sort_values(by=fxn+'1_cf')
    sliced.index=np.arange(len(sliced))
    ax[1,i].plot(np.log10(sliced[fxn+'1_cf']),'.',c='black',alpha=0.8,label='rep1',lw=1,ms=1)
    ax[1,i].plot(np.log10(sliced[fxn+'2_cf']),c='blue',alpha=0.3,label='rep2')
    ax[1,i].plot(np.log10(sliced[fxn+'3_cf']),c='aqua',alpha=0.3,label='rep3')
ax[1,0].text(0,4.5,'total')
ax[1,1].text(0,4.5,'supernatant')
ax[1,2].text(0,4.5,'pellet')
ax[1,1].set_xlabel('protein index')
ax[1,0].set_ylabel('log abundance')
ax[1,0].legend(loc='lower right')

bins=np.linspace(0,6,100)
for i, fxn in enumerate(['T','S','P']):
    sliced = np.log10(data[[fxn+'1_cf',fxn+'2_cf',fxn+'3_cf']].dropna())
    ax[0,i].hist((sliced[fxn+'1_cf']),bins=bins,histtype='step',color='black',alpha=1,label='rep1',lw=2)
    ax[0,i].hist((sliced[fxn+'2_cf']),bins=bins,histtype='step',color='blue',alpha=1,label='rep2',lw=2)
    ax[0,i].hist((sliced[fxn+'3_cf']),bins=bins,histtype='step',color='aqua',alpha=1,label='rep3',lw=2)
ax[0,0].text(0,160,'total')
ax[0,1].text(0,160,'supernatant')
ax[0,2].text(0,160,'pellet')
ax[0,1].set_xlabel('log abundance')
ax[0,0].set_ylabel('protein count')
ax[0,0].legend(loc='upper right')
plt.savefig('Fig.S1.svg')

#%% correlation S/Ts (Fig. S2)

fig,ax = plt.subplots(1,3,figsize=[15,5],sharey=True,sharex=True)
ax[0].scatter(data['S1_cf']/data['T1_cf'],data['S2_cf']/data['T2_cf'],edgecolor='k',alpha=0.4,label='rep 1 vs 2')
ax[1].scatter(data['S1_cf']/data['T1_cf'],data['S3_cf']/data['T3_cf'],edgecolor='k',alpha=0.4,label='rep 1 vs 3')
ax[2].scatter(data['S2_cf']/data['T2_cf'],data['S3_cf']/data['T3_cf'],edgecolor='k',alpha=0.4,label='rep 2 vs 3')
ax[0].plot([0,2],[0,2],'--',c='k')
ax[1].plot([0,2],[0,2],'--',c='k')
ax[2].plot([0,2],[0,2],'--',c='k')
ax[0].set_xlim(0,1.6)
ax[0].set_ylim(0,1.6)
ax[0].set_xticks([0,0.5,1,1.5])
ax[0].legend()
ax[1].legend()
ax[2].legend()
ax[1].set_xlabel('S/T')
ax[0].set_ylabel('S/T')
plt.savefig('Fig.S2.svg')

#%% S/(S+P) vs S/T (Fig. 2B)

fig,ax=plt.subplots(figsize=[5,5],sharex=False,sharey=False)
sliced=data[['S/T_mean','S/SP_mean']].dropna()
ax.scatter(sliced['S/T_mean'],sliced['S/SP_mean'],s=10,edgecolor='k',lw=0.5,alpha=0.5)
r=pearsonr(sliced['S/T_mean'],sliced['S/SP_mean']).statistic
ax.plot([0,1],[0,1],'--',c='k')
ax.set_ylim(0,1)
ax.set_xlim(0,1)
ax.set_xlabel('S/T')
ax.set_ylabel('S/(S+P)')
plt.savefig('Fig.2B.svg')

#%% S/P volcano plot (Fig. 2C)

fig,ax = plt.subplots(figsize=(5,5))
ax.scatter(np.log2(data['S/P_mean']),data['-logPval_S/P'],s=10,lw=0.5,edgecolor='k',alpha=0.5)
ax.plot([0,0],[-1,8],'--',c='k')
ax.plot([-8,3],[1.3,1.3],'--',c='k')
ax.set_xlabel('$log_2(FC)$ (S/P)')
ax.set_ylabel('$-log_{10}$ (adj. Pval)')
ax.set_xlim(-7,3)
ax.set_ylim(0,5)
plt.savefig('Fig.2C.svg')

#%% S/T histogram (Fig. 2D)

quants = data['S/T_mean'].quantile(np.linspace(0,1,9)) 
fig,ax = plt.subplots(figsize=[7,7])
bgs = plt.cm.bwr_r(np.linspace(0,1,5))
bins=np.linspace(0,1,50)
ax.hist(data['S/T_mean'],bins=bins,histtype='step',fill=True,lw=3,color='purple',alpha=0.6)
ax.hist(data['S/T_mean'],bins=bins,histtype='step',lw=1,color='k')
for quant in quants.values[1:-1]:
    ax.vlines(quant,0,1000,color='grey',ls='--',zorder=0)
ax.set_xlim(0,1)
ax.set_ylabel('protein count')
ax.set_xlabel('resolubility')
ax.set_ylim(0,400)
data['S/T_mean'].median()
plt.savefig('Fig.2D.svg')