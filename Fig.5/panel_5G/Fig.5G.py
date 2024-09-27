# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:55:40 2024

@author: ssukenik
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_ind

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20
mpl.rcParams['svg.fonttype'] = 'none'
#%%

data = pd.read_csv('Fig.5G.csv')
samples = ['Total','sup']
slopes = np.zeros([len(samples),3])
errs = np.zeros([len(samples),3])
for idx,sample in enumerate(samples):
    for rep in [1,2,3]:
        x=data[(data['sample']==sample)&(data['repitition']==rep)]['mass (ug)']
        y=data[(data['sample']==sample)&(data['repitition']==rep)]['slope']
        fit,err = np.polyfit(x,y,1,cov=True)
        err = np.sqrt(np.diag(err))
        slopes[idx,rep-1]=fit[0]
        errs[idx,rep-1]=err[0]
fig,ax = plt.subplots(figsize=[5,5])
ax.bar(0,slopes.mean(axis=1)[0]/slopes.mean(axis=1)[0],edgecolor='k',lw=1,color='grey',alpha=0.85)
ax.bar(1,slopes.mean(axis=1)[1]/slopes.mean(axis=1)[0],edgecolor='k',lw=1,color='royalblue',alpha=0.85)
ax.scatter([0,0.1,-0.1],slopes[0,:]/slopes[0,1],c='k',s=20,zorder=3)
ax.scatter([1,1.1,0.9],slopes[1,:]/slopes[0,1],c='k',s=20,zorder=3)
ax.errorbar(0,1,yerr=slopes.std(axis=1)[0]/slopes.mean(axis=1)[0],capsize=5,color='k',lw=2)
ax.errorbar(1,slopes.mean(axis=1)[1]/slopes.mean(axis=1)[0],lw=2,yerr=slopes.std(axis=1)[1]/slopes.mean(axis=1)[0],capsize=5,color='k')
ax.set_xticks([0,1])
ax.set_xticklabels(['total','supernatant'])
ax.set_ylim(0,1.5)
# ax.set_yticks([0,0.5,1,1.5])
pval=ttest_ind(slopes[0,:],slopes[1,:]).pvalue
# Add annotation
ax.text(x=0.5, y=1.35, s='p-value < 0.007', ha='center', va='bottom', fontsize=20, color='black')

# Draw a line to denote the comparison
ax.plot([0, 0, 1, 1], [1.22, 1.3, 1.3, 1.22], lw=1.5, color='black')
ax.set_ylabel('CS activity / $\mu$g protein')
plt.savefig('Fig.5G.svg')
        
    