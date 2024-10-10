# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 21:20:27 2024

@author: ssukenik
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20

#%%
xl = pd.ExcelFile('../SI/Fig.3.LiP.xlsx')

#%%
categories = ['pI','CCT binding sites','SSB binding sites']
cmap = plt.cm.rainbow(np.linspace(0,1,len(categories)))
fig,ax=plt.subplots(1,len(categories),figsize=[4.5*len(categories),5],sharey=True)
for i,sheet_name in enumerate(categories):
    data = xl.parse(sheet_name,index_col=0)
    # ticks = np.a
    ax[i].bar(np.arange(len(data)),data['structurally perturbed peptide fraction'],
              color=cmap[i],alpha=0.4,edgecolor='k',lw=2)
    ax[i].set_xlabel(sheet_name)
    if sheet_name == 'Cellular location':
        ax[i].set_xticks(np.arange(len(data)))
        ax[i].set_xticklabels(data.index.astype('string'),rotation=45,ha='right',fontsize=12)
    else:
        ax[i].set_xticks(np.arange(len(data))[::2])
        ax[i].set_xticklabels(data.index.astype('string').str.split('-').str[0][::2])
ax[0].set_ylabel('fraction of \n perturbed peptides')
plt.tight_layout()
fig.savefig('Fig.S4A-C.svg')

fig,ax=plt.subplots(figsize=[10,7])
data = xl.parse('Cellular location',index_col=0)
ax.bar(np.arange(len(data)),data['structurally perturbed peptide fraction'],
          color=cmap[i],alpha=0.4,edgecolor='k',lw=2,facecolor='blue')
ax.set_xticks(np.arange(len(data)))
ax.set_xticklabels(data.index.astype('string'),rotation=45,ha='right',fontsize=20)
ax.set_ylabel('fraction of \n perturbed peptides')
plt.tight_layout()
fig.savefig('Fig.S4D.svg')
