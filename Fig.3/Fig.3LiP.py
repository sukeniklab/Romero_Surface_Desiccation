# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 09:23:19 2024

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
mpl.rcParams['svg.fonttype'] = 'none'

#%%
xl = pd.ExcelFile('./Fig.3.LiP.xlsx')
cmap = ['green','purple','orange','red','blue']
categories = ['log(abundance)','length','Domains','% disordered','log(interactors)']
fig,ax=plt.subplots(1,len(categories),figsize=[4*len(categories),5],sharey=True)
for i,sheet_name in enumerate(categories):
    data = xl.parse(sheet_name,index_col=0)
    # ticks = np.a
    ax[i].bar(np.arange(len(data)),data['structurally perturbed peptide fraction'],
              color=cmap[i],alpha=0.4)
    ax[i].set_xlabel(sheet_name)
    ax[i].set_xticks(np.arange(len(data))[::4])
    ax[i].set_xticklabels(data.index.astype('string').str.split('-').str[0][::4])
ax[0].set_ylabel('fraction of \n perturbed peptides')
plt.tight_layout()
fig.savefig('Fig.3G-K.svg')

