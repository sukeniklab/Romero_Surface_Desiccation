# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:06:42 2024

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


data = pd.read_csv('Fig.4H.csv')
fig,ax = plt.subplots(figsize=[5,5])
for charge in data['Charge'].unique():
    sliced = data[data['Charge']==charge]['ratio']
    plt.scatter(charge,sliced.mean(),s=300,linewidth=1,edgecolor='k',c='limegreen')
    plt.errorbar(charge,sliced.mean(),yerr=sliced.std(),capsize=5,lw=2,color='k')
ax.plot([0,0],[0,2],'--',c='k')
ax.set_ylim(0.3,1.1)
ax.set_xlabel('GFP charge')
ax.set_ylabel('$F_{sup}\ /\ F_{total}$')
fig.savefig('Fig.4X.svg')