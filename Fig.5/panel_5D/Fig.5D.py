# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:22:42 2024

@author: poret
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns 
from scipy import stats
import numpy as np

# Style parameters
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20
mpl.rcParams['svg.fonttype'] = 'none'
#%%

data= pd.read_csv('Fig.5D.csv')

#%%

fig,ax=plt.subplots(figsize=[5,5])
means = data.groupby(['Fraction','Date']).mean()
stds = data.groupby(['Fraction','Date']).std()
fractions = ['before','after']
for date in means.index.levels[1]:
    sliced=means.xs(date,axis=0,level=1)['Osmolarity'].sort_index(ascending=False)
    plt.plot([0,1],sliced,lw=4,zorder=0)
    plt.scatter([0,1],sliced,s=500,edgecolor='k')
for i,frac in enumerate(fractions):
    sliced=data[data['Fraction']==frac]
    mean_sliced = means.xs(frac,axis=0,level=0)
    ax.scatter(np.ones(len(sliced))*i+np.random.randn(len(sliced))*0.01,sliced['Osmolarity'],
               s=100,c='white',edgecolor='k',alpha=0.7)
ax.set_xlim(-0.5,1.5)
ax.set_xticks([0,1])
ax.set_xticklabels(['total','supernatant'])
ax.set_ylabel('$\Pi$ (mOsm)')
plt.savefig('Fig.5D.svg')