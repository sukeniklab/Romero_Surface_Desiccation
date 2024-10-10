# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:19:16 2024

@author: ssukenik
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_ind, kstest

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20
mpl.rcParams['svg.fonttype'] = 'none'

#%%

df = pd.read_csv('../SI/Table_S10.csv',index_col=0)
df = df[(df['N']>100)&(df['Structrually Altered']>0)]
df['effect']=np.log2(df['Structrually Altered']/df['Expected Structrually Altered'])
df['-logpval']=-np.log10(df['Chi Sq P-value'])
fig,ax=plt.subplots(figsize=[10,5])
ax.scatter(df['effect'],df['-logpval'],s=150,edgecolor='k',c='grey')
selected = df[(df['Chi Sq P-value']<0.01)&(df['effect']>1)]
ax.scatter(selected['effect'],selected['-logpval'],s=150,edgecolor='k',c='red')
selected = df[(df['Chi Sq P-value']<0.01)&(df['effect']<-1)]
ax.scatter(selected['effect'],selected['-logpval'],s=150,edgecolor='k',c='blue')
selected = df[(df['Chi Sq P-value']<0.01)&(df['effect'].abs()>1)]
for ecod in selected.index:
    ax.annotate(ecod,[selected.loc[ecod,'effect'],
                      selected.loc[ecod,'-logpval']],fontsize=16)
ax.plot([0,0],[-1,20],'--',c='k')
ax.plot([-1,-1],[-1,20],'--',c='k')
ax.plot([1,1],[-1,20],'--',c='k')
ax.plot([-4,4],-np.log10(np.ones(2)*(0.01)),'--',c='k')
ax.set_xlim(-3,3)
ax.set_ylim(0,13)
ax.set_xlabel('$log_2(FC)$ (observed/expected perturbed peptides)')
ax.set_ylabel('$-log_{10}\ \chi^2$ p-value')
fig.savefig('Fig.S11.svg')