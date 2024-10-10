# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:19:16 2024

@author: ssukenik
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_ind, kstest, false_discovery_control

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20
mpl.rcParams['svg.fonttype'] = 'none'

#%%

df = pd.read_csv('../SI/Table_S7.csv',index_col=0)
TSP = pd.read_csv('../SI/Table_S1.csv',index_col=0)

#%%
df = df[df['Count']>50]
colors = plt.cm.coolwarm_r(np.linspace(0,1,len(df)+1))
for i,grp in enumerate(df.index):
    df.loc[grp,'group_median'] = np.median(pd.to_numeric(df.loc[grp,'S/T_mean'][1:-2].split(',')))
df=df.sort_values('group_median')
groups=['Entire dataset']
ref=TSP['S/T_mean'].dropna().values
medianprops={'linewidth':2,'color':'k'}
boxprops={'linewidth':2,'color':'k'}
fig,ax=plt.subplots(figsize=[4,10])
box = ax.boxplot(ref,positions=[0],vert=False,showfliers=False,widths=0.5,
                 medianprops=medianprops,boxprops=boxprops,whiskerprops=boxprops,capprops=boxprops)
ax.scatter(ref,np.random.randn(len(ref))*.05,s=20,edgecolor='k',lw=0.5,alpha=0.4,color='grey')
pvals=[]
for i,grp in enumerate(df.index):
    thisGrp = df.loc[grp,'S/T_mean']
    ST = pd.to_numeric(thisGrp[1:-2].split(','))
    pvals.append(kstest(ref,ST).pvalue)
    ax.boxplot(ST,positions=[i+1],vert=False,showfliers=False,widths=0.5,
                     medianprops=medianprops,boxprops=boxprops,whiskerprops=boxprops,capprops=boxprops)
    ax.scatter(ST,i+1+np.random.randn(len(ST))*.05,s=20,edgecolor='k',
               lw=0.5,alpha=0.4,color=colors[i])
    groups.append(grp)
sig = np.where(false_discovery_control(pvals)<0.01)[0]
for i in sig:
    ax.text(0.9,i+1.45,'$*$',va='top',ha='center',fontsize=32)

ax.plot(np.ones(2)*np.median(ref),[-1,50],'--',c='r',lw=2)
ax.set_yticks(np.arange(i+2))
ax.set_yticklabels(groups)
ax.set_xlim(0,1)
ax.set_ylim(-1,i+2)
ax.set_xlim(-0.05,1.1)
ax.set_xticks([0,0.5,1])
ax.set_xlabel('resolubility')
fig.savefig('Fig.5B.svg')
