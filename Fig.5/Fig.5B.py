# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:19:16 2024

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

df = pd.read_csv('../SI/Table_S7.csv',index_col=0)
TSP = pd.read_csv('../SI/Table_S1.csv',index_col=0)

#%%

fig,ax=plt.subplots(figsize=[3,7])
groups=['Entire dataset']
ref=TSP['S/T_mean'].dropna().values
medianprops={'linewidth':2,'color':'red'}
boxprops={'linewidth':2,'color':'k'}
box = ax.boxplot(ref,positions=[0],vert=False,showfliers=False,widths=0.5,
                 medianprops=medianprops,boxprops=boxprops,whiskerprops=boxprops,capprops=boxprops)
    
ax.scatter(ref,np.random.randn(len(ref))*.05,s=20,edgecolor='k',lw=0.5,alpha=0.4)
x=1
for i in df.columns:
    N = len(df[i].dropna())
    if N > 50:
        thisGrp = df[i].dropna().values
        pval = ttest_ind(ref,thisGrp).pvalue
        if pval < 0.8:
            groups.append(str.replace(i,'\'','').split('/')[0])
            ax.boxplot(thisGrp,positions=[x],vert=False,showfliers=False,widths=0.5,
                             medianprops=medianprops,boxprops=boxprops,whiskerprops=boxprops,capprops=boxprops)
            ax.scatter(thisGrp,x+np.random.randn(len(thisGrp))*.05,s=20,edgecolor='k',lw=0.5,alpha=0.4)
            if pval < 0.01:
                ax.text(1,x,'$*$',va='center',ha='center',fontsize=32)
            x+=1

ax.plot([np.median(ref),np.median(ref)],[-1,len(groups)],'--',c='k')
ax.set_yticks(range(x))
ax.set_yticklabels(groups)
plt.set_cmap('rainbow')
ax.set_ylim(-0.5,len(groups)-0.5)
ax.set_xlim(-0.05,1.1)
ax.set_xticks([0,0.5,1])
ax.set_xlabel('resolubility')
plt.savefig('xgroups.svg')
