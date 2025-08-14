# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 10:53:39 2024

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

#%% lysate BCA assays
data = pd.read_csv('Fig.1F_BCA_results.csv',index_col=0)
f_S = data['S/S+P'].mean()
f_S_err = data['S/S+P'].std()
f_P = data['P/S+P'].mean()
f_P_err = data['P/S+P'].std()
# S_P = [0.357986029, 0.642013971]
# S_P_err = [0.040775224, 0.040775224]

fig,ax=plt.subplots(figsize=[2,5])
ax.bar(0,f_P,edgecolor='k',width=0.75,color='red')
ax.bar(0,f_S,bottom=f_P,edgecolor='k',width=0.75,color='lightblue')
ax.scatter(0+np.random.randn(len(data['P/S+P']))*0.05,data['P/S+P'],c='grey',edgecolor='k')
ax.scatter(0+np.random.randn(len(data['S/S+P']))*0.05,f_P+data['S/S+P'],c='grey',edgecolor='k')
ax.errorbar(0,f_P,yerr=f_P_err,capsize=5,c='k')
ax.errorbar(0,f_S+f_P,yerr=f_S_err,capsize=5,c='k')
ax.text(0,f_P/2,'pellet',fontsize=20,ha='center')
ax.text(0,f_P+f_S/2,'sup',fontsize=20,ha='center')
ax.spines[['right', 'top']].set_visible(False)
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_ylabel('fraction of total mass')
fig.savefig('Fig.1F.svg')
# %%
