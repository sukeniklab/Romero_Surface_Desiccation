# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:51:36 2024

@author: ssukenik
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from mpl_toolkits.axes_grid1 import Divider, Size

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20

#%%
data=pd.read_csv('../SI/Table_S6.csv')
cmaps=['Blues','Reds']

for col, subset in enumerate(['TOP','BOT']):
    for i in ['BP','MF','KEGG']:
        sliced = data[(data['subset']==subset)&(data['db']==i)]
        fig = plt.figure(figsize=[5,5])
        h = [Size.Fixed(0), Size.Fixed(2)]
        v = [Size.Fixed(0), Size.Fixed(len(sliced)*.4)]
        divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
        ax = fig.add_axes(divider.get_position(),axes_locator=divider.new_locator(nx=1, ny=1))
        ax.hlines(np.arange(len(sliced)),np.zeros(len(sliced)),np.ones(len(sliced))*4.5,
                  color='k',linestyle='--',zorder=0)
        g = sns.scatterplot(x='Fold change',y='term_name',data=sliced,hue='negative_log10_of_adjusted_p_value',
                   size='term_size',sizes=(150,800),
                   legend='brief',edgecolor='k',palette=cmaps[col],
                   hue_norm=(-7,30))
        g.legend(bbox_to_anchor=(1.2,1))
        ax.set_yticks(np.arange(len(sliced)))
        ax.set_yticklabels(sliced['term_name'],fontsize=20)
        ax.set_ylim(-0.5,len(sliced)-0.5)
        ax.set_xlim(0,4.5)
        ax.set_xticks([0,2,4])
        ax.set_xlabel('fold enrichment',fontsize=20)
        ax.set_ylabel('')
        plt.savefig(subset+'_'+i+'.svg')