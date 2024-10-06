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
xl = pd.ExcelFile('Fig.3.LiP.2.xlsx')
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
fig.savefig('LiP_distros2.svg')
#%%
categories = ['Cellular location','pI','CCT binding sites','SSB binding sites']
cmap = plt.cm.rainbow(np.linspace(0,1,len(categories)))
fig,ax=plt.subplots(1,len(categories),figsize=[4*len(categories),5])
for i,sheet_name in enumerate(categories):
    data = xl.parse(sheet_name,index_col=0)
    # ticks = np.a
    ax[i].bar(np.arange(len(data)),data['structurally perturbed peptide fraction'],
              color=cmap[i],alpha=0.4)
    ax[i].set_xlabel(sheet_name)
    if sheet_name == 'Cellular location':
        ax[i].set_xticks(np.arange(len(data)))
        ax[i].set_xticklabels(data.index.astype('string'),rotation=45,ha='right',fontsize=12)
    else:
        ax[i].set_xticks(np.arange(len(data))[::2])
        ax[i].set_xticklabels(data.index.astype('string').str.split('-').str[0][::2])
ax[0].set_ylabel('fraction of \n perturbed peptides')
plt.tight_layout()

#%%
N_bins=8
cmap = plt.cm.Purples(np.linspace(0.2,1,N_bins))
data = pd.read_csv('Fig.3.LiP.csv',index_col=0)
categories = ['Median Abundance','Isoelectric Point','Percent Disordered','Length','CCT sites','Log(Interactors)','pSup_42C',]
fig,ax = plt.subplots(1,len(categories),figsize=[5*len(categories),5],sharey=True)
for i,cat in enumerate(categories):
    sliced=data.dropna(subset=cat)
    sliced.loc[:,cat] = pd.to_numeric(sliced[cat], errors='coerce').dropna()
    bins = sliced[cat].quantile(np.linspace(1/N_bins,1,N_bins)).values
    total = sliced.loc[sliced[cat]<bins[0],'No. of Valid Peptides'].sum()
    modified = sliced.loc[sliced[cat]<bins[0],'No. of Significant Peptides (Adj. P-value)'].sum()
    ax[i].bar(0,modified/total,color=cmap[i],edgecolor='k',lw=1)
    for i_bin in np.arange(1,len(bins)):
        ax[i].plot([i_bin-0.5,i_bin-0.5],[0,1],'--',color='grey')
        total = sliced.loc[(sliced[cat]>bins[i_bin-1])&(sliced[cat]<bins[i_bin]),'No. of Valid Peptides'].sum()
        modified = sliced.loc[(sliced[cat]>bins[i_bin-1])&(sliced[cat]<bins[i_bin]),'No. of Significant Peptides (Adj. P-value)'].sum()
        ax[i].text((i_bin)/(N_bins),1,str("{0:.2g}".format(bins[i_bin-1])),rotation=90,ha='right',va='top',fontsize=16,transform = ax[i].transAxes)
        ax[i].bar(i_bin,modified/total,color=cmap[i],edgecolor='k',lw=1)
    ax[i].set_xticks(np.arange(i_bin+1))
    ax[i].set_xticklabels(np.arange(i_bin+1)+1)
    ax[i].set_xlabel(cat)
    ax[i].set_xlim(-0.5,N_bins-0.5)
    ax[i].set_ylim(0,0.40)
ax[0].set_ylabel('fraction of sig. peptides')
#%%

cmap = ['lightgreen','purple','orange','red','blue']
data = pd.read_csv('Fig.3.LiP.csv',index_col=0)
categories = ['Median Abundance','Length','Number of Domains','Percent Disordered','Interactors']
minmax = [[0,9000,8],[1,2000,15],[-1,5,7],[0,100,15],[0,150,9]]
fig,ax = plt.subplots(1,len(categories),figsize=[3.5*len(categories),4],sharey=False)
fig.subplots_adjust(hspace=0.2, wspace=0.7)
for i,cat in enumerate(categories):
    print(cat)
    sliced=data.dropna(subset=cat)
    sliced.loc[:,cat] = pd.to_numeric(sliced[cat], errors='coerce').dropna()
    bins = np.linspace(minmax[i][0],minmax[i][1],minmax[i][2]).astype('float')
    width = 0.8*(bins[1]-bins[0])
    total = sliced.loc[sliced[cat]<=bins[0],'No. of Valid Peptides'].sum()
    modified = sliced.loc[sliced[cat]<=bins[0],'No. of Significant Peptides (Adj. P-value)'].sum()
    print(bins[0])
    print(total)
    ax[i].bar(bins[0],modified/total,color=cmap[i],edgecolor='k',lw=1,width=width,alpha=0.4)
    for i_bin in np.arange(1,len(bins)):
        total = sliced.loc[(sliced[cat]>bins[i_bin-1])&(sliced[cat]<=bins[i_bin]),'No. of Valid Peptides'].sum()
        modified = sliced.loc[(sliced[cat]>bins[i_bin-1])&(sliced[cat]<=bins[i_bin]),'No. of Significant Peptides (Adj. P-value)'].sum()
        print(bins[i_bin])
        print(total)
        # ax[i].text(1/(N_bins-1)*(i_bin-0.5),1,str(int(total)),rotation=90,ha='center',va='top',fontsize=16,transform = ax[i].transAxes)
        ax[i].bar(bins[i_bin],modified/total,color=cmap[i],edgecolor='k',lw=1,width=width,alpha=0.4)
        # ax[i].plot(np.ones(2)*bins[i_bin]+width*.65,[0,1],'--',color='grey')
    # ax[i].set_xticks(np.arange(i_bin+1))
    # ax[i].set_xticklabels(np.arange(i_bin+1)+1)
    ax[i].set_xlabel(cat)
    ax[i].set_xlim(bins.min()+width*.65,bins.max()+width*.65)
    # ax[i].set_ylim(0,0.5)
ax[0].set_ylabel('fraction of \n perturbed peptides')
# plt.tight_layout()
