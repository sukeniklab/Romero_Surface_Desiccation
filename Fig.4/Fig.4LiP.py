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
xl = pd.ExcelFile('Fig.4.LiP.DODO.xlsx')

#%%
fig,ax = plt.subplots(1,2,figsize=[12,2.5],sharey=True)
fig.subplots_adjust(wspace=0.15)
data = xl.parse('percentAGP',index_col=0)
ax[1].plot(np.arange(len(data)),data['structurally perturbed peptide fraction'],'--',c='purple')
ax[1].scatter(np.arange(len(data)),data['structurally perturbed peptide fraction'],s=300,edgecolor='k',c='purple')
ax[1].set_xticks(np.arange(0,len(data))[::2])
ax[1].set_xticklabels(data.index[::2],rotation=45)
ax[1].set_ylabel('')

data = xl.parse('percentDE',index_col=0)
ax[0].plot(np.arange(len(data)),data['structurally perturbed peptide fraction'],'--',c='r')
ax[0].scatter(np.arange(len(data)),data['structurally perturbed peptide fraction'],s=300,edgecolor='k',c='r')
data = xl.parse('percentKR',index_col=0)
ax[0].plot(np.arange(len(data)),data['structurally perturbed peptide fraction'],'--',c='b')
ax[0].scatter(np.arange(len(data)),data['structurally perturbed peptide fraction'],s=300,edgecolor='k',c='b')
ax[0].set_xticks(np.arange(1,len(data)))
ax[0].set_xticklabels(data.index[1:],rotation=45)
ax[0].set_ylabel('fraction of\nperturbed peptides')
fig.savefig('Fig.4F.svg')
#%%
cmap = ['green','purple','orange','red','blue','pink']
categories = xl.sheet_names
#categories = ['log(abundance)','length','Domains','% disordered','log(interactors)']
fig,ax=plt.subplots(1,len(categories),figsize=[4*len(categories),5],sharey=True)
for i,sheet_name in enumerate(categories):
    data = xl.parse(sheet_name,index_col=0)
    # ticks = np.a
    ax[i].bar(np.arange(len(data)),data['structurally perturbed peptide fraction'],
              color=cmap[i],alpha=0.4)
    ax[i].set_xlabel(sheet_name)
    ax[i].set_xticks([0,int(len(data)/2),len(data)-1])
    ax[i].text(0.1,0.9,'DODO',transform=ax[i].transAxes)
    ax[i].set_xticklabels(data.index.astype('string').str.split('-').str[0][[0,int(len(data)/2),len(data)-1]])
ax[0].set_ylabel('fraction of \n perturbed peptides')
fig.savefig('Fig.S7F.svg')
#%%
xl = pd.ExcelFile('Fig.4.LiP.chainsaw.xlsx')
#%%
cmap = ['green','purple','orange','red','blue','pink']
categories = xl.sheet_names
#categories = ['log(abundance)','length','Domains','% disordered','log(interactors)']
fig,ax=plt.subplots(1,len(categories),figsize=[4*len(categories),5],sharey=True)
for i,sheet_name in enumerate(categories):
    data = xl.parse(sheet_name,index_col=0)
    # ticks = np.a
    ax[i].bar(np.arange(len(data)),data['structurally perturbed peptide fraction'],
              color=cmap[i],alpha=0.4)
    ax[i].set_xlabel(sheet_name)
    ax[i].set_xticks([0,int(len(data)/2),len(data)-1])
    ax[i].set_xticklabels(data.index.astype('string').str.split('-').str[0][[0,int(len(data)/2),len(data)-1]])
    ax[i].text(0.1,0.9,'chainsaw',transform=ax[i].transAxes)
ax[0].set_ylabel('fraction of \n perturbed peptides')
ax[i].text(0.1,0.9,'chainsaw',transform=ax[0].transAxes)
plt.tight_layout()
fig.savefig('Fig.S7E.svg')
#%%
fig,ax = plt.subplots(1,2,figsize=[12,2.5],sharey=True)
fig.subplots_adjust(wspace=0.15)
data = xl.parse('percentAGP',index_col=0)
ax[1].plot(np.arange(len(data)),data['structurally perturbed peptide fraction'],'--',c='purple')
ax[1].scatter(np.arange(len(data)),data['structurally perturbed peptide fraction'],s=300,edgecolor='k',c='purple')
ax[1].set_xticks(np.arange(0,len(data))[::2])
ax[1].set_xticklabels(data.index[::2],rotation=45)
ax[1].text(0.05,0.05,'chainsaw',transform=ax[1].transAxes)
ax[1].set_ylabel('')

data = xl.parse('percentDE',index_col=0)
ax[0].plot(np.arange(len(data)),data['structurally perturbed peptide fraction'],'--',c='r')
ax[0].scatter(np.arange(len(data)),data['structurally perturbed peptide fraction'],s=300,edgecolor='k',c='r')
data = xl.parse('percentKR',index_col=0)
ax[0].plot(np.arange(len(data)),data['structurally perturbed peptide fraction'],'--',c='b')
ax[0].scatter(np.arange(len(data)),data['structurally perturbed peptide fraction'],s=300,edgecolor='k',c='b')
ax[0].set_xticks(np.arange(1,len(data)))
ax[0].set_xticklabels(data.index[1:],rotation=45)
ax[0].text(0.05,0.85,'chainsaw',transform=ax[0].transAxes)
ax[0].set_ylabel('fraction of\nperturbed peptides')
fig.savefig('Fig.S7G.svg')