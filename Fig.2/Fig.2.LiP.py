# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:40:52 2024

@author: ssukenik
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_ind, pearsonr, kstest

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20

#%% Peptide level LiP Volcano plot
data = pd.read_csv('../SI/Table_S2.csv')

fig,ax = plt.subplots(figsize=[5,5])
HT = data[data['Half Tryptic']==True]
HT['color']='lightcoral'
FT = data[data['Half Tryptic']==False]
FT['color']='slategrey'
df = pd.concat([HT,FT],axis=0).sample(frac=1)
ax.scatter(df['Log2 FC'],df['-Log10 Adj. P-value'],s=3,alpha=0.4,c=df['color'])
# ax.scatter(HT['Log2 FC'],HT['-Log10 Adj. P-value'],s=3,alpha=0.4,c='k',label='half-tryptic')
# ax.scatter(FT['Log2 FC'],FT['-Log10 Adj. P-value'],s=3,alpha=0.4,c='violet',label='full-tryptic',
#            zorder=np.random.choice([0,1],len(FT)))
ax.vlines(0,-5,7,color='k',linestyle='--')
ax.hlines(1.3,-30,30,color='k',linestyle='--')
ax.set_ylim(-0.2,6)
ax.set_xlim(-18,18)
ax.set_xlabel('normalized $log_2(FC)$ (S/T)')
ax.set_ylabel('-log$_{10}$ adj. p-value')
plt.savefig('Fig.2F.svg')
plt.savefig('Fig.2F.png')
#%%

TSP = pd.read_csv('../SI/Table_S1.csv',index_col=0)
LiP = pd.read_csv('../SI/Table_S3.csv',index_col=0)
quants = TSP['S/T_mean'].quantile(np.linspace(0,1,9)) 

#%% Protein level LiP retained/perturbed histogram

bins = np.linspace(0,1,50)
retained = LiP[LiP['Retained']==True]
perturbed = LiP[LiP['Perturbed']==True]
ttest_ind(retained['S/T Avg'],perturbed['S/T Avg'])
fig,ax = plt.subplots(figsize=[5,5])
ax.hist(retained['S/T Avg'],bins=bins,histtype='step',density=False,lw=3,color='b',label='retained')
ax.hist(perturbed['S/T Avg'],bins=bins,histtype='step',density=False,lw=3,color='r',label='perturbed')
ax.legend(title='strucutre:',fontsize=16)
ax.vlines(retained['S/T Avg'].median(),0,100,color='blue',ls='--',zorder=0,alpha=0.5,lw=3)
ax.vlines(perturbed['S/T Avg'].median(),0,100,color='red',ls='--',zorder=0,alpha=0.5,lw=3)
# ax.set_ylim(0,6.5)
ax.set_xlabel('resolubility')
ax.set_ylim(0,90)
ax.set_ylabel('protein count')
pval = kstest(retained['S/T Avg'],perturbed['S/T Avg']).pvalue
print('pval = %s\nretained median = %s\nperturbed median = %s' % (pval,
                                                                  retained['S/T Avg'].median(),perturbed['S/T Avg'].median()))

plt.savefig('Fig.2G.svg')
#%% Protein level LiP coverage

fix,ax = plt.subplots(figsize=[3,5])
ax.bar(0,len(LiP)-len(retained)-len(perturbed),color='white',edgecolor='k',linewidth=2,label = 'no significant change')
ax.bar(0,len(retained),bottom=len(LiP)-len(retained)-len(perturbed),color='b',edgecolor='k',linewidth=2, label = 'structure retained')
ax.bar(0,len(perturbed),bottom=len(LiP)-len(perturbed),color='red',edgecolor='k',linewidth=2, label = 'structure altered')
ax.bar(0,len(TSP) - len(LiP),bottom=len(LiP),color='darkgrey',edgecolor='k',linewidth=2, label = 'not detected')
ax.legend(bbox_to_anchor=(1, 0.8),fontsize=12)
ax.set_xticks([])
ax.set_ylabel('protein count')
#%% Protein level LiP correlation with resolubility

fig,ax=plt.subplots(figsize=[5,5])
plt.subplots_adjust(wspace=0.3)
resol=[]
for qidx in range(len(quants)-1):
    print(quants.iloc[qidx],quants.iloc[qidx+1])
    sliced = LiP[(LiP['S/T Avg'] < quants.iloc[qidx+1])&(LiP['S/T Avg'] > quants.iloc[qidx])]
    refoldable = len(sliced[sliced['Retained']==True])
    nonrefoldable = len(sliced[sliced['Perturbed']==True])
    ax.bar(qidx + 0.5,refoldable/(nonrefoldable+refoldable),width=0.8,edgecolor='k',color='b',label='retained')
    ax.bar(qidx + 0.5,1-refoldable/(nonrefoldable+refoldable),bottom=refoldable/(nonrefoldable+refoldable),width=0.8,color='red',label='perturbed')
    ax.bar(qidx + 0.5,refoldable/(nonrefoldable+refoldable),width=0.8,edgecolor='k',color='None')
    ax.bar(qidx + 0.5,1-refoldable/(nonrefoldable+refoldable),bottom=refoldable/(nonrefoldable+refoldable),width=0.8,edgecolor='k',color='None')
    ax.text(qidx + 0.5, 1.02, str(refoldable+nonrefoldable),fontsize=14,ha='center')
    resol.append(refoldable/(nonrefoldable+refoldable))
    # ax.bar(quant - 0.1,(norm-refoldable-nonrefoldable)/norm,bottom=(nonrefoldable+refoldable)/norm,edgecolor='k',width=0.19,color='grey',label='structure altered')
ax.set_xlabel('resolubility quantile')
ax.set_xticks(np.linspace(0,8,5))
ax.set_xticklabels([0,0.25,0.5,0.75,1])
ax.set_ylabel('structural retention')
ax.set_ylim(0,1.1)
pr = pearsonr(np.arange(len(quants)-1),np.array(resol))
plt.savefig('Fig.2H.svg')
