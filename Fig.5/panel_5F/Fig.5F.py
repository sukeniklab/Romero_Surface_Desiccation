# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:58:31 2024

@author: poret
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns 
from scipy import stats

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20
mpl.rcParams['svg.fonttype'] = 'none'
#%%

ATPase= pd.read_csv('Fig.5F.csv')


#%%

tot = pd.to_numeric(ATPase['Relative ATPase activity'][ATPase['Fraction'] == 'Total']).dropna()
sup = pd.to_numeric(ATPase['Relative ATPase activity'][ATPase['Fraction'] == 'Sup']).dropna()

# Perform a t-test
t_stat, p_value = stats.ttest_ind(tot, sup)

bar_colors = ['gray','royalblue']

fig, ax = plt.subplots(figsize=[5,5])

# Create bar plot
sns.barplot(data=ATPase, x='Fraction', y='Relative ATPase activity', edgecolor='k',palette=bar_colors, capsize=0.1, alpha=0.85)
sns.swarmplot(data=ATPase, x="Fraction", y="Relative ATPase activity", dodge=False, color='black')

# Calculate means and errors
means = ATPase.groupby('Fraction')['Relative ATPase activity'].mean()
errors = ATPase.groupby('Fraction')['Relative ATPase activity'].std()

# Define the position for the p-value annotation
y_max = ATPase['Relative ATPase activity'].max() + errors.max()+.05 # Position above the highest bar
y_margin = .2  # Additional margin above the highest bar

# Set the y-axis limit to ensure space for the annotation
ax.set_ylim(bottom=0, top=y_max + y_margin)

# Add annotation
ax.text(x=0.5, y=y_max, s='p-value < $10^{-7}$', ha='center', va='bottom', fontsize=20, color='black')

# Draw a line to denote the comparison
ax.plot([0, 0, 1, 1], [y_max - 0.05, y_max, y_max, y_max - 0.05], lw=1.5, color='black')

# Customize the plot
ax.set_xticks([0,1])
ax.set_xticklabels(['total','supernatant'])
ax.set_xlabel('')
ax.set_ylabel('ATPase activity / μg protein')
plt.savefig('Fig.5F.svg')

