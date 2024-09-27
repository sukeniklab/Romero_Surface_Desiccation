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

# Style parameters
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20


#%%

yeast_via= pd.read_csv('Desiccated Yeast Viability .csv')

#%% osmolarity desiccation event 

sns.barplot(data=yeast_via, x='Growth phase', y='Viability')
plt.yscale('log')

#%%

y_log = yeast_via[yeast_via['Growth phase'] == 'logarithmic']['Viability'].dropna()
y_stat = yeast_via[yeast_via['Growth phase'] == 'stationary']['Viability'].dropna()
# Perform a t-test
t_stat, p_value = stats.ttest_ind(y_log, y_stat)
print(p_value)

bar_colors = ['gold','lightskyblue']

fig, ax = plt.subplots(figsize=[5,5])

# Create bar plot
sns.barplot(data=yeast_via, x='Growth phase', y='Viability', 
            linewidth=2,edgecolor='k',capsize=0.1, palette=bar_colors, alpha=0.85)
# Create swarm plot
sns.scatterplot(data=yeast_via, x="Growth phase", y="Viability", color='black')

# Calculate means and errors
means = yeast_via.groupby('Growth phase')['Viability'].mean()
errors = yeast_via.groupby('Growth phase')['Viability'].sem()

# Define the position for the p-value annotation
y_max = yeast_via['Viability'].max() + errors.max() + 0.75 # Position above the highest bar
y_margin = 3.5  # Additional margin above the highest bar

# Set the y-axis limit to ensure space for the annotation
ax.set_ylim(bottom=yeast_via['Viability'].min(), top=y_max + y_margin)

# Add annotation
ax.text(x=0.5, y=y_max+0.5, s='p-value < 0.002', ha='center', va='bottom', fontsize=16, color='black')

# Draw a line to denote the comparison
ax.plot([0, 0, 1, 1], [y_max - 0.5, y_max, y_max, y_max - 0.5], lw=1.5, color='black')

# Customize the plot
ax.set_xlabel('')
ax.set_ylabel('$CFU_{after}\ /\ CFU_{before}$')
ax.set_yscale('log')


