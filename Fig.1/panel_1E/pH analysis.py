# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:26:12 2024

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

pH_all= pd.read_csv('pH.csv')

#%% pH desiccation event 

sns.barplot(data=pH_all, x='Fraction', y='pH', hue='Date')
plt.legend(loc='upper right',bbox_to_anchor=(1.55, 1))

#%% aggregated pH measurments per experiment date

# Values to aggregate by
dates=pH_all['Date'].unique()
fractions= pH_all['Fraction'].unique()

# Define a data frame for aggregated data
pH_aggregated_df = pd.DataFrame(columns=['Date','Fraction','pH_mean','pH_err'])

# Aggregate data based on criteria and populate new dataframe
for date in dates:
    for fraction in fractions:
        pH_mean = pH_all[(pH_all.Date==date) & 
                              (pH_all.Fraction==fraction)].pH.mean()
        pH_err = pH_all[(pH_all.Date==date) & 
                             (pH_all.Fraction==fraction)].pH.std()
        pH_aggregated_df = pH_aggregated_df.append({'Date':date,'Fraction':fraction,
                                                            'pH_mean':pH_mean,
                                                            'pH_err':pH_err}, ignore_index=True)

#%% bar plot

pH_tot= pd.to_numeric(pH_aggregated_df.pH_mean[pH_aggregated_df.Fraction=='total'])
pH_sup= pd.to_numeric(pH_aggregated_df.pH_mean[pH_aggregated_df.Fraction=='supernatant'])

# Perform a t-test
t_stat, p_value = stats.ttest_ind(pH_tot, pH_sup)

bar_colors = ['gray', 'blue']
point_colors= ['dimgray', 'navy']

fig, ax = plt.subplots(figsize=[5,5])
 
sns.barplot(data=pH_aggregated_df, x='Fraction', y='pH_mean', capsize=0.1, palette=bar_colors, alpha= 0.85)
sns.swarmplot(data=pH_aggregated_df, x="Fraction", y="pH_mean", hue="Date", dodge=False, color= 'black',size=10,alpha=0.5)

pH_means= pH_aggregated_df.pH_mean.mean()
pH_errors=pH_aggregated_df.pH_mean.sem()

# Define the position for the p-value annotation
pH_y_max = pH_means + pH_errors + 1.5  # position above the highest bar
y_margin = 1  # Additional margin above the highest bar

# Set y-axis limit
ax.set_ylim(0, pH_y_max + y_margin)

# Add annotation
ax.text(x=0.5, y=pH_y_max +0.05, s='n.s.', ha='center', va='bottom')

# Optionally draw a line between bars to denote the comparison
ax.plot([0, 0, 1, 1], [pH_y_max - 0.5, pH_y_max, pH_y_max, pH_y_max - 0.5], lw=1.5, color='black')
ax.legend_.remove()
ax.set_xlabel('')
ax.set_ylabel('pH')
