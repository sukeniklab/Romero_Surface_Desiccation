# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 12:33:19 2023

@author: poret
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns 
import numpy as np
from scipy.optimize import curve_fit

# Style parameters
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['font.size'] = 20


#%%

WLoss_all= pd.read_csv('Weight loss_all exp.csv')

t= WLoss_all['time']
w= WLoss_all['Weight  %']

#%%

def exponential_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

df_mean = WLoss_all.groupby('time').mean().reset_index()
x_data = df_mean['time'].values
y_data = df_mean['Weight  %'].values

# Initial guess for the parameters
initial_guess = [3, 1, 0]

# Perform the curve fitting
params, covariance = curve_fit(exponential_decay, x_data, y_data, p0=initial_guess)
err=np.sqrt(np.diag(covariance))
# Extract the fitted parameters
a, b, c = params
print(f"Fitted parameters: a={a}, b={b}, c={c}")

x1_data= np.arange(0, 43, 0.1)

# Generate y values using the fitted parameters
y_fitted = exponential_decay(x1_data, *params)
y_top = exponential_decay(x1_data,*(params+err))
y_bot = exponential_decay(x1_data,*(params-err))

fig, ax1 = plt.subplots(figsize=[5,5])

# Plot the data
ax1.scatter(WLoss_all['time'], WLoss_all['Weight  %'], label='Original Data', alpha=0.7, 
            s=160,linewidth=1,edgecolor='k',zorder=0)
#ax1.scatter(x_data, y_data, label='Averaged Data', color='blue')
ax1.plot(x1_data, y_fitted, label='Fitted Curve', color='darkblue')
ax1.fill_between(x1_data,y_top,y_bot,color='darkblue',alpha=0.2)
plt.xlabel('time (h)')
plt.ylabel('weight %')

means = 11.1875
y1= means-1.178747004
y2= means+1.178747004

ax1.axhline(y=means, color='r',linestyle='--')
ax1.fill_between(x1_data,y1,y2,color='r',alpha=0.2)
# ax1.axhline(y=y2, color='r', linewidth=.2)
ax1.text(0.8, 4, 'desiccated water content', fontsize=18, color='red')

ax1.fill_between(np.arange(-3, 43, 0.1), y1, y2, color='lightcoral', alpha=0.5)

ax1.set_ylim(0, 100)
ax1.set_xlim(0, 43)

#%%

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

WLoss_all = pd.read_csv('Weight loss_all exp.csv')

t = WLoss_all['time']
w = WLoss_all['Weight  %']

def exponential_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

df_mean = WLoss_all.groupby('time').mean().reset_index()
x_data = df_mean['time'].values
y_data = df_mean['Weight  %'].values

# Initial guess for the parameters
initial_guess = [3, 1, 0]

# Perform the curve fitting
params, covariance = curve_fit(exponential_decay, x_data, y_data, p0=initial_guess)

# Extract the fitted parameters
a, b, c = params
print(f"Fitted parameters: a={a}, b={b}, c={c}")

# Calculate the standard deviation of the parameters
perr = np.sqrt(np.diag(covariance))

# Generate x values for the fitted curve
x1_data = np.arange(0, 43, 0.1)

# Generate y values using the fitted parameters
y_fitted = exponential_decay(x1_data, a, b, c)

# Calculate the confidence interval
delta = perr[0] * np.exp(-b * x1_data) + perr[1] * a * x1_data * np.exp(-b * x1_data)

# Upper and lower bounds of the confidence interval
y_upper = y_fitted + delta
y_lower = y_fitted - delta

fig, ax1 = plt.subplots()

# Plot the data
ax1.scatter(WLoss_all['time'], WLoss_all['Weight  %'], label='Original Data', alpha=0.6)
ax1.plot(x1_data, y_fitted, label='Fitted Curve', color='darkblue')

# Plot the confidence interval
ax1.fill_between(x1_data, y_lower, y_upper, color='lightblue', alpha=0.5, label='Confidence Interval')

plt.xlabel('time (h)')
plt.ylabel('Weight %')

means = 11.1875
y1 = means - 1.178747004
y2 = means + 1.178747004

ax1.axhline(y=means, color='r')
ax1.axhline(y=y1, color='r', linewidth=.2)
ax1.axhline(y=y2, color='r', linewidth=.2)
ax1.text(-1, 13.2, 'Final water content 11.2%', fontsize=12, color='red')

ax1.fill_between(np.arange(-3, 43, 0.1), y1, y2, color='lightcoral', alpha=0.5)

ax1.set_xlim(-3, 43)
#ax1.legend()

plt.show()
