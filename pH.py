#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters, Parameter, report_fit,  fit_report
from scipy.integrate import odeint
from scipy import interpolate
import scipy
import time
import pickle
from scipy.signal import savgol_filter


###Data tigr15_1_SLMB
# OTR = pd.read_excel("OTR_result.xlsx",sheet_name='tigr15_1_SLMB')
# OTR.keys()
# OTR_raw = OTR['OTR_raw']

  
# ls = np.transpose([OTR['OTR_raw']]).flatten()
# OTR_sg = savgol_filter(ls, window_length = 144, polyorder = 2, deriv=0)


# fig, ax = plt.subplots()
# ax.plot(OTR['Day'], OTR_raw, label = 'Raw')
# ax.plot(OTR['Day'], OTR_sg, label = 'SG')
# ax.set_title('OTR',fontsize=15)
# ax.set_xlabel('Time (Day)',fontsize=15)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)    
# leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
# plt.show()

#%%
###Data tigr15_1
pH_tigr15_1 = pd.read_excel("pH_result.xlsx",sheet_name='tigr15_1')
pH_tigr15_1.keys()
pH_tigr15_1_raw = pH_tigr15_1['pH_raw']


ls = np.transpose([pH_tigr15_1['pH_raw']]).flatten()
pH_tigr15_1_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr15_1['Day'], pH_tigr15_1_raw, label = 'Raw')
ax.plot(pH_tigr15_1['Day'], pH_tigr15_1_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr15_1',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr15_2
pH_tigr15_2 = pd.read_excel("pH_result.xlsx",sheet_name='tigr15_2')
pH_tigr15_2.keys()
pH_tigr15_2_raw = pH_tigr15_2['pH_raw']


ls = np.transpose([pH_tigr15_2['pH_raw']]).flatten()
pH_tigr15_2_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr15_2['Day'], pH_tigr15_2_raw, label = 'Raw')
ax.plot(pH_tigr15_2['Day'], pH_tigr15_2_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr15_2',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr15_3
pH_tigr15_3 = pd.read_excel("pH_result.xlsx",sheet_name='tigr15_3')
pH_tigr15_3.keys()
pH_tigr15_3_raw = pH_tigr15_3['pH_raw']

ls = np.transpose([pH_tigr15_3['pH_raw']]).flatten()
pH_tigr15_3_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr15_3['Day'], pH_tigr15_3_raw, label = 'Raw')
ax.plot(pH_tigr15_3['Day'], pH_tigr15_3_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr15_3',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr20_6
pH_tigr20_6 = pd.read_excel("pH_result.xlsx",sheet_name='tigr20_6')
pH_tigr20_6.keys()
pH_tigr20_6_raw = pH_tigr20_6['pH_raw']


ls = np.transpose([pH_tigr20_6['pH_raw']]).flatten()
pH_tigr20_6_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr20_6['Day'], pH_tigr20_6_raw, label = 'Raw')
ax.plot(pH_tigr20_6['Day'], pH_tigr20_6_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr20_6',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr20_7
pH_tigr20_7 = pd.read_excel("pH_result.xlsx",sheet_name='tigr20_7')
pH_tigr20_7.keys()
pH_tigr20_7_raw = pH_tigr20_7['pH_raw']


ls = np.transpose([pH_tigr20_7['pH_raw']]).flatten()
pH_tigr20_7_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr20_7['Day'], pH_tigr20_7_raw, label = 'Raw')
ax.plot(pH_tigr20_7['Day'], pH_tigr20_7_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr20_7',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr20_8
pH_tigr20_8 = pd.read_excel("pH_result.xlsx",sheet_name='tigr20_8')
pH_tigr20_8.keys()
pH_tigr20_8_raw = pH_tigr20_8['pH_raw']

  
ls = np.transpose([pH_tigr20_8['pH_raw']]).flatten()
pH_tigr20_8_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr20_8['Day'], pH_tigr20_8_raw, label = 'Raw')
ax.plot(pH_tigr20_8['Day'], pH_tigr20_8_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr20_8',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()


# %%
###Data tigr21_10
pH_tigr21_10 = pd.read_excel("pH_result.xlsx",sheet_name='tigr21_10')
pH_tigr21_10.keys()
pH_tigr21_10_raw = pH_tigr21_10['pH_raw']


ls = np.transpose([pH_tigr21_10['pH_raw']]).flatten()
pH_tigr21_10_sg = savgol_filter(ls, window_length = 144*60, polyorder = 2, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr21_10['Day'], pH_tigr21_10_raw, label = 'Raw')
ax.plot(pH_tigr21_10['Day'], pH_tigr21_10_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr21_10',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()


# %%
###Data tigr21_11
pH_tigr21_11 = pd.read_excel("pH_result.xlsx",sheet_name='tigr21_11')
pH_tigr21_11.keys()
pH_tigr21_11_raw = pH_tigr21_11['pH_raw']


ls = np.transpose([pH_tigr21_11['pH_raw']]).flatten()
pH_tigr21_11_sg = savgol_filter(ls, window_length = 144*60, polyorder = 2, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr21_11['Day'], pH_tigr21_11_raw, label = 'Raw')
ax.plot(pH_tigr21_11['Day'], pH_tigr21_11_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr21_11',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()


# %%
###Data tigr21_12
pH_tigr21_12 = pd.read_excel("pH_result.xlsx",sheet_name='tigr21_12')
pH_tigr21_12.keys()
pH_tigr21_12_raw = pH_tigr21_12['pH_raw']


ls = np.transpose([pH_tigr21_12['pH_raw']]).flatten()
pH_tigr21_12_sg = savgol_filter(ls, window_length = 144*60, polyorder = 2, deriv=0)


fig, ax = plt.subplots()
ax.plot(pH_tigr21_12['Day'], pH_tigr21_12_raw, label = 'Raw')
ax.plot(pH_tigr21_12['Day'], pH_tigr21_12_sg, label = 'SG', linewidth=3)
ax.set_title('pH_tigr21_12',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
ax.set_ylim([6.6, 7.8])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%% 
### Combine all the plot into 3 * 3 plot, for tigr21_10, tigr21_11, tigr21_12, save every 60 points
# List of datasets and titles
datasets = [
    (pH_tigr15_1['Day'][pH_tigr15_1['Day'] <= 14], pH_tigr15_1_raw[pH_tigr15_1['Day'] <= 14], pH_tigr15_1_sg[pH_tigr15_1['Day'] <= 14], 'pH_tigr15_1'),
    (pH_tigr15_2['Day'][pH_tigr15_2['Day'] <= 14], pH_tigr15_2_raw[pH_tigr15_2['Day'] <= 14], pH_tigr15_2_sg[pH_tigr15_2['Day'] <= 14], 'pH_tigr15_2'),
    (pH_tigr15_3['Day'][pH_tigr15_3['Day'] <= 14], pH_tigr15_3_raw[pH_tigr15_3['Day'] <= 14], pH_tigr15_3_sg[pH_tigr15_3['Day'] <= 14], 'pH_tigr15_3'),
    (pH_tigr20_6['Day'][pH_tigr20_6['Day'] <= 14], pH_tigr20_6_raw[pH_tigr20_6['Day'] <= 14], pH_tigr20_6_sg[pH_tigr20_6['Day'] <= 14], 'pH_tigr20_6'),
    (pH_tigr20_7['Day'][pH_tigr20_7['Day'] <= 14], pH_tigr20_7_raw[pH_tigr20_7['Day'] <= 14], pH_tigr20_7_sg[pH_tigr20_7['Day'] <= 14], 'pH_tigr20_7'),
    (pH_tigr20_8['Day'][pH_tigr20_8['Day'] <= 14], pH_tigr20_8_raw[pH_tigr20_8['Day'] <= 14], pH_tigr20_8_sg[pH_tigr20_8['Day'] <= 14], 'pH_tigr20_8'),
    (pH_tigr21_10['Day'][::60][pH_tigr21_10['Day'][::60] <= 14], pH_tigr21_10_raw[::60][pH_tigr21_10['Day'][::60] <= 14], pH_tigr21_10_sg[::60][pH_tigr21_10['Day'][::60] <= 14], 'pH_tigr21_10'),
    (pH_tigr21_11['Day'][::60][pH_tigr21_11['Day'][::60] <= 14], pH_tigr21_11_raw[::60][pH_tigr21_11['Day'][::60] <= 14], pH_tigr21_11_sg[::60][pH_tigr21_11['Day'][::60] <= 14], 'pH_tigr21_11'),
    (pH_tigr21_12['Day'][::60][pH_tigr21_12['Day'][::60] <= 14], pH_tigr21_12_raw[::60][pH_tigr21_12['Day'][::60] <= 14], pH_tigr21_12_sg[::60][pH_tigr21_12['Day'][::60] <= 14], 'pH_tigr21_12'),
]

#%%
import matplotlib.colors as mcolors

# Use colorblind-friendly palette
colors = ['#999999', '#B22222']  # gray for raw, dark orange for smoothed

# Define colors for smoothed lines by row
smoothed_colors = ['blue', 'red', 'green']  # blue, red, green

fig, axs = plt.subplots(3, 3, figsize=(12, 9))
for idx, (day, raw, sg, _) in enumerate(datasets):
    row, col = divmod(idx, 3)
    ax = axs[row, col]
    ax.plot(day, raw, label='Raw', color=colors[0], linewidth=3)
    ax.plot(day, sg, label='Smoothed', color=smoothed_colors[row], linewidth=3)
    ax.set_xlabel('Time (day)', fontsize=15)
    ax.set_ylim([6.6, 7.8])
    ax.set_xlim([-0.2, 14.2])
    ax.set_xticks(np.arange(0, 16, 1))
    ax.legend(loc='best', fontsize=12)
    ax.tick_params(axis='both', labelsize=12)
    ax.text(0.01, 0.98, f'({chr(97+idx)})', transform=ax.transAxes, fontsize=15, fontweight='bold', va='top')
    ax.grid(alpha=0.3)  # add grid
    # Only show y-tick values for first column
    if col != 0:
        ax.set_yticklabels([])
for ax in axs[0:3, 0]:
    ax.set_ylabel('pH', fontsize=15)
plt.tight_layout()
plt.savefig('pH_Profile.png')
plt.show()



#%%
days = np.arange(0, 15)
def feed_profile(days, late_day12):
    return np.piecewise(days,
                        [days < 3, (days >= 3) & (days <= 5),
                         (days >= 6) & (days <= 7),
                         (days >= 8) & (days <= 9),
                         (days >= 10) & (days <= 11),
                         days >= 12],
                        [0, 3, 4, 5, 4, late_day12])

feed_A = feed_profile(days, 3)
feed_B = feed_profile(days, 4)
feed_C = feed_A.copy()

cases = [('Case A', feed_A, 'blue'),
         ('Case B', feed_B, 'red'),
         ('Case C', feed_C, 'green')]

fig, axs = plt.subplots(1, 3, figsize=(10, 3), sharey=True, sharex=True)
for ax, (label, feed, color) in zip(axs, cases):
    ax.step(days, feed, where='post', color=color, linewidth=3)
    #ax.set_title(label, fontsize=12, pad=6)
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.2, 14.2)
    ax.set_xticks(np.arange(0, 15, 1))
    ax.set_ylim(0, 6)
    ax.tick_params(labelsize=10)
    ax.set_xlabel('Time (day)', fontsize=12)

axs[0].set_ylabel('Feed: Cell Boost 7a (Cytiva)', fontsize=12)
#axs[1].set_xlabel('Culture Day', fontsize=12)
#fig.suptitle('Feeding Strategies by Case', fontsize=14, fontweight='bold', y=1.05)
fig.tight_layout()
for idx, ax in enumerate(axs):
    ax.text(0.02, 0.95, f'({chr(97+idx)})', transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')
plt.savefig('Feeding_Strategies.png', bbox_inches='tight')
plt.show()


# %%
days = np.arange(0, 15)
def feed_profile(days, late_day12):
    return np.piecewise(days,
                        [days < 3, (days >= 3) & (days <= 5),
                         (days >= 6) & (days <= 7),
                         (days >= 8) & (days <= 9),
                         (days >= 10) & (days <= 11),
                         days >= 12],
                        [0, 3, 4, 5, 4, late_day12])

feed_A = feed_profile(days, 3)
feed_B = feed_profile(days, 4)
feed_C = feed_A.copy()

cases = [('Case A', feed_A, 'blue'),
         ('Case B', feed_B, 'red'),
         ('Case C', feed_C, 'green')]

fig, axs = plt.subplots(1, 3, figsize=(10, 3), sharey=True, sharex=True)
bar_width = 0.4
for ax, (label, feed, color) in zip(axs, cases):
    ax.bar(days, feed, width=bar_width, color=color, edgecolor='black')
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.2, 14.4)
    ax.set_xticks(np.arange(0, 15, 1))
    ax.set_ylim(0, 6)
    ax.tick_params(labelsize=10)
    ax.set_xlabel('Time (day)', fontsize=12)

axs[0].set_ylabel('Feed Volume (% v/v)\nCell Boost 7a (Cytiva)', fontsize=12)
fig.tight_layout()
for idx, ax in enumerate(axs):
    ax.text(0.02, 0.95, f'({chr(97+idx)})', transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')
plt.savefig('Feeding_Strategies_Bar.png', bbox_inches='tight')
plt.show()


# %%
### save the smoothed data with day to excel, for tigr21_10, tigr21_11, tigr21_12, save every 60 points

with pd.ExcelWriter('pH_result_sg.xlsx') as writer:
    # tigr15_1
    df = pd.DataFrame({'Day': pH_tigr15_1['Day'], 'pH_sg': pH_tigr15_1_sg})
    df.to_excel(writer, sheet_name='tigr15_1', index=False)

    # tigr15_2
    df = pd.DataFrame({'Day': pH_tigr15_2['Day'], 'pH_sg': pH_tigr15_2_sg})
    df.to_excel(writer, sheet_name='tigr15_2', index=False)

    # tigr15_3
    df = pd.DataFrame({'Day': pH_tigr15_3['Day'], 'pH_sg': pH_tigr15_3_sg})
    df.to_excel(writer, sheet_name='tigr15_3', index=False)

    # tigr20_6
    df = pd.DataFrame({'Day': pH_tigr20_6['Day'], 'pH_sg': pH_tigr20_6_sg})
    df.to_excel(writer, sheet_name='tigr20_6', index=False)

    # tigr20_7
    df = pd.DataFrame({'Day': pH_tigr20_7['Day'], 'pH_sg': pH_tigr20_7_sg})
    df.to_excel(writer, sheet_name='tigr20_7', index=False)

    # tigr20_8
    df = pd.DataFrame({'Day': pH_tigr20_8['Day'], 'pH_sg': pH_tigr20_8_sg})
    df.to_excel(writer, sheet_name='tigr20_8', index=False)

    # tigr21_10 (every 60th point)
    df = pd.DataFrame({'Day': pH_tigr21_10['Day'][::60].reset_index(drop=True), 'pH_sg': pH_tigr21_10_sg[::60]})
    df.to_excel(writer, sheet_name='tigr21_10', index=False)

    # tigr21_11 (every 60th point)
    df = pd.DataFrame({'Day': pH_tigr21_11['Day'][::60].reset_index(drop=True), 'pH_sg': pH_tigr21_11_sg[::60]})
    df.to_excel(writer, sheet_name='tigr21_11', index=False)

    # tigr21_12 (every 60th point)
    df = pd.DataFrame({'Day': pH_tigr21_12['Day'][::60].reset_index(drop=True), 'pH_sg': pH_tigr21_12_sg[::60]})
    df.to_excel(writer, sheet_name='tigr21_12', index=False)
# %%
