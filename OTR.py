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
###Data tigr15_1_SLMB
OTR_tigr15_1 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr15_1_SLMB')
OTR_tigr15_1.keys()
OTR_tigr15_1_raw = OTR_tigr15_1['OTR_raw']

  
ls = np.transpose([OTR_tigr15_1['OTR_raw']]).flatten()
OTR_tigr15_1_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr15_1['Day'], OTR_tigr15_1_raw, label = 'Raw')
ax.plot(OTR_tigr15_1['Day'], OTR_tigr15_1_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr15_1',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr15_2_SLMB
OTR_tigr15_2 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr15_2_SLMB')
OTR_tigr15_2.keys()
OTR_tigr15_2_raw = OTR_tigr15_2['OTR_raw']

  
ls = np.transpose([OTR_tigr15_2['OTR_raw']]).flatten()
OTR_tigr15_2_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr15_2['Day'], OTR_tigr15_2_raw, label = 'Raw')
ax.plot(OTR_tigr15_2['Day'], OTR_tigr15_2_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr15_2',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr15_3_SLMB
OTR_tigr15_3 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr15_3_SLMB')
OTR_tigr15_3.keys()
OTR_tigr15_3_raw = OTR_tigr15_3['OTR_raw']

  
ls = np.transpose([OTR_tigr15_3['OTR_raw']]).flatten()
OTR_tigr15_3_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr15_3['Day'], OTR_tigr15_3_raw, label = 'Raw')
ax.plot(OTR_tigr15_3['Day'], OTR_tigr15_3_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr15_3',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr20_6_SLMB
OTR_tigr20_6 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr20_6_SLMB')
OTR_tigr20_6.keys()
OTR_tigr20_6_raw = OTR_tigr20_6['OTR_raw']

  
ls = np.transpose([OTR_tigr20_6['OTR_raw']]).flatten()
OTR_tigr20_6_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr20_6['Day'], OTR_tigr20_6_raw, label = 'Raw')
ax.plot(OTR_tigr20_6['Day'], OTR_tigr20_6_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr20_6',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr20_7_SLMB
OTR_tigr20_7 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr20_7_SLMB')
OTR_tigr20_7.keys()
OTR_tigr20_7_raw = OTR_tigr20_7['OTR_raw']

  
ls = np.transpose([OTR_tigr20_7['OTR_raw']]).flatten()
OTR_tigr20_7_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr20_7['Day'], OTR_tigr20_7_raw, label = 'Raw')
ax.plot(OTR_tigr20_7['Day'], OTR_tigr20_7_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr20_7',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

#%%
###Data tigr20_8_SLMB
OTR_tigr20_8 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr20_8_SLMB')
OTR_tigr20_8.keys()
OTR_tigr20_8_raw = OTR_tigr20_8['OTR_raw']

  
ls = np.transpose([OTR_tigr20_8['OTR_raw']]).flatten()
OTR_tigr20_8_sg = savgol_filter(ls, window_length = 144, polyorder = 1, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr20_8['Day'], OTR_tigr20_8_raw, label = 'Raw')
ax.plot(OTR_tigr20_8['Day'], OTR_tigr20_8_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr20_8',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()


# %%
###Data tigr21_10_SLMB
OTR_tigr21_10 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr21_10_SLMB')
OTR_tigr21_10.keys()
OTR_tigr21_10_raw = OTR_tigr21_10['OTR_raw']


ls = np.transpose([OTR_tigr21_10['OTR_raw']]).flatten()
OTR_tigr21_10_sg = savgol_filter(ls, window_length = 144*60, polyorder = 3, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr21_10['Day'], OTR_tigr21_10_raw, label = 'Raw')
ax.plot(OTR_tigr21_10['Day'], OTR_tigr21_10_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr21_10',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()


# %%
###Data tigr21_11_SLMB
OTR_tigr21_11 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr21_11_SLMB')
OTR_tigr21_11.keys()
OTR_tigr21_11_raw = OTR_tigr21_11['OTR_raw']


ls = np.transpose([OTR_tigr21_11['OTR_raw']]).flatten()
OTR_tigr21_11_sg = savgol_filter(ls, window_length = 144*60, polyorder = 3, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr21_11['Day'], OTR_tigr21_11_raw, label = 'Raw')
ax.plot(OTR_tigr21_11['Day'], OTR_tigr21_11_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr21_11',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

# %%
###Data tigr21_12_SLMB
OTR_tigr21_12 = pd.read_excel("OTR_result.xlsx",sheet_name='tigr21_12_SLMB')
OTR_tigr21_12.keys()
OTR_tigr21_12_raw = OTR_tigr21_12['OTR_raw']


ls = np.transpose([OTR_tigr21_12['OTR_raw']]).flatten()
OTR_tigr21_12_sg = savgol_filter(ls, window_length = 144*60, polyorder = 3, deriv=0)


fig, ax = plt.subplots()
ax.plot(OTR_tigr21_12['Day'], OTR_tigr21_12_raw, label = 'Raw')
ax.plot(OTR_tigr21_12['Day'], OTR_tigr21_12_sg, label = 'SG', linewidth=3)
ax.set_title('OTR_tigr21_12',fontsize=15)
ax.set_xlabel('Time (Day)',fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)    
leg = ax.legend(loc = 'best')
# plt.savefig('OTR.svg')
plt.show()

# %%
