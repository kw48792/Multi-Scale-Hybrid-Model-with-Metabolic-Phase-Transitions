#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters, Parameter, report_fit,  fit_report
from scipy.integrate import odeint
from scipy import interpolate
from scipy.optimize import minimize as sp_minimize
import scipy
import time
import pickle
from lmfit import minimize, Parameters  
from scipy.stats import mode
import math
import os
import random

#%%
#Dry weight of CHO cells is 252.3 pg/cell 
#Biomass of 1 million cells/mL is 252.3 ng/mL
biomass = 252.3 #ng/mL

Delta_t = 0.1 #Days

#%%
###Data
data_extra = pd.read_excel("Tigr6to20.xlsx", sheet_name=None)

data_extra.keys()
for key in data_extra.keys():
    print(key)

selected_columns = ['Day',
    "Glucose_Measure (mM)", "Glutamine_Measure (mM)", "Glutamic Acid_Measure (mM)",
    "Aspartic Acid_Measure (mM)", "Alanine_Measure (mM)", "Lactate_Measure (mM)",
    "Serine_Measure (mM)", "Leucine_Measure (mM)", "Isoleucine_Measure (mM)",
    "Valine_Measure (mM)", "NH3_Measure (mM)", "IgG_Measure (mM)", "VCD_Measure (million cells/mL)"
]


data_extra_cont =  {k:data_extra[k][selected_columns] for k in ['tigr15_1', 'tigr15_2', 'tigr15_3', 'tigr20_6', 'tigr20_7', 'tigr20_8', 
                                                                'tigr21_10', 'tigr21_11', 'tigr21_12']}


### Concentration Changing due to Feeding and Dilution
feeding_cont = pd.read_excel("Tigr6to20_Feeding.xlsx", sheet_name=None)
feeding_cont =  {k:feeding_cont[k][selected_columns] for k in ['tigr15_1', 'tigr15_2', 'tigr15_3', 'tigr20_6', 'tigr20_7', 'tigr20_8', 
                                                                'tigr21_10', 'tigr21_11', 'tigr21_12']}

for key in ['tigr15_1', 'tigr15_2', 'tigr15_3', 'tigr20_6', 'tigr20_7', 'tigr20_8', 'tigr21_10', 'tigr21_11', 'tigr21_12']:
    feeding_cont[key]['Day'] = np.ceil(feeding_cont[key]['Day']) + Delta_t/100 # Feeding happens at the end of the day
    data_extra_cont[key]['Day'] = np.ceil(data_extra_cont[key]['Day'])

### Drop the rows when Glucose_Measure (mM) = 0 in feeding_cont
for key in ['tigr15_1', 'tigr15_2', 'tigr15_3', 'tigr20_6', 'tigr20_7', 'tigr20_8', 'tigr21_10', 'tigr21_11', 'tigr21_12']:
    feeding_cont[key] = feeding_cont[key].loc[feeding_cont[key]['Glucose_Measure (mM)'] != 0, :]


# Reorder and rename feeding and measurement data to the requested extracellular/metabolite keys
_meta_order = ['EGLC', 'EGLN', 'EGLU', 'EASP', 'EALA', 'ELAC',
                'ESER', 'ELEU', 'EILE', 'EVAL', 'ENH4', 'EANTI', 'VCD_0']

# Map meta keys to the original column names present in the loaded sheets
_col_map = {
    'EGLC': 'Glucose_Measure (mM)',
    'EGLN': 'Glutamine_Measure (mM)',
    'EGLU': 'Glutamic Acid_Measure (mM)',
    'EASP': 'Aspartic Acid_Measure (mM)',
    'EALA': 'Alanine_Measure (mM)',
    'ELAC': 'Lactate_Measure (mM)',
    'ESER': 'Serine_Measure (mM)',
    'ELEU': 'Leucine_Measure (mM)',
    'EILE': 'Isoleucine_Measure (mM)',
    'EVAL': 'Valine_Measure (mM)',
    'ENH4': 'NH3_Measure (mM)',
    'EANTI': 'IgG_Measure (mM)',
    'VCD_0': 'VCD_Measure (million cells/mL)',
    'VCD': 'VCD_Measure (million cells/mL)'
}

def _reorder_and_rename(df):
    # ensure Day column kept, create new DF with ordered keys (add NaN if missing)
    out = pd.DataFrame()
    out['Day'] = df['Day'].values
    for key in _meta_order:
        orig = _col_map.get(key)
        if orig in df.columns:
            out[key] = df[orig].values
        else:
            out[key] = np.nan
    return out

# Apply to all batches
data_extra_cont = {k: _reorder_and_rename(v).copy() for k, v in data_extra_cont.items()}
feeding_cont = {k: _reorder_and_rename(v).copy() for k, v in feeding_cont.items()}

#%%
### Loading day specific OUR data
tigr_our_all = pd.read_excel("Rate_v2.xlsx", sheet_name='OUR') 
# replacing the column names with ['tigr15_1', 'tigr15_2', 'tigr15_3', 'tigr20_6', 'tigr20_7', 'tigr20_8']
tigr_our_all.columns = ['Day', 'tigr15_1', 'tigr15_2', 'tigr15_3', 'tigr20_6', 'tigr20_7', 'tigr20_8', 'tigr21_10', 'tigr21_11', 'tigr21_12']
# Remove the duplicate Day column  
tigr_our_all = tigr_our_all.loc[~tigr_our_all['Day'].duplicated(), :]
### Loading pH data
tigr_pH_all = pd.read_excel("pH_result_sg.xlsx", sheet_name=None)

for key in tigr_pH_all.keys():
    print(key)


#%%
###Stoichiometry Matrix
N = pd.read_excel('Simulator.xlsx', sheet_name="N")
N = N.fillna(0).values
N = N[:,1:]
N = N.astype(float).reshape(32,30)

# Create N_3 as a copy of N and zero out the 30th column (index 29)
N_3 = N.copy()
N_3[:, 29] = 0

#%%
### Plot the metabolism data
for i in data_extra_cont['tigr15_1'].columns[2:]:
    fig, ax = plt.subplots()
    for key in data_extra_cont.keys():    
        ax.plot(data_extra_cont[key]['Day'], data_extra_cont[key][i], label=key)
        ax.set_title(i,fontsize=15)
    ax.set_xlabel('Time (Day)',fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)    
    leg = ax.legend(loc = 'best')
    # plt.savefig(str(meta_full[i]) + '(' + dataset[test_index] + ')' + '.svg')
    plt.show()


#%%
### Plot the OUR data in one figure 
fig, ax = plt.subplots(figsize=(8, 5))
for i in tigr_our_all.columns[1:]:
    ax.plot(tigr_our_all['Day'], tigr_our_all[i], marker='o', label=i)
ax.set_title('OUR for All Batches', fontsize=15)
ax.set_xlabel('Time (Day)', fontsize=15)
ax.set_ylabel('OUR', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.legend(loc='best')
plt.tight_layout()
plt.show()


#%%# Plot the pH data in one figure
fig, ax = plt.subplots(figsize=(8, 5))
for i in tigr_pH_all.keys():
    ax.plot(tigr_pH_all[i]['Day'], tigr_pH_all[i]['pH_sg'], marker='o', label=i)
ax.set_title('pH for All Batches', fontsize=15)
ax.set_xlabel('Time (Day)', fontsize=15)
ax.set_ylabel('pH', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.legend(loc='best')
plt.tight_layout()
plt.show()


#%%
### Based on the input time and batch to calculate/interpolate the OUR
def interpolate_our(t, batch):
    """
    Interpolates the OUR at time t for the given key.

    Args:
        t (float): Time.
        batch (str): Batch.

    Returns:
        float: Interpolated OUR at time t.
    """
    # Interpolates the OUR at time t for the given key
    our = tigr_our_all[batch].dropna()
    time = tigr_our_all['Day'][our.index]
    f = interpolate.interp1d(time, our, kind='linear', bounds_error=False, fill_value='extrapolate')
    return f(t)

# interpolate_our(3.5,'tigr15_1')


### Based on the input time and batch to calculate/interpolate the pH
def interpolate_pH(t, batch):
    """
    Interpolates the pH at time t for the given batch.

    Args:
        t (float): Time.
        batch (str): Batch.

    Returns:
        float: Interpolated pH at time t.
    """
    pH = tigr_pH_all[batch].dropna()
    time = tigr_pH_all[batch]['Day'][pH.index]
    f = interpolate.interp1d(time, pH['pH_sg'], kind='linear', bounds_error=False, fill_value='extrapolate')
    return f(t)

# interpolate_pH(3.5,'tigr15_1')


#%%
meta_list = ['AcCoA', 'AKG', 'CO2', 'G6P', 'GLU', 'SUC', 'MAL', 'OAA', 'PYR', 'ALA', 'ASP', 'LAC', 'GLN', 
             'LEU', 'ILE', 'VAL', 'SER', 'NH4', 'ANTI', 'EGLC', 'EGLN', 'EGLU', 'EASP', 'EALA', 'ELAC', 
             'ESER', 'ELEU', 'EILE', 'EVAL', 'ENH4', 'EANTI', 'BIOM_0', 'BIOM_1', 'BIOM_2', 'BIOM_3', 'VCD_0', 
             'VCD_1', 'VCD_2', 'VCD_3']


target_label_mapper = dict(zip(meta_list, range(len(meta_list))))



#%%
##############################################################################################################
###############################################Flux Rate Estimation###########################################
#############################################################################################################
def f_m(t, xs, ps, mode, batch):
    """
    Bio_Kinetic Model.
    
    This function calculates the rates of biochemical reactions in a bio-kinetic model.
    
    Parameters:
        t (float): Time.
        xs (list): List of state variables.
        ps (dict): Dictionary of parameter values.
        mode (binary): store fluxes (1) or not (0).
        
    Returns:
        tuple: Tuple containing the rates of biochemical reactions in the growth phase and stationary phase.
    """ 

    ##Early Growth Phase
    vmax2_0 = ps['vmax2_0'].value
    vmax3f_0 = ps['vmax3f_0'].value
    vmax3r_0 = ps['vmax3r_0'].value
    vmax5f_0 = ps['vmax5f_0'].value
    vmax5r_0 = ps['vmax5r_0'].value
    vmax7_0 = ps['vmax7_0'].value
    vmax16f_0 = ps['vmax16f_0'].value
    vmax16r_0 = ps['vmax16r_0'].value
    vmax18f_0 = ps['vmax18f_0'].value
    vmax18r_0 = ps['vmax18r_0'].value
    vmax19_0 = ps['vmax19_0'].value
    # vmax19r_0 = ps['vmax19r_0'].value
    vmax21_0 = ps['vmax21_0'].value
    vmax23_0 = ps['vmax23_0'].value
    vmax25_0 = ps['vmax25_0'].value
    vmax27f_0 = ps['vmax27f_0'].value
    vmax27r_0 = ps['vmax27r_0'].value
    vmax28_0 = ps['vmax28_0'].value
    vmax30_0 = ps['vmax30_0'].value
    
    K_mEGLC_0 = ps['K_mEGLC_0'].value
    K_iLACtoHK_0 = ps['K_iLACtoHK_0'].value
    K_iLACtoGLNS_0 = ps['K_iLACtoGLNS_0'].value
    K_aGLN_0 = ps['K_aGLN_0'].value
    K_mELAC_0 = ps['K_mELAC_0'].value
    K_mEALA_0 = ps['K_mEALA_0'].value
    K_mESER_0 = ps['K_mESER_0'].value
    K_mEGLN_0 = ps['K_mEGLN_0'].value
    K_mNH4_0 = ps['K_mNH4_0'].value
    K_mENH4_0 = ps['K_mENH4_0'].value
    K_mEGLU_0 = ps['K_mEGLU_0'].value
    K_mEASP_0 = ps['K_mEASP_0'].value
    K_mELEU_0 = ps['K_mELEU_0'].value
    K_mEILE_0 = ps['K_mEILE_0'].value
    K_mEVAL_0 = ps['K_mEVAL_0'].value


    K_mEGLNtoIgG_0 = ps['K_mEGLNtoIgG_0'].value
    K_mEGLUtoIgG_0 = ps['K_mEGLUtoIgG_0'].value
    K_mEASPtoIgG_0 = ps['K_mEASPtoIgG_0'].value
    K_mEALAtoIgG_0 = ps['K_mEALAtoIgG_0'].value
    K_mESERtoIgG_0 = ps['K_mESERtoIgG_0'].value
    K_mELEUtoIgG_0 = ps['K_mELEUtoIgG_0'].value
    K_mEILEtoIgG_0 = ps['K_mEILEtoIgG_0'].value
    K_mEVALtoIgG_0 = ps['K_mEVALtoIgG_0'].value

    
    K_mEGLCtoBIO_0 = ps['K_mEGLCtoBIO_0'].value
    K_mEGLNtoBIO_0 = ps['K_mEGLNtoBIO_0'].value 
    K_mEGLUtoBIO_0 = ps['K_mEGLUtoBIO_0'].value
    K_mEASPtoBIO_0 = ps['K_mEASPtoBIO_0'].value
    K_mEALAtoBIO_0 = ps['K_mEALAtoBIO_0'].value
    K_mESERtoBIO_0 = ps['K_mESERtoBIO_0'].value
    K_mELEUtoBIO_0 = ps['K_mELEUtoBIO_0'].value
    K_mEILEtoBIO_0 = ps['K_mEILEtoBIO_0'].value
    K_mEVALtoBIO_0 = ps['K_mEVALtoBIO_0'].value
   
    k_0 = ps['k_0'].value
    K1_0 = ps['K1_0'].value
    K2_0 = ps['K2_0'].value
    
    ##Late Growth Phase
    vmax2_1 = ps['vmax2_1'].value
    vmax3f_1 = ps['vmax3f_1'].value
    vmax3r_1 = ps['vmax3r_1'].value
    vmax5f_1 = ps['vmax5f_1'].value
    vmax5r_1 = ps['vmax5r_1'].value
    vmax7_1 = ps['vmax7_1'].value
    vmax16f_1 = ps['vmax16f_1'].value
    vmax16r_1 = ps['vmax16r_1'].value
    vmax18f_1 = ps['vmax18f_1'].value
    vmax18r_1 = ps['vmax18r_1'].value
    vmax19_1 = ps['vmax19_1'].value
    # vmax19r_1 = ps['vmax19r_1'].value
    vmax21_1 = ps['vmax21_1'].value
    vmax23_1 = ps['vmax23_1'].value
    vmax25_1 = ps['vmax25_1'].value
    vmax27f_1 = ps['vmax27f_1'].value
    vmax27r_1 = ps['vmax27r_1'].value
    vmax28_1 = ps['vmax28_1'].value
    vmax30_1 = ps['vmax30_1'].value

    K_mEGLC_1 = ps['K_mEGLC_1'].value
    K_iLACtoHK_1 = ps['K_iLACtoHK_1'].value
    K_iLACtoGLNS_1 = ps['K_iLACtoGLNS_1'].value
    K_aGLN_1 = ps['K_aGLN_1'].value
    K_mELAC_1 = ps['K_mELAC_1'].value
    K_mEALA_1 = ps['K_mEALA_1'].value
    K_mESER_1 = ps['K_mESER_1'].value
    K_mEGLN_1 = ps['K_mEGLN_1'].value
    K_mNH4_1 = ps['K_mNH4_1'].value
    K_mENH4_1 = ps['K_mENH4_1'].value
    K_mEGLU_1 = ps['K_mEGLU_1'].value
    K_mEASP_1 = ps['K_mEASP_1'].value
    K_mELEU_1 = ps['K_mELEU_1'].value
    K_mEILE_1 = ps['K_mEILE_1'].value
    K_mEVAL_1 = ps['K_mEVAL_1'].value

    K_mEGLNtoIgG_1 = ps['K_mEGLNtoIgG_1'].value
    K_mEGLUtoIgG_1 = ps['K_mEGLUtoIgG_1'].value
    K_mEASPtoIgG_1 = ps['K_mEASPtoIgG_1'].value
    K_mEALAtoIgG_1 = ps['K_mEALAtoIgG_1'].value
    K_mESERtoIgG_1 = ps['K_mESERtoIgG_1'].value
    K_mELEUtoIgG_1 = ps['K_mELEUtoIgG_1'].value
    K_mEILEtoIgG_1 = ps['K_mEILEtoIgG_1'].value
    K_mEVALtoIgG_1 = ps['K_mEVALtoIgG_1'].value

    
    K_mEGLCtoBIO_1 = ps['K_mEGLCtoBIO_1'].value
    K_mEGLNtoBIO_1 = ps['K_mEGLNtoBIO_1'].value
    K_mEGLUtoBIO_1 = ps['K_mEGLUtoBIO_1'].value
    K_mEASPtoBIO_1 = ps['K_mEASPtoBIO_1'].value
    K_mEALAtoBIO_1 = ps['K_mEALAtoBIO_1'].value
    K_mESERtoBIO_1 = ps['K_mESERtoBIO_1'].value
    K_mELEUtoBIO_1 = ps['K_mELEUtoBIO_1'].value
    K_mEILEtoBIO_1 = ps['K_mEILEtoBIO_1'].value
    K_mEVALtoBIO_1 = ps['K_mEVALtoBIO_1'].value

    ##Stationary Phase
    vmax2_2 = ps['vmax2_2'].value
    vmax3f_2 = ps['vmax3f_2'].value
    vmax3r_2 = ps['vmax3r_2'].value
    vmax5f_2 = ps['vmax5f_2'].value
    vmax5r_2 = ps['vmax5r_2'].value
    vmax7_2 = ps['vmax7_2'].value
    vmax16f_2 = ps['vmax16f_2'].value
    vmax16r_2 = ps['vmax16r_2'].value
    vmax18f_2 = ps['vmax18f_2'].value
    vmax18r_2 = ps['vmax18r_2'].value
    vmax19_2 = ps['vmax19_2'].value
    # vmax19r_2 = ps['vmax19r_2'].value
    vmax21_2 = ps['vmax21_2'].value
    vmax23_2 = ps['vmax23_2'].value
    vmax25_2 = ps['vmax25_2'].value
    vmax27f_2 = ps['vmax27f_2'].value
    vmax27r_2 = ps['vmax27r_2'].value
    vmax28_2 = ps['vmax28_2'].value
    vmax30_2 = ps['vmax30_2'].value

    K_mEGLC_2 = ps['K_mEGLC_2'].value
    K_iLACtoHK_2 = ps['K_iLACtoHK_2'].value
    K_iLACtoGLNS_2 = ps['K_iLACtoGLNS_2'].value
    K_aGLN_2 = ps['K_aGLN_2'].value
    K_mELAC_2 = ps['K_mELAC_2'].value
    K_mEALA_2 = ps['K_mEALA_2'].value
    K_mESER_2 = ps['K_mESER_2'].value
    K_mEGLN_2 = ps['K_mEGLN_2'].value
    K_mNH4_2 = ps['K_mNH4_2'].value
    K_mENH4_2 = ps['K_mENH4_2'].value
    K_mEGLU_2 = ps['K_mEGLU_2'].value
    K_mEASP_2 = ps['K_mEASP_2'].value
    K_mELEU_2 = ps['K_mELEU_2'].value
    K_mEILE_2 = ps['K_mEILE_2'].value
    K_mEVAL_2 = ps['K_mEVAL_2'].value

    K_mEGLNtoIgG_2 = ps['K_mEGLNtoIgG_2'].value
    K_mEGLUtoIgG_2 = ps['K_mEGLUtoIgG_2'].value
    K_mEASPtoIgG_2 = ps['K_mEASPtoIgG_2'].value
    K_mEALAtoIgG_2 = ps['K_mEALAtoIgG_2'].value
    K_mESERtoIgG_2 = ps['K_mESERtoIgG_2'].value
    K_mELEUtoIgG_2 = ps['K_mELEUtoIgG_2'].value
    K_mEILEtoIgG_2 = ps['K_mEILEtoIgG_2'].value
    K_mEVALtoIgG_2 = ps['K_mEVALtoIgG_2'].value

    
    K_mEGLCtoBIO_2 = ps['K_mEGLCtoBIO_2'].value
    K_mEGLNtoBIO_2 = ps['K_mEGLNtoBIO_2'].value
    K_mEGLUtoBIO_2 = ps['K_mEGLUtoBIO_2'].value
    K_mEASPtoBIO_2 = ps['K_mEASPtoBIO_2'].value
    K_mEALAtoBIO_2 = ps['K_mEALAtoBIO_2'].value
    K_mESERtoBIO_2 = ps['K_mESERtoBIO_2'].value
    K_mELEUtoBIO_2 = ps['K_mELEUtoBIO_2'].value
    K_mEILEtoBIO_2 = ps['K_mEILEtoBIO_2'].value
    K_mEVALtoBIO_2 = ps['K_mEVALtoBIO_2'].value


    ##Decline Phase
    vmax2_3 = ps['vmax2_3'].value
    vmax3f_3 = ps['vmax3f_3'].value
    vmax3r_3 = ps['vmax3r_3'].value
    vmax5f_3 = ps['vmax5f_3'].value
    vmax5r_3 = ps['vmax5r_3'].value
    vmax7_3 = ps['vmax7_3'].value
    vmax16f_3 = ps['vmax16f_3'].value
    vmax16r_3 = ps['vmax16r_3'].value
    vmax18f_3 = ps['vmax18f_3'].value
    vmax18r_3 = ps['vmax18r_3'].value
    vmax19_3 = ps['vmax19_3'].value
    # vmax19r_2 = ps['vmax19r_2'].value
    vmax21_3 = ps['vmax21_3'].value
    vmax23_3 = ps['vmax23_3'].value
    vmax25_3 = ps['vmax25_3'].value
    vmax27f_3 = ps['vmax27f_3'].value
    vmax27r_3 = ps['vmax27r_3'].value
    vmax28_3 = ps['vmax28_3'].value
    vmax30_3 = ps['vmax30_3'].value

    K_mEGLC_3 = ps['K_mEGLC_3'].value
    K_iLACtoHK_3 = ps['K_iLACtoHK_3'].value
    K_iLACtoGLNS_3 = ps['K_iLACtoGLNS_3'].value
    K_aGLN_3 = ps['K_aGLN_3'].value
    K_mELAC_3 = ps['K_mELAC_3'].value
    K_mEALA_3 = ps['K_mEALA_3'].value
    K_mESER_3 = ps['K_mESER_3'].value
    K_mEGLN_3 = ps['K_mEGLN_3'].value
    K_mNH4_3 = ps['K_mNH4_3'].value
    K_mENH4_3 = ps['K_mENH4_3'].value
    K_mEGLU_3 = ps['K_mEGLU_3'].value
    K_mEASP_3 = ps['K_mEASP_3'].value
    K_mELEU_3 = ps['K_mELEU_3'].value
    K_mEILE_3 = ps['K_mEILE_3'].value
    K_mEVAL_3 = ps['K_mEVAL_3'].value

   
    K_mEGLNtoIgG_3 = ps['K_mEGLNtoIgG_3'].value
    K_mEGLUtoIgG_3 = ps['K_mEGLUtoIgG_3'].value
    K_mEASPtoIgG_3 = ps['K_mEASPtoIgG_3'].value
    K_mEALAtoIgG_3 = ps['K_mEALAtoIgG_3'].value
    K_mESERtoIgG_3 = ps['K_mESERtoIgG_3'].value
    K_mELEUtoIgG_3 = ps['K_mELEUtoIgG_3'].value
    K_mEILEtoIgG_3 = ps['K_mEILEtoIgG_3'].value
    K_mEVALtoIgG_3 = ps['K_mEVALtoIgG_3'].value

    
    K_mEGLCtoBIO_3 = ps['K_mEGLCtoBIO_3'].value
    K_mEGLNtoBIO_3 = ps['K_mEGLNtoBIO_3'].value
    K_mEGLUtoBIO_3 = ps['K_mEGLUtoBIO_3'].value
    K_mEASPtoBIO_3 = ps['K_mEASPtoBIO_3'].value
    K_mEALAtoBIO_3 = ps['K_mEALAtoBIO_3'].value
    K_mESERtoBIO_3 = ps['K_mESERtoBIO_3'].value
    K_mELEUtoBIO_3 = ps['K_mELEUtoBIO_3'].value
    K_mEILEtoBIO_3 = ps['K_mEILEtoBIO_3'].value
    K_mEVALtoBIO_3 = ps['K_mEVALtoBIO_3'].value


    AcCoA, AKG, CO2, G6P, GLU, SUC, MAL, OAA, PYR, ALA, ASP, LAC, GLN, LEU, ILE, VAL, SER, NH4, ANTI, \
        EGLC, EGLN, EGLU, EASP, EALA, ELAC, ESER, ELEU, EILE, EVAL, ENH4, EANTI, BIOM_0, BIOM_1, BIOM_2, BIOM_3, VCD_0, VCD_1, VCD_2, VCD_3 = xs
    
    ##Early Growth Phase
    v2_0 = k_0 / (1.0 + 10.0 ** (-interpolate_pH(t,batch)) / K1_0 + K2_0 / 10.0 ** (-interpolate_pH(t,batch))) * vmax2_0 * EGLC/(K_mEGLC_0 + EGLC)  * K_iLACtoHK_0/(K_iLACtoHK_0 + ELAC)
    # print(interpolate_pH(t,batch), k_0 / (1.0 + 10.0 ** (-interpolate_pH(t,batch)) / K1_0 + K2_0 / 10.0 ** (-interpolate_pH(t,batch))))
    # v2_0 = vmax2_0 * EGLC/(K_mEGLC_0 + EGLC)  * K_iLACtoHK_0/(K_iLACtoHK_0 + ELAC)

    # v3_0 = vmax3f_0 * EGLC/(K_mEGLC_0 * (1 + K_aGLN_0/EGLN) + EGLC) - vmax3r_0 * ELAC/(K_mELAC_0 + ELAC) 
    
    v3_0 = k_0 / (1.0 + 10.0 ** (-interpolate_pH(t,batch)) / K1_0 + K2_0 / 10.0 ** (-interpolate_pH(t,batch))) * vmax3f_0 * EGLC/(K_mEGLC_0 + EGLC) - vmax3r_0 * ELAC/(K_mELAC_0 + ELAC)

    v5_0 = k_0 / (1.0 + 10.0 ** (-interpolate_pH(t,batch)) / K1_0 + K2_0 / 10.0 ** (-interpolate_pH(t,batch))) * vmax5f_0 * EGLC/(K_mEGLC_0 + EGLC) - vmax5r_0 * EALA/(K_mEALA_0 + EALA)
        
    v7_0 = vmax7_0 * ESER/(K_mESER_0 + ESER)
    
    v16_0 = vmax16f_0 * EGLN/(K_mEGLN_0 + EGLN) * K_iLACtoGLNS_0/(K_iLACtoGLNS_0 + ELAC) - vmax16r_0 * EGLU/(K_mEGLU_0 + EGLU) # * NH4/(K_mNH4_0 + NH4)

    v18_0 = vmax18f_0 * EGLU/(K_mEGLU_0 + EGLU) - vmax18r_0 * EGLC/(K_mEGLC_0 + EGLC)

    v19_0 = vmax19_0 * EASP/(K_mEASP_0 + EASP) #- vmax19r_0 * EGLU/(K_mEGLU_0 + EGLU) * NH4/(K_mNH4_0 + NH4)
    
    v21_0 = vmax21_0 * ELEU/(K_mELEU_0 + ELEU)
    
    v23_0 = vmax23_0 * EVAL/(K_mEVAL_0 + EVAL)
    
    v25_0 = vmax25_0 * EILE/(K_mEILE_0 + EILE)

    v27_0 = vmax27f_0 * NH4/(K_mNH4_0 + NH4) - vmax27r_0 * ENH4/(K_mENH4_0 + ENH4)
    
    v28_0 = vmax28_0 * EGLN/(K_mEGLNtoIgG_0 + EGLN) * EGLU/(K_mEGLUtoIgG_0 + EGLU) * EALA/(K_mEALAtoIgG_0 + EALA) * EASP/(K_mEASPtoIgG_0 + EASP) * ESER/(K_mESERtoIgG_0 + ESER) * ELEU/(K_mELEUtoIgG_0 + ELEU) * EILE/(K_mEILEtoIgG_0 + EILE) * EVAL/(K_mEVALtoIgG_0 + EVAL)
    
    # v28_0 = 0.0008

    v30_0 = vmax30_0 * EGLC/(K_mEGLCtoBIO_0 + EGLC) * EGLN/(K_mEGLNtoBIO_0 + EGLN) * EGLU/(K_mEGLUtoBIO_0 + EGLU) * EALA/(K_mEALAtoBIO_0 + EALA) * EASP/(K_mEASPtoBIO_0 + EASP) * ESER/(K_mESERtoBIO_0 + ESER) * ELEU/(K_mELEUtoBIO_0 + ELEU) * EILE/(K_mEILEtoBIO_0 + EILE) * EVAL/(K_mEVALtoBIO_0 + EVAL)
    
    # v30_0 = 0.035

    v1_0 = v2_0
    
    v4_0 = v3_0
    
    v6_0 = v5_0
    
    v8_0 = v7_0
    
    v14_0 = v18_0 + v19_0 + v23_0 + v25_0
    
    v9_0 = 2 * v2_0 + v7_0 + v14_0 - v3_0 - v5_0
    
    v10_0 = v9_0 + 3 * v21_0 + v25_0
            
    v11_0 = v10_0 + v18_0
    
    v12_0 = v11_0 + v23_0 + v25_0
    
    v13_0 = v12_0 - v14_0
    
    v15_0 = v16_0
    
    v17_0 = v16_0 - v18_0 - v5_0 + v21_0 + v19_0 + v23_0 + v25_0
    
    v20_0 = v19_0
    
    v22_0 = v21_0
    
    v24_0 = v23_0
    
    v26_0 = v25_0
    
    v29_0 = v28_0
    
    # v30_0 = v29_0
    
    v_0= [v1_0, v2_0, v3_0, v4_0, v5_0, v6_0, v7_0, v8_0, v9_0, v10_0, v11_0, v12_0, v13_0, v14_0, v15_0, v16_0, v17_0, v18_0,	
        v19_0, v20_0, v21_0, v22_0, v23_0, v24_0, v25_0, v26_0, v27_0,v28_0, v29_0, v30_0]


    ##Late Exponential Phase
    v2_1 = vmax2_1 * EGLC/(K_mEGLC_1 + EGLC)  * K_iLACtoHK_1/(K_iLACtoHK_1 + ELAC)
    
    # v3_1 = vmax3f_1 * EGLC/(K_mEGLC_1 * (1 + K_aGLN_1/EGLN) + EGLC) - vmax3r_1 * ELAC/(K_mELAC_1 + ELAC) 

    v3_1 = vmax3f_1 * EGLC/(K_mEGLC_1 + EGLC) - vmax3r_1 * ELAC/(K_mELAC_1 + ELAC) 
      
    v5_1 = vmax5f_1 * EGLC/(K_mEGLC_1 + EGLC) - vmax5r_1 * EALA/(K_mEALA_1 + EALA)
     
    v7_1 = vmax7_1 * ESER/(K_mESER_1 + ESER)
     
    v16_1 = vmax16f_1 * EGLN/(K_mEGLN_1 + EGLN) * K_iLACtoGLNS_1/(K_iLACtoGLNS_1 + ELAC) - vmax16r_1 * EGLU/(K_mEGLU_1 + EGLU) # * NH4/(K_mNH4_1 + NH4)

    v18_1 = vmax18f_1 * EGLU/(K_mEGLU_1 + EGLU) - vmax18r_1 * EGLC/(K_mEGLC_1 + EGLC)

    v19_1 = vmax19_1 * EASP/(K_mEASP_1 + EASP) #- vmax19r_1 * EGLU/(K_mEGLU_1 + EGLU) * NH4/(K_mNH4_1 + NH4)
     
    v21_1 = vmax21_1 * ELEU/(K_mELEU_1 + ELEU)
     
    v23_1 = vmax23_1 * EVAL/(K_mEVAL_1 + EVAL)
     
    v25_1 = vmax25_1 * EILE/(K_mEILE_1 + EILE)

    v27_1 = vmax27f_1 * NH4/(K_mNH4_1 + NH4) - vmax27r_1 * ENH4/(K_mENH4_1 + ENH4)

    v28_1 = vmax28_1 * EGLN/(K_mEGLNtoIgG_1 + EGLN) * EGLU/(K_mEGLUtoIgG_1 + EGLU) * EALA/(K_mEALAtoIgG_1 + EALA) * EASP/(K_mEASPtoIgG_1 + EASP) * ESER/(K_mESERtoIgG_1 + ESER) * ELEU/(K_mELEUtoIgG_1 + ELEU) * EILE/(K_mEILEtoIgG_1 + EILE) * EVAL/(K_mEVALtoIgG_1 + EVAL)

    # v28_1 = 0.00042

    v30_1 = vmax30_1 * EGLC/(K_mEGLCtoBIO_1 + EGLC) * EGLN/(K_mEGLNtoBIO_1 + EGLN) * EGLU/(K_mEGLUtoBIO_1 + EGLU) * EALA/(K_mEALAtoBIO_1 + EALA) * EASP/(K_mEASPtoBIO_1 + EASP) * ESER/(K_mESERtoBIO_1 + ESER) * ELEU/(K_mELEUtoBIO_1 + ELEU) * EILE/(K_mEILEtoBIO_1 + EILE) * EVAL/(K_mEVALtoBIO_1 + EVAL)
    
    # v30_1 = 0.023

    v1_1 = v2_1
    
    v4_1 = v3_1
    
    v6_1 = v5_1
    
    v8_1 = v7_1
    
    v14_1 = v18_1 + v19_1 + v23_1 + v25_1
    
    v9_1 = 2 * v2_1 + v7_1 + v14_1 - v3_1 - v5_1
    
    v10_1 = v9_1 + 3 * v21_1 + v25_1
            
    v11_1 = v10_1 + v18_1
    
    v12_1 = v11_1 + v23_1 + v25_1
    
    v13_1 = v12_1 - v14_1
    
    v15_1 = v16_1
    
    v17_1 = v16_1 - v18_1 - v5_1 + v21_1 + v19_1 + v23_1 + v25_1
    
    v20_1 = v19_1
    
    v22_1 = v21_1
    
    v24_1 = v23_1
    
    v26_1 = v25_1
    
    v29_1 = v28_1
    
    # v30_1 = v29_1
    
    v_1= [v1_1, v2_1, v3_1, v4_1, v5_1, v6_1, v7_1, v8_1, v9_1, v10_1, v11_1, v12_1, v13_1, v14_1, v15_1, v16_1, v17_1, v18_1,	
        v19_1, v20_1, v21_1, v22_1, v23_1, v24_1, v25_1, v26_1, v27_1, v28_1, v29_1, v30_1]


    ##Stationary Phase
    v2_2 = vmax2_2 * EGLC/(K_mEGLC_2 + EGLC)  * K_iLACtoHK_2/(K_iLACtoHK_2 + ELAC)
    
    # v3_2 = vmax3f_2 * EGLC/(K_mEGLC_2 * (1 + K_aGLN_2/EGLN) + EGLC) - vmax3r_2 * ELAC/(K_mELAC_2 + ELAC) 

    v3_2 = vmax3f_2 * EGLC/(K_mEGLC_2 + EGLC) - vmax3r_2 * ELAC/(K_mELAC_2 + ELAC) 
     
    v5_2 = vmax5f_2 * EGLC/(K_mEGLC_2 + EGLC) - vmax5r_2 * EALA/(K_mEALA_2 + EALA)
     
    v7_2 = vmax7_2 * ESER/(K_mESER_2 + ESER)
     
    v16_2 = vmax16f_2 * EGLN/(K_mEGLN_2 + EGLN) * K_iLACtoGLNS_2/(K_iLACtoGLNS_2 + ELAC) - vmax16r_2 * EGLU/(K_mEGLU_2 + EGLU) # * NH4/(K_mNH4_2 + NH4)

    v18_2 = vmax18f_2 * EGLU/(K_mEGLU_2 + EGLU)  - vmax18r_2 * EGLC/(K_mEGLC_2 + EGLC)

    v19_2 = vmax19_2 * EASP/(K_mEASP_2 + EASP) #- vmax19r_2 * EGLU/(K_mEGLU_2 + EGLU) * NH4/(K_mNH4_2 + NH4)
     
    v21_2 = vmax21_2 * ELEU/(K_mELEU_2 + ELEU)
     
    v23_2 = vmax23_2 * EVAL/(K_mEVAL_2 + EVAL)
     
    v25_2 = vmax25_2 * EILE/(K_mEILE_2 + EILE)

    v27_2 = vmax27f_2 * NH4/(K_mNH4_2 + NH4) - vmax27r_2 * ENH4/(K_mENH4_2 + ENH4)

    v28_2 = vmax28_2 * EGLN/(K_mEGLNtoIgG_2 + EGLN) * EGLU/(K_mEGLUtoIgG_2 + EGLU) * EALA/(K_mEALAtoIgG_2 + EALA) * EASP/(K_mEASPtoIgG_2 + EASP) * ESER/(K_mESERtoIgG_2 + ESER) * ELEU/(K_mELEUtoIgG_2 + ELEU) * EILE/(K_mEILEtoIgG_2 + EILE) * EVAL/(K_mEVALtoIgG_2 + EVAL)

    # v28_2 = 0.0004

    v30_2 = vmax30_2 * EGLC/(K_mEGLCtoBIO_2 + EGLC) * EGLN/(K_mEGLNtoBIO_2 + EGLN) * EGLU/(K_mEGLUtoBIO_2 + EGLU) * EALA/(K_mEALAtoBIO_2 + EALA) * EASP/(K_mEASPtoBIO_2 + EASP) * ESER/(K_mESERtoBIO_2 + ESER) * ELEU/(K_mELEUtoBIO_2 + ELEU) * EILE/(K_mEILEtoBIO_2 + EILE) * EVAL/(K_mEVALtoBIO_2 + EVAL)
    
    # v30_2 = 0.008

    v1_2 = v2_2
    
    v4_2 = v3_2
    
    v6_2 = v5_2
    
    v8_2 = v7_2
    
    v14_2 = v18_2 + v19_2 + v23_2 + v25_2
    
    v9_2 = 2 * v2_2 + v7_2 + v14_2 - v3_2 - v5_2
    
    v10_2 = v9_2 + 3 * v21_2 + v25_2
            
    v11_2 = v10_2 + v18_2
    
    v12_2 = v11_2 + v23_2 + v25_2
    
    v13_2 = v12_2 - v14_2
    
    v15_2 = v16_2
    
    v17_2 = v16_2 - v18_2 - v5_2 + v21_2 + v19_2 + v23_2 + v25_2
    
    v20_2 = v19_2
    
    v22_2 = v21_2
    
    v24_2 = v23_2
    
    v26_2 = v25_2
    
    v29_2 = v28_2
    
    # v30_2 = v29_2
    
    v_2= [v1_2, v2_2, v3_2, v4_2, v5_2, v6_2, v7_2, v8_2, v9_2, v10_2, v11_2, v12_2, v13_2, v14_2, v15_2, v16_2, v17_2, v18_2,	
         v19_2, v20_2, v21_2, v22_2, v23_2, v24_2, v25_2, v26_2, v27_2, v28_2, v29_2, v30_2]
    
    
    ##Decline Phase
    v2_3 = vmax2_3 * EGLC/(K_mEGLC_3 + EGLC)  * K_iLACtoHK_3/(K_iLACtoHK_3 + ELAC)
    
    # v3_3 = vmax3f_3 * EGLC/(K_mEGLC_3 * (1 + K_aGLN_3/EGLN) + EGLC) - vmax3r_3 * ELAC/(K_mELAC_3 + ELAC) 

    v3_3 = vmax3f_3 * EGLC/(K_mEGLC_3 + EGLC) - vmax3r_3 * ELAC/(K_mELAC_3 + ELAC)

    v5_3 = vmax5f_3 * EGLC/(K_mEGLC_3 + EGLC) - vmax5r_3 * EALA/(K_mEALA_3 + EALA)

    v7_3 = vmax7_3 * ESER/(K_mESER_3 + ESER)

    v16_3 = vmax16f_3 * EGLN/(K_mEGLN_3 + EGLN) * K_iLACtoGLNS_3/(K_iLACtoGLNS_3 + ELAC) - vmax16r_3 * EGLU/(K_mEGLU_3 + EGLU) # * NH4/(K_mNH4_3 + NH4)

    v18_3 = vmax18f_3 * EGLU/(K_mEGLU_3 + EGLU) - vmax18r_3 * EGLC/(K_mEGLC_3 + EGLC)

    v19_3 = vmax19_3 * EASP/(K_mEASP_3 + EASP) #- vmax19r_3 * EGLU/(K_mEGLU_3 + EGLU) * NH4/(K_mNH4_3 + NH4)
     
    v21_3 = vmax21_3 * ELEU/(K_mELEU_3 + ELEU)
     
    v23_3 = vmax23_3 * EVAL/(K_mEVAL_3 + EVAL)
     
    v25_3 = vmax25_3 * EILE/(K_mEILE_3 + EILE)

    v27_3 = vmax27f_3 * NH4/(K_mNH4_3 + NH4) - vmax27r_3 * ENH4/(K_mENH4_3 + ENH4)
     
    v28_3 = vmax28_3 * EGLN/(K_mEGLNtoIgG_3 + EGLN) * EGLU/(K_mEGLUtoIgG_3 + EGLU) * EALA/(K_mEALAtoIgG_3 + EALA) * EASP/(K_mEASPtoIgG_3 + EASP) * ESER/(K_mESERtoIgG_3 + ESER) * ELEU/(K_mELEUtoIgG_3 + ELEU) * EILE/(K_mEILEtoIgG_3 + EILE) * EVAL/(K_mEVALtoIgG_3 + EVAL)
    
    # v28_3 = 0.0003

    v30_3 = vmax30_3 * EGLC/(K_mEGLCtoBIO_3 + EGLC) * EGLN/(K_mEGLNtoBIO_3 + EGLN) * EGLU/(K_mEGLUtoBIO_3 + EGLU) * EALA/(K_mEALAtoBIO_3 + EALA) * EASP/(K_mEASPtoBIO_3 + EASP) * ESER/(K_mESERtoBIO_3 + ESER) * ELEU/(K_mELEUtoBIO_3 + ELEU) * EILE/(K_mEILEtoBIO_3 + EILE) * EVAL/(K_mEVALtoBIO_3 + EVAL)
    
    # v30_3 = -0.0055

    v1_3 = v2_3
    
    v4_3 = v3_3
    
    v6_3 = v5_3
    
    v8_3 = v7_3
    
    v14_3 = v18_3 + v19_3 + v23_3 + v25_3
    
    v9_3 = 2 * v2_3 + v7_3 + v14_3 - v3_3 - v5_3
    
    v10_3 = v9_3 + 3 * v21_3 + v25_3
            
    v11_3 = v10_3 + v18_3
    
    v12_3 = v11_3 + v23_3 + v25_3
    
    v13_3 = v12_3 - v14_3
    
    v15_3 = v16_3
    
    v17_3 = v16_3 - v18_3 - v5_3 + v21_3 + v19_3 + v23_3 + v25_3
    
    v20_3 = v19_3
    
    v22_3 = v21_3
    
    v24_3 = v23_3
    
    v26_3 = v25_3
    
    v29_3 = v28_3
    
    # v30_3 = v29_3
    
    v_3= [v1_3, v2_3, v3_3, v4_3, v5_3, v6_3, v7_3, v8_3, v9_3, v10_3, v11_3, v12_3, v13_3, v14_3, v15_3, v16_3, v17_3, v18_3,	
         v19_3, v20_3, v21_3, v22_3, v23_3, v24_3, v25_3, v26_3, v27_3, v28_3, v29_3, v30_3]

    du_0 = (N @ v_0 * VCD_0)[:-1]  # Exclude Biomass_0 change for du_0
    du_1 = (N @ v_1 * VCD_1)[:-1]  # Exclude Biomass_1 change for du_1
    du_2 = (N @ v_2 * VCD_2)[:-1]  # Exclude Biomass_2 change for du_2
    du_3 = (N_3 @ v_3 * VCD_3)[:-1]  # Exclude Biomass_3 change for du_3
    dx_0 = v30_0*VCD_0
    dx_1 = v30_1*VCD_1
    dx_2 = v30_2*VCD_2
    dx_3 = v30_3*VCD_3

    # dx = v30_0*VCD_0 + v30_1*VCD_1
    if mode == 0:
        # return v_0, v_1, v_2
        return du_0, du_1, du_2, du_3, dx_0, dx_1, dx_2, dx_3
    if mode == 1:
        # store fluxes to Excel: create if missing, otherwise append rows
        flux_file = 'flux.xlsx'
        # column names: Time, then v1..v30 for each phase suffix _0, _1, _2, _3
        cols = ['Time'] + [f'v{i+1}_{phase}' for phase in (0, 1, 2, 3) for i in range(30)]
        # prepare row values: include current time t then all fluxes
        flux_vals = list(v_0) + list(v_1) + list(v_2) + list(v_3)
        row_vals = [float(t)] + [float(x) for x in flux_vals]
        df_row = pd.DataFrame([row_vals], columns=cols)

        try:
            if not os.path.exists(flux_file):
                df_row.to_excel(flux_file, index=False)
            else:
                df_existing = pd.read_excel(flux_file)
                df_concat = pd.concat([df_existing, df_row], ignore_index=True)
                df_concat.to_excel(flux_file, index=False)
        except Exception:
            # fallback: attempt to write a fresh file if anything goes wrong
            df_row.to_excel(flux_file, index=False)
        return du_0, du_1, du_2, du_3, dx_0, dx_1, dx_2, dx_3
    # return np.hstack([du[:2], np.repeat(0,3), du[-40:], dx]).tolist()



#%%
# def transition_prob_01(t, qO_t, v_ELAC, v_EGLN, v_NH4, beta):
def transition_prob_01(t, batch, beta):
    """
    Calculate the transition probability from state 0 to state 1 at time t.

    Parameters:
    t (float): The time at which the transition probability is calculated.
    beta (dict): A dictionary containing the beta coefficients.

    Returns:
    float: The transition probability from state 0 to state 1 at time t.
    """
    # Calculate p^01[s_t_h]
    p01_t_denominator = 1 + np.exp(-(beta['beta_01_0'] + beta['beta_01_1'] * t #+ 
                                         + beta['beta_01_2'] * interpolate_our(t, batch) 
                                         + beta['beta_01_3'] * interpolate_pH(t,batch)                                         
                                         #+ beta['beta_01_3'] * v_ELAC + 
                                         # beta['beta_01_4'] * v_EGLN + beta['beta_01_5'] * v_NH4
                                         ))
    p01_t = 1 / p01_t_denominator
    
    return p01_t

# def transition_prob_02(t, qO_t, v_ELAC, beta):
def transition_prob_12(t, batch, beta):
    """
    Calculate the transition probability from state 1 to state 2 at time t.

    Parameters:
    t (float): The time at which the transition probability is calculated.
    beta (dict): A dictionary containing the beta values.

    Returns:
    float: The transition probability from state 0 to state 1 at time t.
    """
    # Calculate p^01[s_t_h]
    p12_t_denominator = 1 + np.exp(-(beta['beta_12_0'] + beta['beta_12_1'] * t #+ 
                                         + beta['beta_12_2'] * interpolate_our(t, batch) 
                                         + beta['beta_12_3'] * interpolate_pH(t,batch)
                                         #qO_t + beta['beta_12_3'] * v_ELAC
                                         ))
    p12_t = 1 / p12_t_denominator

    return p12_t


def transition_prob_23(t, batch, beta):
    """
    Calculate the transition probability from state 2 to state 3 at time t.

    Parameters:
    t (float): The time at which the transition probability is calculated.
    beta (dict): A dictionary containing the beta values.

    Returns:
    float: The transition probability from state 2 to state 3 at time t.
    """
    # Calculate p^23[s_t_h]
    p23_t_denominator = 1 + np.exp(-(beta['beta_23_0'] + beta['beta_23_1'] * t #+ 
                                         + beta['beta_23_2'] * interpolate_our(t, batch) 
                                         + beta['beta_23_3'] * interpolate_pH(t,batch)
                                         #qO_t + beta['beta_23_3'] * v_ELAC
                                         ))
    p23_t = 1 / p23_t_denominator

    return p23_t


#%%

def model_prediction_determ(t, x0, ps, mode, batch, t_eval=None, growth_prop=0, metabolite_prop=0):
    """
    Solution with initial condition x(0) = x0
    Args:
        t (array): Time points where the solution is evaluated.
        x0 (array): Initial condition.
        ps (dict): Model parameters.
        mode (int): The mode indicating the mode function to use.
        batch (int): Batch number.
        t_eval (array, optional): Specific time points to evaluate the solution. Defaults to None.
        growth_prop (float, optional): Growth proportion. Defaults to 0.
        metabolite_prop (float, optional): Metabolite proportion. Defaults to 0.
    """
    # def sde(t, xs, ps):
    #    du_0_t, du_1_t, du_2_t, du_3_t, dx_0_t, dx_1_t, dx_2_t, dx_3_t = f_m(t, xs, ps, model) # deterministic part
    #    g_val = 0.01 * np.sqrt(np.abs(f_val)) # stochastic part (scaled by a small factor)
    #    return f_val, g_val

    num_steps = len(t)
    num_vars = len(x0)
    x = np.zeros((num_vars, num_steps))
    x[:, 0] = x0

    for i in range(1, num_steps):
        # if i % 10 == 0:
        #     print(f"Time step {i}/{num_steps-1}, Time: {t[i]:.2f}")
        if t[i] in feeding_cont[batch]['Day'].values: # Check if feeding occurs at time i
            # Post-feeding adjustement
            x[target_label_mapper['EGLC']:target_label_mapper['BIOM_0'], i] = x[target_label_mapper['EGLC']:target_label_mapper['BIOM_0'], i-1] + feeding_cont[batch][feeding_cont[batch]['Day'] == t[i]].values[0][1:-1]  # Add feed concentration at feeding times
            x[target_label_mapper['BIOM_0']:, i] = x[target_label_mapper['BIOM_0']:, i-1] # Biomass and VCD remain unchanged during feeding

        else:
            du_0_t, du_1_t, du_2_t, du_3_t, dx_0_t, dx_1_t, dx_2_t, dx_3_t = f_m(t[i], x[:, i-1], ps, mode, batch)  # sde(t[i-1], x[:, i-1], ps)
            # prevent large negative changes when concentrations are (near) zero
            small_mask = x[:target_label_mapper['BIOM_0'], i-1] < 1e-3
            if np.any(small_mask):
                for du in (du_0_t, du_1_t, du_2_t, du_3_t):
                    # only scale negative rates for metabolites that are already very small
                    neg_mask = small_mask & (du < 0)
                    if np.any(neg_mask):
                        du[neg_mask] *= 0.001

            # dw = np.random.normal(0, np.sqrt(Delta_t), size=num_vars)
            x[:target_label_mapper['BIOM_0'],i] = x[:target_label_mapper['BIOM_0'],i-1] + (du_0_t + du_1_t + du_2_t + du_3_t) * 24 * Delta_t # Metabolite

            ###Phase transition probabilities
            # p01_t = transition_prob_01(t, qO_t, v_ELAC_0, v_EGLN_0, v_NH4_0, beta_01)
            p01_t = transition_prob_01(t[i], batch, ps)
            # p01_t = 0.0
            p00_t = 1 - p01_t
            # p12_t = transition_prob_02(t, qO_t, v_ELAC_1, beta)
            p12_t = transition_prob_12(t[i], batch, ps)
            # p12_t = 0.0
            p11_t = 1 - p12_t

            # p23_t = transition_prob_23(t, qO_t, v_ELAC_2, beta)
            p23_t = transition_prob_23(t[i], batch, ps)
            p22_t = 1 - p23_t

            P_t = np.array([[p00_t, p01_t, 0, 0],
                        [0, p11_t, p12_t, 0],
                        [0, 0, p22_t, p23_t],
                        [0, 0, 0, 1]])
        

            x[target_label_mapper['VCD_0'],i] = x[target_label_mapper['VCD_0'],i-1] + dx_0_t * 24 * Delta_t # VCD_0
            x[target_label_mapper['VCD_1'],i] = x[target_label_mapper['VCD_1'],i-1] + dx_1_t * 24 * Delta_t # VCD_1
            x[target_label_mapper['VCD_2'],i] = x[target_label_mapper['VCD_2'],i-1] + dx_2_t * 24 * Delta_t # VCD_2
            x[target_label_mapper['VCD_3'],i] = x[target_label_mapper['VCD_3'],i-1] + dx_3_t * 24 * Delta_t # VCD_3

            x[target_label_mapper['VCD_0']:target_label_mapper['VCD_3']+1,i] = x[target_label_mapper['VCD_0']:target_label_mapper['VCD_3']+1,i] @ P_t  # Apply phase transition probabilities

            x[target_label_mapper['BIOM_0'],i] = x[target_label_mapper['VCD_0'],i] * biomass
            x[target_label_mapper['BIOM_1'],i] = x[target_label_mapper['VCD_1'],i] * biomass
            x[target_label_mapper['BIOM_2'],i] = x[target_label_mapper['VCD_2'],i] * biomass
            x[target_label_mapper['BIOM_3'],i] = x[target_label_mapper['VCD_3'],i] * biomass
        


    return x


#%%
def model_prediction_stochastic(t, x0, ps, mode, batch, sigma, t_eval=None, growth_prop=0, metabolite_prop=0):
    """
    Solution with initial condition x(0) = x0
    Args:
        t (array): Time points where the solution is evaluated.
        x0 (array): Initial condition.
        ps (dict): Model parameters.
        mode (int): The mode indicating the mode function to use.
        batch (int): Batch number.
        t_eval (array, optional): Specific time points to evaluate the solution. Defaults to None.
        growth_prop (float, optional): Growth proportion. Defaults to 0.
        metabolite_prop (float, optional): Metabolite proportion. Defaults to 0.
    """
    # def sde(t, xs, ps):
    #    du_0_t, du_1_t, du_2_t, du_3_t, dx_0_t, dx_1_t, dx_2_t, dx_3_t = f_m(t, xs, ps, model) # deterministic part
    #    g_val = 0.01 * np.sqrt(np.abs(f_val)) # stochastic part (scaled by a small factor)
    #    return f_val, g_val

    num_steps = len(t)
    num_vars = len(x0)
    x = np.zeros((num_vars, num_steps))
    x[:, 0] = x0

    for i in range(1, num_steps):
        # if i % 10 == 0:
        #     print(f"Time step {i}/{num_steps-1}, Time: {t[i]:.2f}")
        if t[i] in feeding_cont[batch]['Day'].values: # Check if feeding occurs at time i
            # Post-feeding adjustement
            x[target_label_mapper['EGLC']:target_label_mapper['BIOM_0'], i] = x[target_label_mapper['EGLC']:target_label_mapper['BIOM_0'], i-1] + feeding_cont[batch][feeding_cont[batch]['Day'] == t[i]].values[0][1:-1]  # Add feed concentration at feeding times
            x[target_label_mapper['BIOM_0']:, i] = x[target_label_mapper['BIOM_0']:, i-1] # Biomass and VCD remain unchanged during feeding

        else:
            du_0_t, du_1_t, du_2_t, du_3_t, dx_0_t, dx_1_t, dx_2_t, dx_3_t = f_m(t[i], x[:, i-1], ps, mode, batch)  # sde(t[i-1], x[:, i-1], ps)
            # generate noise per element (each du_* is an array, dx_* are scalars) and keep as a Python list
            loc_list = [du_0_t, du_1_t, du_2_t, du_3_t, dx_0_t, dx_1_t, dx_2_t, dx_3_t]            
            d_stoc_t = [np.random.normal(loc=val, scale=sigma[i, :]*np.sqrt(np.abs(val))) for i, val in enumerate(loc_list)]  # Stochastic part (list of arrays/scalars)
            
            # prevent large negative changes when concentrations are (near) zero
            small_mask = x[:target_label_mapper['BIOM_0'], i-1] < 1e-3
            if np.any(small_mask):
                for du in (d_stoc_t[0], d_stoc_t[1], d_stoc_t[2], d_stoc_t[3]):
                    # only scale negative rates for metabolites that are already very small
                    neg_mask = small_mask & (du < 0)
                    if np.any(neg_mask):
                        du[neg_mask] *= 0.01

            x[:target_label_mapper['BIOM_0'],i] = x[:target_label_mapper['BIOM_0'],i-1] + (d_stoc_t[0] + d_stoc_t[1] + d_stoc_t[2] + d_stoc_t[3]) * 24 * Delta_t # Metabolite

            ###Phase transition probabilities
            # p01_t = transition_prob_01(t, qO_t, v_ELAC_0, v_EGLN_0, v_NH4_0, beta_01)
            p01_t = transition_prob_01(t[i], batch, ps)
            # p01_t = 0.0
            p00_t = 1 - p01_t
            # p12_t = transition_prob_02(t, qO_t, v_ELAC_1, beta)
            p12_t = transition_prob_12(t[i], batch, ps)
            # p12_t = 0.0
            p11_t = 1 - p12_t

            # p23_t = transition_prob_23(t, qO_t, v_ELAC_2, beta)
            p23_t = transition_prob_23(t[i], batch, ps)
            p22_t = 1 - p23_t

            P_t = np.array([[p00_t, p01_t, 0, 0],
                        [0, p11_t, p12_t, 0],
                        [0, 0, p22_t, p23_t],
                        [0, 0, 0, 1]])
        

            x[target_label_mapper['VCD_0'],i] = x[target_label_mapper['VCD_0'],i-1] + d_stoc_t[4][0] * 24 * Delta_t # VCD_0
            x[target_label_mapper['VCD_1'],i] = x[target_label_mapper['VCD_1'],i-1] + d_stoc_t[5][0] * 24 * Delta_t # VCD_1
            x[target_label_mapper['VCD_2'],i] = x[target_label_mapper['VCD_2'],i-1] + d_stoc_t[6][0] * 24 * Delta_t # VCD_2
            x[target_label_mapper['VCD_3'],i] = x[target_label_mapper['VCD_3'],i-1] + d_stoc_t[7][0] * 24 * Delta_t # VCD_3

            x[target_label_mapper['VCD_0']:target_label_mapper['VCD_3']+1,i] = x[target_label_mapper['VCD_0']:target_label_mapper['VCD_3']+1,i] @ P_t  # Apply phase transition probabilities

            x[target_label_mapper['BIOM_0'],i] = x[target_label_mapper['VCD_0'],i] * biomass
            x[target_label_mapper['BIOM_1'],i] = x[target_label_mapper['VCD_1'],i] * biomass
            x[target_label_mapper['BIOM_2'],i] = x[target_label_mapper['VCD_2'],i] * biomass
            x[target_label_mapper['BIOM_3'],i] = x[target_label_mapper['VCD_3'],i] * biomass
        


    return x


#%%

def calc_wape_for_batch(batch, params, mode_var, determ=True):
    """
    Calculate WAPE per extracellular key and overall WAPE for a given batch.
    Uses globals: Delta_t, data_extra_cont, feeding_cont, target_label_mapper, biomass, model_prediction_determ.
    Returns: (wape_dict, overall_wape, t, x) where x is the prediction matrix (states x times).
    WAPE is returned as percentage (0-100). If sum(|observed|)==0, WAPE for that key is np.nan.
    """
    # construct time vector including feeding and measurement days
    last_day = data_extra_cont[batch]['Day'].iloc[-1]
    t = np.arange(0, last_day + Delta_t, Delta_t)
    t = np.unique(np.concatenate((t, feeding_cont[batch]['Day'].values, data_extra_cont[batch]['Day'].values)))

    # initial condition (same construction used elsewhere)
    x0 = np.concatenate((np.zeros(target_label_mapper['ANTI'] + 1), data_extra_cont[batch].iloc[0, 1:-1].values))
    x0 = np.hstack((x0, np.array((data_extra_cont[batch].iloc[0, -1] * biomass, 0, 0, 0, data_extra_cont[batch].iloc[0, -1], 0, 0, 0))))

    # get predictions (mode=1 to return concentration derivatives as used)
    # remove existing flux file so f_m will create a fresh one
    flux_path = 'flux.xlsx'
    if os.path.exists(flux_path):
        try:
            os.remove(flux_path)
        except OSError:
            pass
    if determ:
        x = model_prediction_determ(t, x0, params, mode=mode_var, batch=batch, t_eval=None)
    else:
        x = model_prediction_stochastic(t, x0, params, mode=mode_var, batch=batch, t_eval=None)

    extracellular_keys = ['EGLC','EGLN','EGLU','EASP','EALA','ELAC',
                          'ESER','ELEU','EILE','EVAL','ENH4','EANTI','VCD']

    wape_dict = {}
    total_abs = 0.0
    total_obs_abs = 0.0

    meas_days = data_extra_cont[batch]['Day'].values

    for key in extracellular_keys:
        if key == 'VCD':
            pred_ts = (x[target_label_mapper['VCD_0'], :] +
                       x[target_label_mapper['VCD_1'], :] +
                       x[target_label_mapper['VCD_2'], :] +
                       x[target_label_mapper['VCD_3'], :])
            meas_vals = data_extra_cont[batch]['VCD_0'].values
        else:
            idx = target_label_mapper[key]
            pred_ts = x[idx, :]
            meas_vals = data_extra_cont[batch][key].values

        # sample predictions at measurement days
        pred_at_meas = np.interp(meas_days, t, pred_ts)

        valid = ~np.isnan(meas_vals)
        if np.any(valid):
            abs_errs = np.abs(pred_at_meas[valid] - meas_vals[valid])
            obs_abs = np.abs(meas_vals[valid])
            numerator = float(np.sum(abs_errs))
            denominator = float(np.sum(obs_abs))
            if denominator > 0:
                wape = 100.0 * numerator / denominator
            else:
                wape = np.nan
            wape_dict[key] = wape
            total_abs += numerator
            total_obs_abs += denominator
        else:
            wape_dict[key] = np.nan

    overall_wape = 100.0 * total_abs / total_obs_abs if total_obs_abs > 0 else np.nan

    return wape_dict, overall_wape, t, x



#%%
###set parameters incluing bounds
theta = Parameters()

##Early Growth Phase
theta.add('vmax2_0', value=1.10720253209616235, min=0, max=4*1E2)
# theta.add('vmax2_0', value=0.20720253209616235, min=0, max=4*1E2)
theta.add('vmax3f_0', value=1.10073394869049, min=0, max=2*1E2)
# theta.add('vmax3f_0', value=0.2420073394869049, min=0, max=2*1E2)
theta.add('vmax3r_0', value=0.07218660357981577, min=0, max=1*1E2)
theta.add('vmax5f_0', value=0.24, min=0, max=1*1E2)
# theta.add('vmax5f_0', value=0.05, min=0, max=1*1E2)
theta.add('vmax5r_0', value=0.0005, min=0, max=0.4*1E2)
theta.add('vmax7_0', value=0.0262731193747201, min=0, max=1*1E2)
theta.add('vmax16f_0', value=0.17, min=0, max=1*1E2)
theta.add('vmax16r_0', value=0.003, min=0, max=0.1*1E2)
theta.add('vmax18f_0', value=0.000, min=0, max=1*1E2)
theta.add('vmax18r_0', value=0.013, min=0, max=0.1*1E2)
theta.add('vmax19_0', value=6.601569291220244e-10, min=0, max=0.1*1E2)
theta.add('vmax21_0', value=0.00, min=0, max=1*1E2)
theta.add('vmax23_0', value=0.00, min=0, max=1*1E2)
theta.add('vmax25_0', value=0.00, min=0, max=1*1E2)
theta.add('vmax27f_0', value=0.02, min=0, max=0.2*1E2)
theta.add('vmax27r_0', value=0.01, min=0, max=0.1*1E2)
theta.add('vmax28_0', value=0.0008, min=0, max=1*1E2)
theta.add('vmax30_0', value=0.032358314716596226, min=0, max=1*1E5)

theta.add('K_mEGLC_0', value=24.9998285002715, min=0, max=0.5*1E2)
theta.add('K_iLACtoHK_0', value=600, min=0, max=3*1E3)
theta.add('K_iLACtoGLNS_0', value=500, min=0, max=3*1E3)
theta.add('K_aGLN_0', value=2, min=0, max=0.3*1E2)
theta.add('K_mELAC_0', value=6.202813724972499, min=0, max=0.3*1E2)
theta.add('K_mEALA_0', value=1.3, min=0, max=0.04*1E2)
theta.add('K_mESER_0', value=7.316149453202534, min=0, max=0.15*1E2)
theta.add('K_mEGLN_0', value=14, min=0, max=0.3*1E2)
theta.add('K_mNH4_0', value=4, min=0, max=0.2*1E2)
theta.add('K_mENH4_0', value=0.5, min=0, max=0.1*1E2)
theta.add('K_mEGLU_0', value=2, min=0, max=0.3*1E2)
theta.add('K_mEASP_0', value=1.8520762392081802, min=0, max=0.2*1E2)
theta.add('K_mELEU_0', value=6, min=0, max=0.2*1E2)
theta.add('K_mEILE_0', value=3, min=0, max=0.2*1E2)
theta.add('K_mEVAL_0', value=4, min=0, max=0.2*1E2)

theta.add('K_mEGLNtoIgG_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoIgG_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEASPtoIgG_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoIgG_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoIgG_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoIgG_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoIgG_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoIgG_0', value=1E-3, min=1E-5, max=1E2)

theta.add('K_mEGLCtoBIO_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLNtoBIO_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoBIO_0', value=1E-3, min=1E-5, max=1E2)  
theta.add('K_mEASPtoBIO_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoBIO_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoBIO_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoBIO_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoBIO_0', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoBIO_0', value=1E-3, min=1E-5, max=1E2)

theta.add('k_0', value=2.798, min=0, max=1*1E2)
theta.add('K1_0', value=1e-8, min=0, max=1*1E2)
theta.add('K2_0', value=1e-8, min=0, max=1*1E2)

##Late Exponential Phase
theta.add('vmax2_1', value=0.08701946222030127, min=0, max=4*1E2)
theta.add('vmax3f_1', value=0.033156207163731466, min=0, max=2*1E2)
theta.add('vmax3r_1', value=0.06836652381304598, min=0, max=1*1E2)
theta.add('vmax5f_1', value=0.01, min=0, max=1*1E2)
theta.add('vmax5r_1', value=0.002, min=0, max=0.4*1E2)
theta.add('vmax7_1', value=0.002844105585730894, min=0, max=1*1E2)
theta.add('vmax16f_1', value=0.004, min=0, max=1*1E2)
theta.add('vmax16r_1', value=0.012, min=0, max=0.1*1E2)
theta.add('vmax18f_1', value=0.000, min=0, max=1*1E2)
theta.add('vmax18r_1', value=0.002, min=0, max=0.1*1E2)
theta.add('vmax19_1', value=6.640286653869509e-09, min=0, max=0.1*1E2)
theta.add('vmax21_1', value=0.007, min=0, max=1*1E2)
theta.add('vmax23_1', value=0.0045, min=0, max=1*1E2)
theta.add('vmax25_1', value=0.002, min=0, max=1*1E2)
theta.add('vmax27f_1', value=0.02, min=0, max=0.2*1E2)
theta.add('vmax27r_1', value=0.01, min=0, max=0.1*1E2)
theta.add('vmax28_1', value=0.00042, min=0, max=0.1*1E2)
theta.add('vmax30_1', value=0.02373936522464892, min=0, max=1*1E5)

theta.add('K_mEGLC_1', value=25.00084433353726, min=0, max=0.5*1E2)
theta.add('K_iLACtoHK_1', value=600, min=0, max=3*1E3)
theta.add('K_iLACtoGLNS_1', value=500, min=0, max=3*1E3)
theta.add('K_aGLN_1', value=2, min=0, max=0.3*1E2)
theta.add('K_mELAC_1', value=5.365332099429065, min=0, max=0.3*1E2)
theta.add('K_mEALA_1', value=1.3, min=0, max=0.04*1E2)
theta.add('K_mESER_1', value=1.8488330391242271, min=0, max=0.15*1E2)
theta.add('K_mEGLN_1', value=14, min=0, max=0.3*1E2)
theta.add('K_mNH4_1', value=4, min=0, max=0.2*1E2)
theta.add('K_mENH4_1', value=0.5, min=0, max=0.1*1E2)
theta.add('K_mEGLU_1', value=2, min=0, max=0.3*1E2)
theta.add('K_mEASP_1', value=1.7493272823215666, min=0, max=0.2*1E2)
theta.add('K_mELEU_1', value=6, min=0, max=0.2*1E2)
theta.add('K_mEILE_1', value=3, min=0, max=0.2*1E2)
theta.add('K_mEVAL_1', value=4, min=0, max=0.2*1E2)

theta.add('K_mEGLNtoIgG_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoIgG_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEASPtoIgG_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoIgG_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoIgG_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoIgG_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoIgG_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoIgG_1', value=1E-3, min=1E-5, max=1E2)

theta.add('K_mEGLCtoBIO_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLNtoBIO_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoBIO_1', value=1E-3, min=1E-5, max=1E2)  
theta.add('K_mEASPtoBIO_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoBIO_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoBIO_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoBIO_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoBIO_1', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoBIO_1', value=1E-3, min=1E-5, max=1E2)

##Stational Phase
theta.add('vmax2_2', value=0.03201946222030127, min=0, max=4*1E2)
theta.add('vmax3f_2', value=0.033156207163731466, min=0, max=2*1E2)
theta.add('vmax3r_2', value=0.06836652381304598, min=0, max=1*1E2)
theta.add('vmax5f_2', value=0.0043, min=0, max=1*1E2)
theta.add('vmax5r_2', value=0.00, min=0, max=0.4*1E2)
theta.add('vmax7_2', value=2.6313441980896357e-07, min=0, max=1*1E2)
theta.add('vmax16f_2', value=0.0, min=0, max=1*1E2)
theta.add('vmax16r_2', value=0.0063, min=0, max=0.1*1E2)
theta.add('vmax18f_2', value=0.002, min=0, max=1*1E2)
theta.add('vmax18r_2', value=0.001, min=0, max=0.1*1E2)
theta.add('vmax19_2', value=0.003109839922081024, min=0, max=0.1*1E2)
theta.add('vmax21_2', value=0.000, min=0, max=1*1E2)
theta.add('vmax23_2', value=0.000, min=0, max=1*1E2)
theta.add('vmax25_2', value=0.000, min=0, max=1*1E2)
theta.add('vmax27f_2', value=0.02, min=0, max=0.2*1E2)
theta.add('vmax27r_2', value=0.01, min=0, max=0.1*1E2)
theta.add('vmax28_2', value=0.0004, min=0, max=0.1*1E2)
theta.add('vmax30_2', value=0.011665127869209257, min=0, max=1*1E5)

theta.add('K_mEGLC_2', value=25.00084433353726, min=0, max=0.5*1E2)
theta.add('K_iLACtoHK_2', value=600, min=0, max=3*1E3)
theta.add('K_iLACtoGLNS_2', value=500, min=0, max=3*1E3)
theta.add('K_aGLN_2', value=2, min=0, max=0.3*1E2)
theta.add('K_mELAC_2', value=5.365332099429065, min=0, max=0.3*1E2)
theta.add('K_mEALA_2', value=1.3, min=0, max=0.04*1E2)
theta.add('K_mESER_2', value=14.892890909804066, min=0, max=0.15*1E2)
theta.add('K_mEGLN_2', value=14, min=0, max=0.3*1E2)
theta.add('K_mNH4_2', value=4, min=0, max=0.2*1E2)
theta.add('K_mENH4_2', value=0.5, min=0, max=0.1*1E2)
theta.add('K_mEGLU_2', value=2, min=0, max=0.3*1E2)
theta.add('K_mEASP_2', value=2.2489535724871823, min=0, max=0.2*1E2)
theta.add('K_mELEU_2', value=6, min=0, max=0.2*1E2)
theta.add('K_mEILE_2', value=3, min=0, max=0.2*1E2)
theta.add('K_mEVAL_2', value=4, min=0, max=0.2*1E2)

theta.add('K_mEGLNtoIgG_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoIgG_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEASPtoIgG_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoIgG_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoIgG_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoIgG_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoIgG_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoIgG_2', value=1E-3, min=1E-5, max=1E2)

theta.add('K_mEGLCtoBIO_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLNtoBIO_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoBIO_2', value=1E-3, min=1E-5, max=1E2)  
theta.add('K_mEASPtoBIO_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoBIO_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoBIO_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoBIO_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoBIO_2', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoBIO_2', value=1E-3, min=1E-5, max=1E2)


##Decline Phase
theta.add('vmax2_3', value=0.0400705865745511, min=0, max=4*1E2)
theta.add('vmax3f_3', value=0.0006747313147894893, min=0, max=2*1E2)
theta.add('vmax3r_3', value=0.0002637912886849225, min=0, max=1*1E2)
theta.add('vmax5f_3', value=0.0000, min=0, max=1*1E2)
theta.add('vmax5r_3', value=0.002, min=0, max=0.4*1E2)
theta.add('vmax7_3', value=0.011160024868961393, min=0, max=1*1E2)
theta.add('vmax16f_3', value=0.0, min=0, max=1*1E2)
theta.add('vmax16r_3', value=0.0018, min=0, max=0.1*1E2)
theta.add('vmax18f_3', value=0.0088, min=0, max=1*1E2)
theta.add('vmax18r_3', value=0.0045, min=0, max=0.1*1E2)
theta.add('vmax19_3', value=0.003262152166469912, min=0, max=0.1*1E2)
theta.add('vmax21_3', value=0.0038, min=0, max=1*1E2)
theta.add('vmax23_3', value=0.002, min=0, max=1*1E2)
theta.add('vmax25_3', value=0.0028, min=0, max=1*1E2)
# theta.add('vmax21_3', value=0.0038, min=0, max=1*1E2)
# theta.add('vmax23_3', value=0.0032, min=0, max=1*1E2)
# theta.add('vmax25_3', value=0.0024, min=0, max=1*1E2)
theta.add('vmax27f_3', value=0.02, min=0, max=0.2*1E2)
theta.add('vmax27r_3', value=0.01, min=0, max=0.1*1E2)
theta.add('vmax28_3', value=0.0001, min=0, max=0.1*1E2)
theta.add('vmax30_3', value=-0.0055997028483249011, min=-0.04 * 600, max=0.04 * 600)

theta.add('K_mEGLC_3', value=10.3678240759968, min=0, max=0.5*1E2)
theta.add('K_iLACtoHK_3', value=600, min=0, max=3*1E3)
theta.add('K_iLACtoGLNS_3', value=500, min=0, max=3*1E3)
theta.add('K_aGLN_3', value=2, min=0, max=0.3*1E2)
theta.add('K_mELAC_3', value=4.956932922907262, min=0, max=0.3*1E2)
theta.add('K_mEALA_3', value=1.3, min=0, max=0.04*1E2)
theta.add('K_mESER_3', value=14.073645389412535, min=0, max=0.15*1E2)
theta.add('K_mEGLN_3', value=14, min=0, max=0.3*1E2)
theta.add('K_mNH4_3', value=4, min=0, max=0.2*1E2)
theta.add('K_mENH4_3', value=0.5, min=0, max=0.1*1E2)
theta.add('K_mEGLU_3', value=2, min=0, max=0.3*1E2)
theta.add('K_mEASP_3', value=2.064585882053671, min=0, max=0.2*1E2)
theta.add('K_mELEU_3', value=4, min=0, max=0.2*1E2)
theta.add('K_mEILE_3', value=3, min=0, max=0.2*1E2)
theta.add('K_mEVAL_3', value=4, min=0, max=0.2*1E2)

theta.add('K_mEGLNtoIgG_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoIgG_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEASPtoIgG_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoIgG_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoIgG_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoIgG_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoIgG_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoIgG_3', value=1E-3, min=1E-5, max=1E2)

theta.add('K_mEGLCtoBIO_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLNtoBIO_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEGLUtoBIO_3', value=1E-3, min=1E-5, max=1E2)  
theta.add('K_mEASPtoBIO_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEALAtoBIO_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mESERtoBIO_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mELEUtoBIO_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEILEtoBIO_3', value=1E-3, min=1E-5, max=1E2)
theta.add('K_mEVALtoBIO_3', value=1E-3, min=1E-5, max=1E2)

##Transition Parameters
theta.add('beta_01_0', value=-11, min=-100, max=100)
theta.add('beta_01_1', value=3, min=-100, max=100)
theta.add('beta_01_2', value=-0.11560364182264493, min=-100, max=100)
theta.add('beta_01_3', value=-0.07251129594865802, min=-100, max=100)
theta.add('beta_12_0', value=-20, min=-100, max=100) 
theta.add('beta_12_1', value=3, min=-100, max=100)
theta.add('beta_12_2', value=-0.14189451903004624, min=-100, max=100)
theta.add('beta_12_3', value=-0.07249531621353356, min=-100, max=100)
theta.add('beta_23_0', value=-26, min=-100, max=100)
theta.add('beta_23_1', value=3, min=-100, max=100)
theta.add('beta_23_2', value=-0.13162093443973788, min=-100, max=100)
theta.add('beta_23_3', value=-0.07056700786166914, min=-100, max=100)


#%% Deterministic model prediction with the optimized parameters
for b in ['tigr20_6']:
    # use positional indexing to get the last Day value and np.arange for float step sizes
    last_day = data_extra_cont[b]['Day'].iloc[-1]
    t = np.arange(0, last_day + Delta_t, Delta_t)
    t = np.unique(np.concatenate((t, feeding_cont[b]['Day'].values))) # Ensure feeding days are included
    # print(t)

    x0 = np.concatenate((np.zeros(target_label_mapper['ANTI'] + 1), data_extra_cont[b].iloc[0, 1:-1].values))
    x0 = np.hstack((x0, np.array((data_extra_cont[b].iloc[0, -1] * biomass, 0, 0, 0, data_extra_cont[b].iloc[0, -1], 0, 0, 0))))  # Initial VCD and Biomass
    x = model_prediction_determ(t, x0, theta, mode = 0, batch=b, t_eval=None, growth_prop=0, metabolite_prop=0)


    # Plot predicted vs measured for key extracellular metabolites for batch b

    extracellular_keys = ['EGLC','EGLN','EGLU','EASP','EALA','ELAC',
                          'ESER','ELEU','EILE','EVAL','ENH4','EANTI','VCD', 'VCD_0', 'VCD_1', 'VCD_2', 'VCD_3']

    n = len(extracellular_keys)
    ncols = 4
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*2.8), squeeze=False)
    axes_flat = axes.flatten()

    # loop over extracellular keys and plot each on its own subplot
    for i, key in enumerate(extracellular_keys):
        # select the subplot axes for this metabolite/VCD
        ax = axes_flat[i]

        if key == 'VCD':
            # For VCD we show the total cell density as the sum of all phase-specific VCDs
            vcd_total = (
                x[target_label_mapper['VCD_0'], :]
                + x[target_label_mapper['VCD_1'], :]
                + x[target_label_mapper['VCD_2'], :]
                + x[target_label_mapper['VCD_3'], :]
            )
            # plot predicted total VCD time-series
            ax.plot(t, vcd_total, lw=1.5, color='tab:blue', label='Predicted')
            # measured VCD is stored in the VCD_0 column of the measurements
            meas_vals = data_extra_cont[b]['VCD_0'].values
        elif key in ['VCD_0', 'VCD_1', 'VCD_2', 'VCD_3']:
            # plot predicted phase-specific VCD time-series
            idx = target_label_mapper[key]
            ax.plot(t, x[idx, :], lw=1.5, color='tab:blue', label='Predicted')
            # measured VCD phase is stored in the corresponding VCD_x column of the measurements
            idx = target_label_mapper[key]
            ax.plot(t, x[idx, :], lw=1.5, color='tab:blue', label='Predicted')
            # No phase-specific VCD measurements available — use NaNs so nothing is scattered
            meas_vals = np.full(data_extra_cont[b]['Day'].shape, np.nan, dtype=float)
        else:
            # map metabolite key to row index in prediction matrix and plot predicted values
            idx = target_label_mapper[key]
            ax.plot(t, x[idx, :], lw=1.5, color='tab:blue', label='Predicted')
            # get measured values for this metabolite from the batch data
            meas_vals = data_extra_cont[b][key].values

        # measurement time points (days) for scatter overlay
        meas_days = data_extra_cont[b]['Day'].values

        # filter out NaNs in measurements before plotting
        valid = ~np.isnan(meas_vals)
        if np.any(valid):
            # overlay measured points as orange dots with black edge
            ax.scatter(meas_days[valid], meas_vals[valid], color='tab:orange', edgecolor='k', zorder=5, label='Measured')

        # formatting for the subplot
        ax.set_title(key, fontsize=10)
        ax.set_xlabel('Time (day)', fontsize=9)
        ax.set_ylabel('Concentration', fontsize=9)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.legend(loc='best', fontsize=8)

    # hide any unused subplots in the grid (if n < nrows*ncols)
    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.tight_layout()
    plt.show()


#%% Stochastic model prediction repeated 30 times for one batch
ensemble_size = 30


# reproducible random seeds for numpy and python's random
np.random.seed(5)
# 5 for Rep 1 & 2 & 3
# Load sigma matrix
# load sigma as a numeric numpy matrix (first sheet)
sigma_df = pd.read_excel("sigma.xlsx", sheet_name=0, header=None)
sigma = sigma_df.astype(float).values


for b in ['tigr20_8']:
    last_day = data_extra_cont[b]['Day'].iloc[-1]
    t = np.arange(0, last_day + Delta_t, Delta_t)
    t = np.unique(np.concatenate((t, feeding_cont[b]['Day'].values)))  # Ensure feeding days are included

    x0 = np.concatenate((np.zeros(target_label_mapper['ANTI'] + 1), data_extra_cont[b].iloc[0, 1:-1].values))
    x0 = np.hstack((x0, np.array((data_extra_cont[b].iloc[0, -1] * biomass, 0, 0, 0, data_extra_cont[b].iloc[0, -1], 0, 0, 0))))  # Initial VCD and Biomass

    num_vars = len(x0)
    num_steps = len(t)

    # run ensemble
    ensemble = np.zeros((ensemble_size, num_vars, num_steps))
    for rep in range(ensemble_size):
        ensemble[rep] = model_prediction_stochastic(t, x0, theta, mode=0, batch=b, sigma=sigma, t_eval=None, growth_prop=0, metabolite_prop=0)

    # compute median and 95% PI
    median = np.median(ensemble, axis=0)
    lower = np.percentile(ensemble, 2.5, axis=0)
    upper = np.percentile(ensemble, 97.5, axis=0)

    extracellular_keys = ['EGLC','EGLN','EGLU','EASP','EALA','ELAC',
                          'ESER','ELEU','EILE','EVAL','ENH4','EANTI','VCD', 'VCD_0', 'VCD_1', 'VCD_2', 'VCD_3']

    n = len(extracellular_keys)
    ncols = 4
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*2.8), squeeze=False)
    axes_flat = axes.flatten()

    for i, key in enumerate(extracellular_keys):
        ax = axes_flat[i]

        meas_days = data_extra_cont[b]['Day'].values

        if key == 'VCD':
            # ensemble total VCD per replicate
            vcd_ens = ensemble[:, target_label_mapper['VCD_0'], :] + ensemble[:, target_label_mapper['VCD_1'], :] + ensemble[:, target_label_mapper['VCD_2'], :] + ensemble[:, target_label_mapper['VCD_3'], :]
            med = np.median(vcd_ens, axis=0)
            lo = np.percentile(vcd_ens, 2.5, axis=0)
            hi = np.percentile(vcd_ens, 97.5, axis=0)
            # plot individual trajectories lightly
            for rep in range(ensemble_size):
                ax.plot(t, vcd_ens[rep], color='tab:blue', alpha=0.12, lw=0.8)
            ax.plot(t, med, color='tab:blue', lw=1.5, label='Median Predicted')
            ax.fill_between(t, lo, hi, color='tab:blue', alpha=0.18, label='95% PI')
            meas_vals = data_extra_cont[b]['VCD_0'].values
        elif key in ['VCD_0', 'VCD_1', 'VCD_2', 'VCD_3']:
            idx = target_label_mapper[key]
            # per-replicate trajectories
            for rep in range(ensemble_size):
                ax.plot(t, ensemble[rep, idx, :], color='tab:blue', alpha=0.12, lw=0.8)
            ax.plot(t, median[idx], color='tab:blue', lw=1.5, label='Median Predicted')
            ax.fill_between(t, lower[idx], upper[idx], color='tab:blue', alpha=0.18, label='95% PI')
            meas_vals = np.full(data_extra_cont[b]['Day'].shape, np.nan, dtype=float)
        else:
            idx = target_label_mapper[key]
            # per-replicate trajectories
            for rep in range(ensemble_size):
                ax.plot(t, ensemble[rep, idx, :], color='tab:blue', alpha=0.12, lw=0.8)
            ax.plot(t, median[idx], color='tab:blue', lw=1.5, label='Median Predicted')
            ax.fill_between(t, lower[idx], upper[idx], color='tab:blue', alpha=0.18, label='95% PI')
            meas_vals = data_extra_cont[b][key].values

        # overlay measured points
        valid = ~np.isnan(meas_vals)
        if np.any(valid):
            ax.scatter(meas_days[valid], meas_vals[valid], color='tab:orange', edgecolor='k', zorder=5, label='Measured')

        ax.set_title(key, fontsize=10)
        ax.set_xlabel('Time (day)', fontsize=9)
        ax.set_ylabel('Concentration', fontsize=9)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.legend(loc='best', fontsize=8)

    # hide any unused subplots
    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.tight_layout()
    plt.show()


#%% Save the parameter set to a file 
with open('params_B.pkl', 'wb') as f:
    pickle.dump(theta, f)


#%% Draw the plot for VCD, EGLC, ELAC, EGLN, EGLU, IgG, LEU, ILE, VAL. 
extracellular_keys = {'VCD': 'VCD', 'EGLC': 'Glucose','ELAC': 'Lactate', 'EGLN': 'Glutamine','EGLU': 'Glutamate', 'EANTI': 'IgG', 'ELEU': 'Leucine','EILE': 'Isoleucine','EVAL': 'Valine'}

n = len(extracellular_keys)
ncols = 3
nrows = math.ceil(n / ncols)

fig, axs = plt.subplots(int(nrows), int(ncols), figsize=(15, 3 * nrows))
axs = axs.flatten()

# define tick days to use across subplots
all_days = np.int32(np.unique(data_extra_cont[b]['Day'].values))

for i, key in enumerate(extracellular_keys.keys()):
    ax = axs[i]

    meas_days = data_extra_cont[b]['Day'].values

    if key == 'VCD':
        # ensemble total VCD per replicate
        vcd_ens = ensemble[:, target_label_mapper['VCD_0'], :] + ensemble[:, target_label_mapper['VCD_1'], :] + ensemble[:, target_label_mapper['VCD_2'], :] + ensemble[:, target_label_mapper['VCD_3'], :]
        med = np.median(vcd_ens, axis=0)
        lo = np.percentile(vcd_ens, 2.5, axis=0)
        hi = np.percentile(vcd_ens, 97.5, axis=0)
        # plot individual trajectories lightly
        for rep in range(ensemble_size):
            ax.plot(t, vcd_ens[rep], color='tab:blue', alpha=0.12, lw=0.8)
        ax.plot(t, med, color='tab:blue', lw=1.5, label='Prediction')
        ax.fill_between(t, lo, hi, color='tab:blue', alpha=0.18, label='95% PI')
        meas_vals = data_extra_cont[b]['VCD_0'].values
    else:
        idx = target_label_mapper[key]
        # per-replicate trajectories
        for rep in range(ensemble_size):
            ax.plot(t, ensemble[rep, idx, :], color='tab:blue', alpha=0.12, lw=0.8)
        ax.plot(t, median[idx], color='tab:blue', lw=1.5, label='Prediction')
        ax.fill_between(t, lower[idx], upper[idx], color='tab:blue', alpha=0.18, label='95% PI')
        meas_vals = data_extra_cont[b][key].values

    # set y-label based on the type of variable
    if key == 'VCD':
        ax.set_ylabel(f'{extracellular_keys[key]} (million cells/mL)', fontsize=15)
    elif 'EANTI' in key.upper():  # IgG related key
        ax.set_ylabel(f'{extracellular_keys[key]} (g/L)', fontsize=15)
    else:
        ax.set_ylabel(f'{extracellular_keys[key]} Conc. (mM)', fontsize=15)
    ax.tick_params(axis='y', labelsize=12)

    # overlay measured points
    valid = ~np.isnan(meas_vals)
    if np.any(valid):
        ax.scatter(meas_days[valid], meas_vals[valid], color='tab:orange', zorder=5, label='Measurements')

    # show legend (will include Measurements if present)
    if i == 0:
        ax.legend(loc='best', fontsize=10)

    ax.set_xlabel('Time (day)', fontsize=15)
    ax.set_xticks(all_days)
    ax.set_xticklabels(all_days, fontsize=12)
    ax.set_xlim([all_days.min() - 1, all_days.max() + 1])
    ax.set_ylim(bottom=-0.05 * np.nanmax(meas_vals))
    ax.grid(True)
    ax.text(0.92, 0.1, f'({chr(97+i)})', transform=ax.transAxes, fontsize=15, fontweight='bold', va='top')


# hide any unused subplots
for j in range(n, len(axs)):
    axs[j].set_visible(False)
plt.tight_layout()
plt.savefig('Prediction_B3.png')
plt.show()


#%% minimizing WAPE of one metabolite only
def objective_function_single(params, target_key):
    total_wape = 0.0
    batch_count = 0
    for b in ['tigr20_6']:
        wape_dict, overall_wape, t, x = calc_wape_for_batch(b, params, mode_var=0)
        wape_single = wape_dict.get(target_key, np.nan)
        if not np.isnan(wape_single):
            total_wape += wape_single
            batch_count += 1
    if batch_count > 0:
        return total_wape / batch_count
    else:
        return np.inf


#%%
# prepare a copy of theta and make only specified params perturbable
params_single = theta.copy()

# EDIT this list to contain the parameter names you want to vary (exact names as in theta)
# GLC 
# perturbable_params = ['vmax2_0', 'vmax2_1', 'vmax2_2', 'vmax2_3', 'K_mEGLC_0', 'K_mEGLC_1', 'K_mEGLC_2', 'K_iLACtoHK_0', 'K_iLACtoHK_1', 'K_iLACtoHK_2']  # <- replace with your list
# LAC
# perturbable_params = ['vmax3f_0', 'vmax3f_1', 'vmax3f_2','K_mEGLC_0', 'K_mEGLC_1', 'K_mEGLC_2', 'vmax3r_0', 'vmax3r_1', 'vmax3r_2','K_mELAC_0', 'K_mELAC_1', 'K_mELAC_2']  # <- replace with your list
# SER
# perturbable_params = ['vmax7_0', 'vmax7_1', 'vmax7_2', 'vmax7_3', 'K_mESER_0', 'K_mESER_1', 'K_mESER_2', 'K_mESER_3']  # <- replace with your list
# ASP
# perturbable_params = ['vmax19_0', 'vmax19_1', 'vmax19_2', 'vmax19_3', 'K_mEASP_0', 'K_mEASP_1', 'K_mEASP_2', 'K_mEASP_3']  # <- replace with your list
# BIOM
perturbable_params = ['vmax30_0', 'vmax30_1', 'vmax30_2', 'vmax30_3', 'beta_01_2', 'beta_01_3', 'beta_12_2', 'beta_12_3', 'beta_23_2', 'beta_23_3']  # <- replace with your list

if not isinstance(perturbable_params, (list, tuple, set)):
    raise RuntimeError("perturbable_params must be a list/tuple/set of parameter names.")

all_names = list(params_single.keys())
missing = [p for p in perturbable_params if p not in all_names]
if missing:
    raise RuntimeError(f"The following perturbable parameter names are not present in theta: {missing}")

# set .vary only for listed parameters
for name in all_names:
    params_single[name].vary = (name in perturbable_params)

# wrapper to print progress for single metabolite objective
eval_counter = {'n': 0}
def objective_print_single(params):
    val = objective_function_single(params, target_key='VCD')
    eval_counter['n'] += 1
    print(f"Iteration {eval_counter['n']}: VCD error = {val}")
    return val

# run optimizer (only listed parameters will be changed)
try:
    first_msg = "', '".join(perturbable_params)
    print(f"Running optimizer: Nelder-Mead (varying: '{first_msg}')")
    eval_counter['n'] = 0
    res = minimize(objective_print_single, params_single, method='Nelder-Mead', options={'maxiter': 500, 'maxfev': 5000})
    val = objective_function_single(res.params, target_key='VCD')
    print(f"  -> finished Nelder-Mead, VCD value = {val}\n")
    best_res = res
    best_val = val
except Exception as exc:
    print(f"  -> Nelder-Mead failed: {exc}\n")
    best_res = None
    best_val = np.inf

#%%
# print perturbable parameter values (before optimization) and optimized values if available
print("Perturbable parameters (current params_single values):")
for name in perturbable_params:
    p = params_single[name]
    print(f"  {name}: value={p.value}, vary={p.vary}, min={p.min}, max={p.max}")

if 'best_res' in globals() and best_res is not None:
    print("\nOptimized values from best_res.params:")
    for name in perturbable_params:
        pv = best_res.params[name]
        print(f"  {name}: value={pv.value}")
else:
    print("\nNo optimization result available (best_res is None).")


