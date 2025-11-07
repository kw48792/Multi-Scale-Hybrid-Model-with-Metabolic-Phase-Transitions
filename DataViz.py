#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# Load data from the Excel file
# file_path = 'Rate.xlsx'
file_path = 'Rate_v2.xlsx'
# Loading rate data
tigr15_rate = pd.read_excel(file_path, sheet_name='Tigr15_Rate')
tigr20_rate = pd.read_excel(file_path, sheet_name='Tigr20_Rate')
tigr21_rate = pd.read_excel(file_path, sheet_name='Tigr21_Rate')

# Loading concentration data
tigr15_conc = pd.read_excel(file_path, sheet_name='Tigr15_Conc')
tigr20_conc = pd.read_excel(file_path, sheet_name='Tigr20_Conc')
tigr21_conc = pd.read_excel(file_path, sheet_name='Tigr21_Conc')



#%%
def create_combined_concentration_plot(conc_df_list, labels, metabolites):
    num_metabolites = len(metabolites)
    cols = 3
    rows = np.ceil(num_metabolites / cols)

    fig, axs = plt.subplots(int(rows), int(cols), figsize=(15, 3 * rows))
    axs = axs.flatten()

    offset = 0.2
    colors = ['blue', 'red', 'green']
    fmts = ['-o', '-s', '-D']

    # Work on copies to avoid mutating original data
    dfs = [df.copy() for df in conc_df_list]

    # Set Glutamine / Glutamate Mean and Std to NA for tigr21 at day 10 and 11
#    if len(dfs) >= 3:
#        df3 = dfs[2]
#        cols_to_na = ['Glutamine Mean', 'Glutamine Std'] #, 'Glutamate Mean', 'Glutamate Std']
#        if 'Day' in df3.columns:
#            mask = df3['Day'].isin([10, 11])
#            for c in cols_to_na:
#                if c in df3.columns:
#                    df3.loc[mask, c] = np.nan
#            dfs[2] = df3

    # Determine unique days across all datasets for x-ticks
    all_days = np.unique(np.concatenate([df['Day'] for df in dfs]))

    for i, metabolite in enumerate(metabolites):
        ax = axs[i]
        for j, (df, label) in enumerate(zip(dfs, labels)):
            # For tigr21 (Case C) exclude days 10 and 11 only for Glutamine.
            # For other metabolites adjust the day values slightly so the later isin([10,11]) check
            # does not remove those rows (keeps plotting intact).
            if j == 2 and metabolite != 'Glutamine':
                df = df.copy()
                mask = df['Day'].isin([10, 11])
                if mask.any():
                    df.loc[mask, 'Day'] = df.loc[mask, 'Day'].astype(float) + 1e-6
            if j == 2:
                plot_df = df.loc[~df['Day'].isin([10, 11])]
            else:
                plot_df = df
            
            if j == 0:  # Case A
                x = plot_df['Day'] - offset
            elif j == 1:  # Case B
                x = plot_df['Day']
            else:  # Case C
                x = plot_df['Day'] + offset

            ax.errorbar(x, plot_df[f'{metabolite} Mean'],
                        yerr=plot_df[f'{metabolite} Std']/np.sqrt(3), fmt=fmts[j % 3],
                        label=label, color=colors[j % len(colors)], capsize=5)

        ax.set_xlabel('Time (day)', fontsize=15)
        ax.set_xticks(all_days)
        ax.set_xticklabels(all_days, fontsize=12)
        ax.set_xlim([all_days.min() - 1, all_days.max() + 1])
        ax.set_ylim(bottom=-0.05 * np.max(plot_df[f'{metabolite} Mean']))
        if metabolite == 'VCD':
            ax.set_ylabel(f'{metabolite} (million cells/mL)', fontsize=15)
        elif metabolite == 'IgG':
            ax.set_ylabel(f'{metabolite} Conc. (g/L)', fontsize=15)
        else:
            ax.set_ylabel(f'{metabolite} Conc. (mM)', fontsize=15)
        ax.tick_params(axis='y', labelsize=12)
        if i == 0:
            ax.legend()
        ax.grid(True)
        ax.text(0.92, 0.1, f'({chr(97+i)})', transform=ax.transAxes, fontsize=15, fontweight='bold', va='top')

    for ax in axs[num_metabolites:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig('Concentration_combine_v3.png')
    plt.show()

# Get a list of metabolites from the column headers
metabolites = [col.replace(' Mean', '') for col in tigr15_rate.columns if ' Mean' in col]
metabolites_with = ['NH3', 'Alanine', 'Aspartate', 'Serine']
metabolites_without = [metabolite for metabolite in metabolites if metabolite not in metabolites_with]

# Call for all three cases
create_combined_concentration_plot(
    [tigr15_conc, tigr20_conc, tigr21_conc],
    ['Case A', 'Case B', 'Case C'],
    metabolites_without
)


#%%
#################################
# Loading phase specific flux rate data
tigr15_phase = pd.read_excel(file_path, sheet_name='Tigr15_Phase')
tigr20_phase = pd.read_excel(file_path, sheet_name='Tigr20_Phase')
tigr21_phase = pd.read_excel(file_path, sheet_name='Tigr21_Phase')

# Function to create combined plots for concentration data
def create_combined_phase_plot(phase_df1, phase_df2, phase_df3, metabolites):
    num_metabolites = len(metabolites)
    cols = 3  # 3 columns
    rows = np.ceil(num_metabolites / cols)  # calculate rows needed

    fig, axs = plt.subplots(int(rows), int(cols), figsize=(15, 3 * rows))
    axs = axs.flatten()  # Flatten the array for easy iteration

    for i, metabolite in enumerate(metabolites):
        ax = axs[i]
        index = np.arange(len(phase_df1['Phase']))

        ax.bar(index - 0.3, phase_df1[f'{metabolite} Mean'], yerr=phase_df1[f'{metabolite} Std']/np.sqrt(3),
               width=0.3, capsize=5, label='Case A', color='blue')
        ax.bar(index, phase_df2[f'{metabolite} Mean'], yerr=phase_df2[f'{metabolite} Std']/np.sqrt(3),
               width=0.3, capsize=5, label='Case B', color='red')
        ax.bar(index + 0.3, phase_df3[f'{metabolite} Mean'], yerr=phase_df3[f'{metabolite} Std']/np.sqrt(3),
               width=0.3, capsize=5, label='Case C', color='green')
        ax.set_xticks(index)
        ax.set_xticklabels(phase_df1['Phase'], fontsize=12)

        if metabolite == 'VCD':
            ax.set_ylabel('Growth Rate (h$^{-1}$)', fontsize=15)
        elif metabolite == 'IgG':
            ax.set_ylabel(f'{metabolite} Flux \n (mg/million cell·h)', fontsize=15)
        else:
            ax.set_ylabel(f'{metabolite} Flux \n ($\mu$mol/million cell·h)', fontsize=15)
        ax.tick_params(axis='y', labelsize=12)
        if i == 0:
            ax.legend()
        ax.grid(True)

        # Add subplot index (a), (b), (c),...
        ax.text(0.9, 0.32, f'({chr(97+i)})', transform=ax.transAxes, fontsize=15, fontweight='bold', va='top')

    # Hide any unused subplots
    for ax in axs[num_metabolites:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig('Phase_combine_v3.png')
    plt.show()

create_combined_phase_plot(tigr15_phase, tigr20_phase, tigr21_phase, metabolites_without)


#%%
# Function to create combined concentration and phase plot
def create_combined_plot(conc_df_list, phase_df_list, labels, metabolites):
    num_metabolites = len(metabolites)
    cols = 2  # 2 columns for concentration and phase
    rows = num_metabolites  # Each metabolite in a separate row

    fig, axs = plt.subplots(int(rows), int(cols), figsize=(10, 12))
    axs = axs.flatten()

    offset = 0.3
    colors = ['blue', 'red', 'green']
    fmts = ['-o', '-s', '-D']

    # Work on copies to avoid mutating original data
    dfs = [df.copy() for df in conc_df_list]

    # # Set NH3 Mean and Std to NA for tigr21 (assumed index 2) at day 10 and 11
    # if len(dfs) >= 3:
    #     df3 = dfs[2]
    #     cols_to_na = ['NH3 Mean', 'NH3 Std']
    #     if 'Day' in df3.columns:
    #         mask = df3['Day'].isin([10, 11])
    #         for c in cols_to_na:
    #             if c in df3.columns:
    #                 df3.loc[mask, c] = np.nan
    #         dfs[2] = df3

    # Determine unique days across all concentration datasets for x-ticks
    all_days = np.unique(np.concatenate([df['Day'] for df in dfs]))

    for i, metabolite in enumerate(metabolites):
        # Concentration plot
        ax_conc = axs[i * 2]
        for j, (df, label) in enumerate(zip(dfs, labels)):
            # For tigr21 (Case C) exclude days 10 and 11 only for NH3.
            # For other metabolites adjust the day values slightly so the later isin([10,11]) check
            # does not remove those rows (keeps plotting intact).
            if j == 2 and metabolite != 'NH3':
                df = df.copy()
                mask = df['Day'].isin([10, 11])
                if mask.any():
                    df.loc[mask, 'Day'] = df.loc[mask, 'Day'].astype(float) + 1e-6
            if j == 2:
                plot_df = df.loc[~df['Day'].isin([10, 11])]
            else:
                plot_df = df

            if j == 0:
                x = plot_df['Day'] - offset
            elif j == 1:
                x = plot_df['Day']
            else:
                x = plot_df['Day'] + offset
            ax_conc.errorbar(x, plot_df[f'{metabolite} Mean'],
                             yerr=plot_df[f'{metabolite} Std']/np.sqrt(3), fmt=fmts[j % 3],
                             label=label, color=colors[j % len(colors)], capsize=5)
        ax_conc.set_xlabel('Time (day)', fontsize=15)
        ax_conc.set_xticks(all_days)
        ax_conc.set_xticklabels(all_days, fontsize=12)
        ax_conc.set_xlim([all_days.min() - 1, all_days.max() + 1])
        if metabolite == 'VCD':
            ax_conc.set_ylabel(f'{metabolite} (million cells/mL)', fontsize=15)
        elif metabolite == 'NH3':
            ax_conc.set_ylabel('Ammonia Conc. (mM)', fontsize=15)
        else:
            ax_conc.set_ylabel(f'{metabolite} Conc. (mM)', fontsize=15)
        ax_conc.tick_params(axis='y', labelsize=12)
        if i == 0:
            ax_conc.legend()
        ax_conc.grid(True)
        ax_conc.text(0.9, 0.1, f'({chr(97+i)})', transform=ax_conc.transAxes, fontsize=15, fontweight='bold', va='top')
        ax_conc.set_ylim(bottom=-0.05 * np.max(plot_df[f'{metabolite} Mean']))
        # Phase plot
        ax_phase = axs[i * 2 + 1]
        index = np.arange(len(phase_df_list[0]['Phase']))
        for j, (df, label) in enumerate(zip(phase_df_list, labels)):
            if j == 0:
                x = index - offset
            elif j == 1:
                x = index
            else:
                x = index + offset
            ax_phase.bar(x, df[f'{metabolite} Mean'],
                         yerr=df[f'{metabolite} Std']/np.sqrt(3),
                         width=0.3, capsize=5, label=label, color=colors[j % len(colors)])
        ax_phase.set_xlabel('Phase', fontsize=15)
        ax_phase.set_xticks(index)
        ax_phase.set_xticklabels(phase_df_list[0]['Phase'], fontsize=12)
        if metabolite == 'VCD':
            ax_phase.set_ylabel('Growth Rate (h$^{-1}$)', fontsize=15)
        elif metabolite == 'NH3':
            ax_phase.set_ylabel('Ammonia Flux \n ($\mu$mol/million cell·h)', fontsize=15)
        else:
            ax_phase.set_ylabel(f'{metabolite} Flux \n ($\mu$mol/million cell·h)', fontsize=15)
        ax_phase.tick_params(axis='y', labelsize=12)
        if i == 0:
            ax_phase.legend()
        ax_phase.grid(True)
        ax_phase.text(0.9, 0.45, f'({chr(101+i)})', transform=ax_phase.transAxes, fontsize=15, fontweight='bold', va='top')

    # Hide any unused subplots
    for ax in axs[num_metabolites * 2:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig('Combined_plot_v3.png')
    plt.show()

# Call the function with the required arguments
create_combined_plot(
    [tigr15_conc, tigr20_conc, tigr21_conc],
    [tigr15_phase, tigr20_phase, tigr21_phase],
    ['Case A', 'Case B', 'Case C'],
    metabolites_with
)



#%%
#################################
# Loading day specific OUR data
tigr_our_all = pd.read_excel(file_path, sheet_name='Day_OUR')

def create_day_OUR_plot(tigr_our):
    # Ensure the figure can accommodate the number of phases
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Generate index for bar positions
    index = np.arange(len(tigr_our['Day']))
    
    ax.bar(index - 0.3, tigr_our['qO2_Tigr15_Mean'], yerr=tigr_our['qO2_Tigr15_Std']/np.sqrt(3), 
           width=0.3, capsize=5, label='Case A', color='blue')
    ax.bar(index , tigr_our['qO2_Tigr20_Mean'], yerr=tigr_our['qO2_Tigr20_Std']/np.sqrt(3), 
           width=0.3, capsize=5, label='Case B', color='red')
    ax.bar(index + 0.3, tigr_our['qO2_Tigr21_Mean'], yerr=tigr_our['qO2_Tigr21_Std']/np.sqrt(3), 
           width=0.3, capsize=5, label='Case C', color='green')
    # ax.set_xlabel('Phase', fontsize=18)
    ax.set_xticks(index)
    # wrapped_labels = ['\n'.join(textwrap.wrap(label, width=10)) for label in tigr_our['Phase']]
    # ax.set_xticklabels(wrapped_labels, fontsize=15)
    
    ax.set_xticklabels(tigr_our['Day'], fontsize=15)
    ax.set_xlabel('Time (day)', fontsize=24)
    ax.set_ylabel('qO$_2$ (pg/cell·h)', fontsize=24)

    ax.legend(loc='best', fontsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.grid(True)
    plt.savefig('OUR_Day.png')
    plt.show()

# Generate and display plots for OUR data
create_day_OUR_plot(tigr_our_all)










####Discarded code below####

#%%
# Assuming file_path and loading of data as per your existing code

def create_combined_concentration_plot(conc_df1, conc_df2, metabolites):
    num_metabolites = len(metabolites)
    cols = 3 
    rows = np.ceil(num_metabolites / cols)  # Calculate rows needed
    
    fig, axs = plt.subplots(int(rows), int(cols), figsize=(15, 3 * rows))
    axs = axs.flatten()  # Flatten the array for easy iteration
    
    offset = 0.2  # Define an offset for the x-values
    
    # Determine unique days across both datasets for x-ticks
    all_days = np.union1d(conc_df1['Day'], conc_df2['Day'])
    
    for i, metabolite in enumerate(metabolites):
        ax = axs[i]
        # Apply offset to Case A and Case B x-values for clarity
        ax.errorbar(conc_df1['Day'] - offset, conc_df1[f'{metabolite} Mean'], 
                    yerr=conc_df1[f'{metabolite} Std']/np.sqrt(3), fmt='-o', 
                    label='Case A', color='blue', capsize=5)
        ax.errorbar(conc_df2['Day'] + offset, conc_df2[f'{metabolite} Mean'], 
                    yerr=conc_df2[f'{metabolite} Std']/np.sqrt(3), fmt='-s', 
                    label='Case B', color='red', capsize=5)
        
        ax.set_xlabel('Time (day)', fontsize=15)  
        ax.set_xticks(all_days)  # Set x-ticks to all unique days
        ax.set_xticklabels(all_days, fontsize=12)  # Set x-tick labels with fontsize
        ax.set_xlim([all_days.min() - 1, all_days.max() + 1])  # Optional: adjust x-limits for better visualization
        
        if metabolite == 'VCD':
            ax.set_ylabel(f'{metabolite} (million cells/mL)', fontsize=15)
        else:
            ax.set_ylabel(f'{metabolite} Conc. (mM)', fontsize=15)
        ax.tick_params(axis='y', labelsize=12)
        ax.legend()
        ax.grid(True)
        
        # Add subplot index (a), (b), (c),...
        ax.text(0.92, 0.1, f'({chr(97+i)})', transform=ax.transAxes, fontsize=15, fontweight='bold', va='top')
    
    # Hide any unused subplots
    for ax in axs[num_metabolites:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig('Concentration_combine_v2.png')
    plt.show()

# Get a list of metabolites from the column headers
metabolites = [col.replace(' Mean', '') for col in tigr15_rate.columns if ' Mean' in col]
metabolites_with = ['NH3', 'Alanine', 'Aspartate', 'Serine']
metabolites_without = [metabolite for metabolite in metabolites if metabolite not in metabolites_with]

# create_combined_concentration_plot(tigr15_conc, tigr20_conc, metabolites_with)
create_combined_concentration_plot(tigr15_conc, tigr20_conc, metabolites_without)



#%%
#################################
# Loading day specific flux rate data
tigr15_day = pd.read_excel(file_path, sheet_name='Tigr15_Rate')
tigr20_day = pd.read_excel(file_path, sheet_name='Tigr20_Rate')


# Function to create combined plots for concentration data
def create_combined_day_plot(day_df1, day_df2, metabolites):
    # Determine the layout size based on the number of metabolites
    num_metabolites = len(metabolites)
    cols = 3  # 3 columns
    rows = np.ceil(num_metabolites / cols)  # calculate rows needed
    
    fig, axs = plt.subplots(int(rows), int(cols), figsize=(15, 3 * rows))
    axs = axs.flatten()  # Flatten the array for easy iteration
    
    for i, metabolite in enumerate(metabolites):
        ax = axs[i]
        index = np.arange(len(day_df1['Day']))
        
        ax.bar(index - 0.2, day_df1[f'{metabolite} Mean'], yerr=day_df1[f'{metabolite} Std']/np.sqrt(3), 
               width=0.4, capsize=5, label='Case A', color='blue')
        ax.bar(index + 0.2, day_df2[f'{metabolite} Mean'], yerr=day_df2[f'{metabolite} Std']/np.sqrt(3), 
               width=0.4, capsize=5, label='Case B', color='red')
        # ax.set_title(metabolite)
        # ax.set_xlabel('Phase', fontsize=18)
        ax.set_xticks(index)
        ax.set_xticklabels(day_df1['Day'], fontsize=12)
        ax.set_xlabel('Time (day)', fontsize=15)  
        if metabolite == 'VCD':
            ax.set_ylabel('Growth Rate (h$^{-1}$)', fontsize=15)
        else:
            # ax.set_ylabel(f'{metabolite} Uptake/Secretion Flux \n (mmol/cell·h)', fontsize=15)
            ax.set_ylabel(f'{metabolite} Flux \n (mmol/cell·h)', fontsize=15)
        ax.tick_params(axis='y', labelsize=12)
        ax.legend()
        ax.grid(True)
            
        # Add subplot index (a), (b), (c),...
        ax.text(0.9, 0.35, f'({chr(97+i)})', transform=ax.transAxes, fontsize=15, fontweight='bold', va='top')
    
    # Hide any unused subplots
    for ax in axs[num_metabolites:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig('Day_combine_v2.png')
    plt.show()

create_combined_day_plot(tigr15_day, tigr20_day, metabolites_without)




#%%
# Function to create combined concentration and day plot
def create_combined_plot(conc_df1, conc_df2, day_df1, day_df2, metabolites):
    num_metabolites = len(metabolites)
    cols = 2  # 2 columns for concentration and phase
    rows = num_metabolites  # Each metabolite in a separate row
    
    fig, axs = plt.subplots(int(rows), int(cols), figsize=(10, 12))
    axs = axs.flatten()  # Flatten the array for easy iteration
    
    offset = 0.2  # Define an offset for the x-values
    
    # Determine unique days across both datasets for x-ticks
    all_days = np.union1d(conc_df1['Day'], conc_df2['Day'])
    
    for i, metabolite in enumerate(metabolites):
        # Concentration plot
        ax_conc = axs[i * 2]
        ax_conc.errorbar(conc_df1['Day'] - offset, conc_df1[f'{metabolite} Mean'], 
                         yerr=conc_df1[f'{metabolite} Std']/np.sqrt(3), fmt='-o', 
                         label='Case A', color='blue', capsize=5)
        ax_conc.errorbar(conc_df2['Day'] + offset, conc_df2[f'{metabolite} Mean'], 
                         yerr=conc_df2[f'{metabolite} Std']/np.sqrt(3), fmt='-s', 
                         label='Case B', color='red', capsize=5)
        ax_conc.set_xlabel('Time (day)', fontsize=15)  
        ax_conc.set_xticks(all_days)  # Set x-ticks to all unique days
        ax_conc.set_xticklabels(all_days, fontsize=12)  # Set x-tick labels with fontsize
        ax_conc.set_xlim([all_days.min() - 1, all_days.max() + 1])  # Optional: adjust x-limits for better visualization
        if metabolite == 'VCD':
            ax_conc.set_ylabel(f'{metabolite} (million cells/mL)', fontsize=15)
        else:
            ax_conc.set_ylabel(f'{metabolite} Conc. (mM)', fontsize=15)
        ax_conc.tick_params(axis='y', labelsize=12)
        ax_conc.legend()
        ax_conc.grid(True)
        ax_conc.text(0.9, 0.1, f'({chr(97+i)})', transform=ax_conc.transAxes, fontsize=15, fontweight='bold', va='top')
        
        # Day plot
        ax_phase = axs[i * 2 + 1]
        index = np.arange(len(day_df1['Day']))
        ax_phase.bar(index - 0.3, day_df1[f'{metabolite} Mean'], yerr=day_df1[f'{metabolite} Std']/np.sqrt(3), 
                      width=0.4, capsize=5, label='Case A', color='blue')
        ax_phase.bar(index + 0.2, day_df2[f'{metabolite} Mean'], yerr=day_df2[f'{metabolite} Std']/np.sqrt(3), 
                      width=0.4, capsize=5, label='Case B', color='red')
        ax_phase.set_xlabel('Day', fontsize=15)
        ax_phase.set_xticks(index)
        ax_phase.set_xticklabels(day_df1['Day'], fontsize=12)
        if metabolite == 'VCD':
            ax_phase.set_ylabel('Growth Rate (h$^{-1}$)', fontsize=15)
        else:
            ax_phase.set_ylabel(f'{metabolite} Flux \n (mmol/cell·h)', fontsize=15)
        ax_phase.tick_params(axis='y', labelsize=12)
        ax_phase.legend()
        ax_phase.grid(True)
        ax_phase.text(0.9, 0.35, f'({chr(101+i)})', transform=ax_phase.transAxes, fontsize=15, fontweight='bold', va='top')
    
    # Hide any unused subplots
    for ax in axs[num_metabolites * 2:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig('Combined_plot_day_v2.png')
    plt.show()

# Call the function with the required arguments
create_combined_plot(tigr15_conc, tigr20_conc, tigr15_day, tigr20_day, metabolites_with)




#%%
#################################
# Loading phase specific OUR data
tigr_our_all = pd.read_excel(file_path, sheet_name='Phase_OUR')

def create_phase_OUR_plot(tigr_our):
    # Ensure the figure can accommodate the number of phases
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Generate index for bar positions
    index = np.arange(len(tigr_our['Phase']))
    
    ax.bar(index - 0.2, tigr_our['qO2_Tigr15_Mean'], yerr=tigr_our['qO2_Tigr15_Std']/np.sqrt(3), 
           width=0.4, capsize=5, label='Case A', color='blue')
    ax.bar(index + 0.2, tigr_our['qO2_Tigr20_Mean'], yerr=tigr_our['qO2_Tigr20_Std']/np.sqrt(3), 
           width=0.4, capsize=5, label='Case B', color='red')
    
    # ax.set_xlabel('Phase', fontsize=18)
    ax.set_xticks(index)
    # wrapped_labels = ['\n'.join(textwrap.wrap(label, width=10)) for label in tigr_our['Phase']]
    # ax.set_xticklabels(wrapped_labels, fontsize=15)
    
    ax.set_xticklabels(tigr_our['Phase'], fontsize=15)
    
    ax.set_ylabel('qO$_2$ (pg/cell·h)', fontsize=18)
    
    ax.legend(loc='best', fontsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.grid(True)
    plt.savefig('OUR_Phase.png')
    plt.show()

# Generate and display plots for OUR data
create_phase_OUR_plot(tigr_our_all)



#%%
#################################
# Function to create two-level plots using gridspec
def create_two_level_plots(rate_df1, rate_df2, conc_df1, conc_df2, metabolite):
    gs = gridspec.GridSpec(2, 1)
    fig = plt.figure(figsize=(10, 8))

    # Upper plot for concentration data
    ax_conc = fig.add_subplot(gs[0])
    ax_conc.plot(conc_df1['Day'], conc_df1[f'{metabolite} Mean'], label='Tigr15', color='blue', marker='o')
    ax_conc.plot(conc_df2['Day'], conc_df2[f'{metabolite} Mean'], label='Tigr20', color='orange', marker='x')
    if metabolite == 'VCD':
        ax_conc.set_ylabel(f'{metabolite} (million cells/mL)', fontsize=18)
    else:
        ax_conc.set_ylabel(f'{metabolite} Conc. (mM)', fontsize=18)
    ax_conc.legend(loc='upper left')
    ax_conc.grid(True)
    ax_conc.tick_params(axis='y', labelsize=15)
    ax_conc.get_yaxis().set_label_coords(-0.1, 0.5)

    # Lower plot for rate data
    ax_rate = fig.add_subplot(gs[1], sharex=ax_conc)
    ax_rate.bar(rate_df1['Day'] + 0.5 - 0.2, rate_df1[f'{metabolite} Mean'], yerr=rate_df1[f'{metabolite} Std'],
                width=0.4, capsize=5, label='Tigr15', color='blue')
    ax_rate.bar(rate_df2['Day'] + 0.5 + 0.2, rate_df2[f'{metabolite} Mean'], yerr=rate_df2[f'{metabolite} Std'],
                width=0.4, capsize=5, label='Tigr20', color='orange')
    ax_rate.set_xlabel('Day', fontsize=18)
    if metabolite == 'VCD':
        ax_rate.set_ylabel('Rate (h$^{-1}$)', fontsize=18)
    else:
        ax_rate.set_ylabel('Rate (mmol/cell·h)', fontsize=18)
    ax_rate.legend(loc='best')
    ax_rate.tick_params(axis='y', labelsize=15)
    ax_rate.grid(True)
    ax_rate.get_yaxis().set_label_coords(-0.1, 0.5)

    # Set the X-axis to display each day as a step
    all_days = sorted(set(conc_df1['Day']).union(conc_df2['Day']))
    # ax_conc.set_xticks(all_days)
    # ax_rate.set_xticks(all_days)
    plt.xticks(all_days, fontsize=15)

    # Adjust the spacing between the plots
    plt.subplots_adjust(hspace=0.03)

    # Save the figure
    # plt.savefig(f'{metabolite}.png')
    plt.show()


# Generate and display two-level plots using gridspec
for metabolite in metabolites:
    create_two_level_plots(tigr15_rate, tigr20_rate, tigr15_conc, tigr20_conc, metabolite)
    

# Function to create plots for concentration data
def create_concentration_plot(conc_df1, conc_df2, metabolite):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(conc_df1['Day'], conc_df1[f'{metabolite} Mean'], yerr=conc_df1[f'{metabolite} Std']/np.sqrt(3), 
               fmt='-o', label='Tigr15', color='blue', capsize=5)
    ax.errorbar(conc_df2['Day'], conc_df2[f'{metabolite} Mean'], yerr=conc_df2[f'{metabolite} Std']/np.sqrt(3), 
               fmt='-x', label='Tigr20', color='orange', capsize=5)
    if metabolite == 'VCD':
        ax.set_ylabel(f'{metabolite} (million cells/mL)', fontsize=18)
    else:
        ax.set_ylabel(f'{metabolite} Conc. (mM)', fontsize=18)
    ax.set_xlabel('Day', fontsize=18)
    ax.legend(loc='upper left', fontsize=15)
    ax.grid(True)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # plt.title(f'{metabolite} Concentration Over Time', fontsize=20)
    # plt.savefig(f'{metabolite}_Concentration.png')
    plt.show()

# Generate and display plots for concentration data
for metabolite in metabolites:
    create_concentration_plot(tigr15_conc, tigr20_conc, metabolite)








#%%
# Function to create plots for rate data
def create_rate_plot(rate_df1, rate_df2, metabolite):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(rate_df1['Day'] + 0.5 - 0.2, rate_df1[f'{metabolite} Mean'], yerr=rate_df1[f'{metabolite} Std'],
                width=0.4, capsize=5, label='Tigr15', color='blue')
    ax.bar(rate_df2['Day'] + 0.5 + 0.2, rate_df2[f'{metabolite} Mean'], yerr=rate_df2[f'{metabolite} Std'],
                width=0.4, capsize=5, label='Tigr20', color='orange')
    ax.set_xlabel('Day', fontsize=18)
    if metabolite == 'VCD':
        ax.set_ylabel('Growth Rate (h$^{-1}$)', fontsize=18)
    else:
        ax.set_ylabel(f'{metabolite} Uptake/Secretion Flux \n (mmol/cell·h)', fontsize=18)
    ax.legend(loc='best', fontsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.grid(True)
    
    # Set the X-axis to display each day as a step
    all_days = sorted(set(rate_df1['Day']).union(rate_df2['Day']))
    plt.xticks(all_days, fontsize=15)
    plt.yticks(fontsize=15)
    # plt.title(f'{metabolite} Concentration Over Time', fontsize=20)
    # plt.savefig(f'{metabolite}_Rate.png')
    plt.show()

# Generate and display plots for concentration data
for metabolite in metabolites:
    create_rate_plot(tigr15_rate, tigr20_rate, metabolite)



#%%
# Loading rate data
tigr15_phase = pd.read_excel(file_path, sheet_name='Tigr15_Phase')
tigr20_phase = pd.read_excel(file_path, sheet_name='Tigr20_Phase')

def create_phase_plot(phase_df1, phase_df2, metabolite):
    # Ensure the figure can accommodate the number of phases
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Generate index for bar positions
    index = np.arange(len(phase_df1['Phase']))
    
    ax.bar(index - 0.2, phase_df1[f'{metabolite} Mean'], yerr=phase_df1[f'{metabolite} Std']/np.sqrt(3), 
           width=0.4, capsize=5, label='Tigr15', color='blue')
    ax.bar(index + 0.2, phase_df2[f'{metabolite} Mean'], yerr=phase_df2[f'{metabolite} Std']/np.sqrt(3), 
           width=0.4, capsize=5, label='Tigr20', color='orange')
    
    ax.set_xlabel('Phase', fontsize=18)
    ax.set_xticks(index, fontsize=18)
    ax.set_xticklabels(phase_df1['Phase'], rotation=45)
    
    if metabolite == 'VCD':
        ax.set_ylabel('Growth Rate (h$^{-1}$)', fontsize=18)
    else:
        ax.set_ylabel(f'{metabolite} Uptake/Secretion Flux \n (mmol/cell·h)', fontsize=18)
    
    ax.legend(loc='best', fontsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.grid(True)
    # plt.savefig(f'{metabolite}_Phase.png')
    plt.show()

# Generate and display plots for concentration data
for metabolite in metabolites:
    create_phase_plot(tigr15_phase, tigr20_phase, metabolite)