# import required modules
import pandas as pd
import numpy as np
import scipy
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib
import os                           # to get the current working directory / filepath
from fit_DSL import *

import warnings
warnings.filterwarnings("ignore")
pwd = os.path.abspath(os.getcwd())  # present working directory

## Read in the data
CO2_df = pd.read_csv('CO2_df.csv')
N2_df  = pd.read_csv('N2_df.csv')


## Get useful variables and statistics
unique_materials = CO2_df['Name'].unique()
col_headers = CO2_df.columns.values.tolist()
P_key = col_headers[0]
L_key = col_headers[1]
T_key = col_headers[2]
ID_key = col_headers[3]
error_key =  col_headers[4]
name_key = col_headers[5]

## create dictionary for ID number
ID_num = {1: 'UFF_DDEC',
          2: 'UFF_EQeq',
          3: 'UFF_Qeq',
          4: 'UFF_Neutral',
          5: 'DRE_DDEC',
          6: 'DRE_EQeq',
          7: 'DRE_Qeq',
          8: 'DRE_Neutral'}


def quick_plot(q_s1, q_s2, b01, b02, delta_H1, delta_H2, P, T=298):
    site_1 = q_s1 * (P * b01 * np.exp(delta_H1 / (8.314 * T))) / (1 + (P * b01 * np.exp(delta_H1 / (8.314 * T))))
    site_2 = q_s2 * (P * b02 * np.exp(delta_H2 / (8.314 * T))) / (1 + (P * b02 * np.exp(delta_H2 / (8.314 * T))))
    return (site_1 + site_2)


for material in unique_materials:

    material_data = CO2_df[CO2_df['Name'] == material]
    print(material)
    print('CO2')

    ## Statistics
    sources = pd.unique(material_data[ID_key])  # source identifiers (i.e. number of unique ID numbers)
    no_sources = sources.size  # no. of sources
    uT = pd.unique(material_data[T_key])  # unique temperatures
    T = uT[np.argsort(uT)]  # sort from low to high, so we ensure first entry is lowest temperature
    no_T = uT.size  # number of temperatures to evaluate

    # Partition by forcefield type
    by_researcher = []
    for j in sources:
        #     print(j)
        m = material_data[material_data[ID_key] == j]
        by_researcher.append(m)

        # plot and perform regressions
    cols = 2;
    rows = 4
    fig, ax = plt.subplots(rows, cols, figsize=(8, 10))
    CO2_dsl = []

    l = 0
    next_researcher = 0

    for i in range(rows):
        for j in range(cols):
            if len(by_researcher) == l:
                break
            else:
                data_indv = by_researcher[
                    next_researcher]  # get dataframe for an individual ID (forcefield combination)
                data_indv.reset_index(drop=True, inplace=True)

                # partition data by temperature
                T273 = data_indv[data_indv[T_key] == 273]
                T298 = data_indv[data_indv[T_key] == 298]
                T323 = data_indv[data_indv[T_key] == 323]

                ## Get initial parameter guess
                try:
                    init_params = CO2_DSL_regress(T273)
                    params = init_params.params.valuesdict()
                    ## use init_params to instantiate temperature dependent DSL model
                    opt_params = CO2_DSL_regress(data_indv, init_guess=params)
                    x = opt_params.params.valuesdict()

                    CO2_dsl.append([x['qs1'], x['qs2'], x['b01'], x['b02'], x['deltaH1'], x['deltaH2'], material,
                                    ID_num[data_indv[ID_key].iloc[0]]])
                    # plot the data by temperature
                    ax[i][j].plot(np.linspace(0, 10, 100000),
                                  quick_plot(x['qs1'], x['qs2'], x['b01'], x['b02'], x['deltaH1'], x['deltaH2'],
                                             np.linspace(0, 10, 100000), 273), '-b')
                    ax[i][j].plot(T273[P_key], T273[L_key], 'o', c='b', ms=4, label='273K')

                    ax[i][j].plot(np.linspace(0, 10, 100000),
                                  quick_plot(x['qs1'], x['qs2'], x['b01'], x['b02'], x['deltaH1'], x['deltaH2'],
                                             np.linspace(0, 10, 100000), 298), '-m')
                    ax[i][j].plot(T298[P_key], T298[L_key], 'o', c='m', ms=4, label='298K')

                    ax[i][j].plot(np.linspace(0, 10, 100000),
                                  quick_plot(x['qs1'], x['qs2'], x['b01'], x['b02'], x['deltaH1'], x['deltaH2'],
                                             np.linspace(0, 10, 100000), 323), '-r')
                    ax[i][j].plot(T323[P_key], T323[L_key], 'o', c='r', ms=4, label='323K')

                    ax[i][j].set_xscale('log')
                    #             ax[i][j].set_yscale('log')
                    ax[i][j].set(title=material + ', ' + ID_num[data_indv[ID_key].iloc[0]])
                    ax[i][j].legend()
                    ax[i][j].grid()
                except:
                    ## Get initial parameter guess
                    init_params = CO2_DSL_regress_BACKUP(T273, tol=1e-10, norm_factor=1)

                    ## use init_params to instantiate temperature dependent DSL model
                    opt_params = CO2_DSL_regress_BACKUP(data_indv, init_guess=init_params.x, tol=1e-12, norm_factor=50,
                                                        maxiter=200)
                    x = opt_params.x
                    CO2_dsl.append([x[0], x[1], x[2], x[3], x[4], x[5], material, ID_num[data_indv[ID_key].iloc[0]]])

                    # plot the data by temperature
                    ax[i][j].plot(np.linspace(0, 10, 100000),
                                  quick_plot(x[0], x[1], x[2], x[3], x[4], x[5], np.linspace(0, 10, 100000), 273), '-b')
                    ax[i][j].plot(T273[P_key], T273[L_key], 'o', c='b', ms=4, label='273K')

                    ax[i][j].plot(np.linspace(0, 10, 100000),
                                  quick_plot(x[0], x[1], x[2], x[3], x[4], x[5], np.linspace(0, 10, 100000), 298), '-m')
                    ax[i][j].plot(T298[P_key], T298[L_key], 'o', c='m', ms=4, label='298K')

                    ax[i][j].plot(np.linspace(0, 10, 100000),
                                  quick_plot(x[0], x[1], x[2], x[3], x[4], x[5], np.linspace(0, 10, 100000), 323), '-r')
                    ax[i][j].plot(T323[P_key], T323[L_key], 'o', c='r', ms=4, label='323K')

                    ax[i][j].set_xscale('log')
                    #             ax[i][j].set_yscale('log')
                    ax[i][j].set(title=material + ', ' + ID_num[data_indv[ID_key].iloc[0]])
                    ax[i][j].legend()
                    ax[i][j].grid()

                L_lim = np.max(material_data[L_key])  # for determining a common y-axis limit for every graph
                #             ax[i][j].set_ylim([0, L_lim])

                ## Printing output to screen
                #             print(material+', '+ID_num[data_indv[ID_key].iloc[0]])
                #             print('Optimised DSL Model Parameters: \n q_s1 = {}\n q_s2 = {}\n b_01 = {}\n b_02 = {}\n deltaH_1 = {}\n deltaH_2 = {}\n\n'.format(x[0],x[1],x[2],x[3],x[4],x[5]))
                next_researcher += 1
            l += 1
    fig.supxlabel('Pressure [bar]')
    fig.supylabel('$CO_2$ Loading [mmol/g]')
    fig.tight_layout()
    fig.savefig(os.path.join(pwd, 'visualise_fits', material + '_CO2.png'), facecolor='w')
    plt.close()

    # get DSL model params in a dataframe
    col_names = ["q_s1", "q_s2", "b01", "b02", "deltaH1", "deltaH2", 'Material', 'Forcefield']
    CO2_DSL_params = pd.DataFrame(CO2_dsl, columns=col_names)
    CO2_DSL_params.to_csv(os.path.join(pwd, 'DSL_model_values', material + '_CO2_DSLparams.csv'), index=False)

    # ******************************************************* N2 ******************************************
    print('N2')
    material_data = N2_df[N2_df['Name'] == material]
    by_researcher_N2 = []
    for j in sources:
        m = material_data[material_data[ID_key] == j]
        by_researcher_N2.append(m)

    cols = 2;
    rows = 4
    fig, ax = plt.subplots(rows, cols, figsize=(8, 10))
    N2_dsl = []

    l = 0
    next_researcher = 0
    count = 0;
    for i in range(rows):
        for j in range(cols):
            if len(by_researcher_N2) == l:
                break
            else:
                data_indv = by_researcher_N2[
                    next_researcher]  # get dataframe for an individual ID (forcefield combination)
                data_indv.reset_index(drop=True, inplace=True)
                indv_CO2_DSL_params = CO2_DSL_params[CO2_DSL_params['Forcefield'] == ID_num[data_indv['ID'].iloc[0]]]
                # partition data by temperature
                T273 = data_indv[data_indv[T_key] == 273]
                T298 = data_indv[data_indv[T_key] == 298]
                T323 = data_indv[data_indv[T_key] == 323]

                ## Get initial parameter guess
                init_params = N2_DSL_regress(T273, indv_CO2_DSL_params, tol=1e-9)

                ## use init_params to instantiate temperature dependent DSL model
                opt_params = N2_DSL_regress(data_indv, indv_CO2_DSL_params, init_guess=init_params.x, tol=1e-15)
                x = opt_params.x

                N2_dsl.append(
                    [indv_CO2_DSL_params['q_s1'].iloc[0], indv_CO2_DSL_params['q_s2'].iloc[0], x[0], x[1], x[2], x[3],
                     material, ID_num[data_indv[ID_key].iloc[0]]])

                num_pts = 10000
                # plot the data by temperature
                ax[i][j].plot(np.linspace(0, 1, num_pts),
                              quick_plot(indv_CO2_DSL_params['q_s1'].iloc[0], indv_CO2_DSL_params['q_s2'].iloc[0], x[0],
                                         x[1], x[2], x[3], np.linspace(0, 1, num_pts), 273), '-b')
                ax[i][j].plot(T273[P_key], T273[L_key], 'o', c='b', ms=4, label='273K')

                ax[i][j].plot(np.linspace(0, 1, num_pts),
                              quick_plot(indv_CO2_DSL_params['q_s1'].iloc[0], indv_CO2_DSL_params['q_s2'].iloc[0], x[0],
                                         x[1], x[2], x[3], np.linspace(0, 1, num_pts), 298), '-m')
                ax[i][j].plot(T298[P_key], T298[L_key], 'o', c='m', ms=4, label='298K')

                ax[i][j].plot(np.linspace(0, 1, num_pts),
                              quick_plot(indv_CO2_DSL_params['q_s1'].iloc[0], indv_CO2_DSL_params['q_s2'].iloc[0], x[0],
                                         x[1], x[2], x[3], np.linspace(0, 1, num_pts), 323), '-r')
                ax[i][j].plot(T323[P_key], T323[L_key], 'o', c='r', ms=4, label='323K')

                #             ax[i][j].set_xscale('log')
                ax[i][j].set(title=material + ', ' + ID_num[data_indv[ID_key].iloc[0]])
                ax[i][j].legend()
                ax[i][j].grid()

                L_lim = np.max(material_data[L_key]) * 1.1  # for determining a common y-axis limit for every graph
                ax[i][j].set_ylim([0, L_lim])

                ## Printing output to screen
                #             print(material+', '+ID_num[data_indv[ID_key].iloc[0]])
                #             print('Optimised DSL Model Parameters: \n q_s1 = {}\n q_s2 = {}\n b_01 = {}\n b_02 = {}\n deltaH_1 = {}\n deltaH_2 = {}\n\n'.format(x[0],x[1],x[2],x[3],x[4],x[5]))
                next_researcher += 1
            l += 1
    fig.supxlabel('Pressure [bar]')
    fig.supylabel('$N_2$ Loading [mmol/g]')
    fig.tight_layout()
    fig.savefig(os.path.join(pwd, 'visualise_fits', material + '_N2.png'), facecolor='w')
    plt.close()

    # get DSL model params in a dataframe
    col_names = ["q_s1", "q_s2", "b01", "b02", "deltaH1", "deltaH2", 'Material', 'Forcefield']
    N2_DSL_params = pd.DataFrame(N2_dsl, columns=col_names)
    N2_DSL_params.to_csv(os.path.join(pwd, 'DSL_model_values', material + '_N2_DSLparams.csv'), index=False)