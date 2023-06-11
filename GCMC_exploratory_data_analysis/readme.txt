********************************************************************************************************

/atom_charge_csv/ - contains the csv files that compares the charge predictions of Qeq, EQeq, and DDEC 
                    for different atom types in the CRAFTED database

/CIF_FILES/ - contains the cif files with charges pre-computed for CRAFTED

/GCMC_data/ - contains two subdirectories: (1) MIXTURE_ISOTHERMS, and (2) SINGLE_ISOTHERMS. In each dir, 
              separate csv files are used to store the GCMC simulated data of different CRAFTED MOF 
              structures. ID codes are used to differentiate between different forcefields

              - 1 = UFF+DDEC
              - 2 = UFF+EQeq      
              - 3 = UFF+Qeq
              - 4 = UFF+Neutral
              - 5 = Dre+DDEC
              - 6 = Dre+EQeq
              - 7 = Dre+Qeq
              - 8 = Dre+Neutral

********************************************************************************************************

Exploratory_data_analysis.ipynb - jupyter notebook to visualise adsorption data and calculate material-
                                  level correlations between forcefields

structural_data.csv - contains LCD, PLD, void fraction, etc., calculated using zeo++ for CRAFTED MOFs.

neverDOE.txt - MOFs which never meet DOE constraints (CO2 purity and recovery >= 0.9) for any combination
               of the forcefield parameter (CL4 class, as described in the main article section 3.4)

structurally_inconsistent_MOFs.txt - MOFs with structural abnormalities (described in the main article, 
                                     section 2.1.)