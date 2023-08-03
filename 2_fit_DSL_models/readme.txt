********************************************************************************************************

CO2_df.csv - contains the CO2 GCMC simulated data for CRAFTED MOFs across all forcefields and temperatures

N2_df.csv  - contains the N2 GCMC simulated data for CRAFTED MOFs across all forcefields and temperatures

main.py - python script which runs the DSL model fitting subroutine in an automated fashion. The fittings 
          are then printed to a subdirectory called /visualise_fits/ (user must create this directory). 
          The fitted DSL model parameters are then printed to a .csv file in /DSL_model_values/ subdirectory 
          (user must create this directory).

fit_DSL.py - function to perform nonlinear regression of DSL model parameters against GCMC simulated isotherms

********************************************************************************************************