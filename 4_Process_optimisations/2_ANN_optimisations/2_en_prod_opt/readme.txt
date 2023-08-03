********************************************************************************************************

/GA_files/   - functions to perform the genetic optimisation using NSGA-II algorithm

/sorted_material_dataframes/ - feature vectors for CRAFTED MOFs for each molecular forcefield

********************************************************************************************************

net_En_logb0.m      - Energy penalty ANN model 

net_Prod_logb0.m    - Productivity ANN model 

net_Pu_logb0.m      - CO2 purity ANN model 

net_Re_logb0.m 	    - CO2 recovery ANN model 

objective_evaluation.m - function which calculates the CO2 purity and recovery KPIs for NSGA-II optimisa-
 			 tion. 

run_visualisations.m   - script to visualise all of the ANN-generated results (after the optimisations 
			 have been completed) and save to .png files. 


runANNEconOpt.m  - main script. Running this script will optimise the energy penalty - productivity performance 
                      of a material (or all CRAFTED materials, if desired), extract the final Pareto front 
		      as well as the optimised PSA process variables, and save the visualised results in 
		      as a png file.  

visualise_results.m - helper function to visualise the optimisation results for each MOF as they are being 
	   	      completed. 

********************************************************************************************************
