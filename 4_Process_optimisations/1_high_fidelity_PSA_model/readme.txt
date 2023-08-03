********************************************************************************************************

/CycleSteps/ - functions for each of the individual PSA steps in the modified skarstrom cycle

/GA_files/   - functions to perform the genetic optimisation using NSGA-II algorithm

/opt_results/ - ANN-generated optimisation results for CL1 class of MOFs. The optimised PSA variables 
		obtained from the ANN model is used as an initial seed for the high-fidelity PSA process
		model optimisation. 

/sorted_material_dataframes/ - feature vectors for CRAFTED MOFs for each molecular forcefield

********************************************************************************************************

CL1_MOFs.m	    - list of CL1 class MOFs to be evaluated

ProcessInputParameters.m - Define model parameters such as PSA process parameters, material properties,
			   properties of the flue gas, etc. 

ProcessOptValidation.m   - main script! High fidelity PSA process model refinement. 

PSACycle.m		 - main function which constructs the actual PSA cycle from the individual cycle 
	 		   steps contained in /CycleSteps/

PSACycleSimulation.m     - simple wrapper function for the NSGA-II algorithm. 


********************************************************************************************************
