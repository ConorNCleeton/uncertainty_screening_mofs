********************************************************************************************************

/1_Latin_hypercube_sampling/ - Code to generate training data for a surrogate model using a latin hyper-
                               cube sampling of the input phase space (optimisable PSA process variables 
                               and MOF material properties) -> main script = PSA_LHS.m

********************************************************************************************************

/2_Bootstrap_sampling/ - Code used to generate high purity training data points from the 
                         PSA process model. Bootstrap approach used whereby different random seeds 
                         are used to instantiate the genetic algorithm optimisation (described fully 
                         in the SI, section 4.2) -> main script = BootstrapOptimization.m

********************************************************************************************************

/3_training_data/ - training data generated from LHS and Bootstrap optimisation. Both the transformed 
                    data (i.e., testDataEcon_transformed.txt) and raw data (i.e., total_raw_data.mat)
                    provided. Data transformations are used to train the ANN model, and are described 
                    in the SI, section 4.2

********************************************************************************************************