import scipy.optimize as opt
import operator
import random
import pandas as pd
import numpy as np
from scipy.stats import linregress
from lmfit import minimize, Parameters
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def CO2_DSL_regress(df,init_guess=[]):
    """
    Function to determine the initial parameters for a DSL model based on the experimental or GCMC data passed in.
    It is assumed that the maximum loading value is close to saturation.
    
    Inputs:
        - df          = dataframe containing the data to regress against
        - init_guess  = initial guess for model parameters (dictionary, optional)
    Outputs: 
        - result = Object containing the optimized parameters and several goodness-of-fit statistics. 
    """
    col_headers = df.columns.values.tolist()
    P_key = col_headers[0]
    L_key = col_headers[1]
    T_key = col_headers[2]
    ID_key = col_headers[3]
    error_key =  col_headers[4]
    name_key = col_headers[5]
    
    

    # ************************************* INITIAL GUESS FOR PARAMS ****************************************************
    lowT_df = df[df[T_key]==273]
    total_saturation = lowT_df[L_key].max()*1.05       # assume saturation value is the pressure point at 10bar 

    ## saturation capacity parameters
    q_s1_guess = total_saturation*0.5        # initial guess for q_s1 and q_s2 
    q_s2_guess = q_s1_guess 

    ## Heats of adsorption of each site
    # First determine K_H from the linear portion of the isotherm
    nonzero = lowT_df[lowT_df[L_key] != 0]                       # get all non-zero loading values
    delta_P = nonzero[P_key].nsmallest(2).values       # smallest values for P and L (to determine an approximation of K_H)
    delta_L = nonzero[L_key].nsmallest(2).values   

    if len(delta_P) < 2 or len(delta_L) < 2:
        K_H = 0
    else:
        K_H = (delta_L[1]-delta_L[0])/(delta_P[1]-delta_P[0])   # approximation to K_H based on smallest two loading values

    # Create an initial guess for the heat of adsorption based on the henry coefficient, defining appropriate 
    # upper and lower bounds for K_H to use for interpolating values of H1 
    UpperK_H = 400; logUpperK_H = np.log(UpperK_H)
    LowerK_H = 0.0001; logLowerK_H = np.log(LowerK_H) 
    interpK_H = min(max(LowerK_H,K_H),UpperK_H); loginterpK_H = np.log(interpK_H)
    H1_guess = 17000 +((interpK_H-LowerK_H)*((44000-17000)/(UpperK_H-LowerK_H)))        
#     H1_guess = 17000 +((loginterpK_H-logLowerK_H)*((44000-17000)/(logUpperK_H-logLowerK_H)))        
    H2_guess = 0.9*H1_guess

    b01_guess = K_H*0.5/(np.exp(H1_guess/(8.314*lowT_df[T_key].iloc[0]))*q_s1_guess)
    b02_guess = K_H*0.5/(np.exp(H2_guess/(8.314*lowT_df[T_key].iloc[0]))*q_s2_guess)
    
    # create a Parameters() object
    params = Parameters()
    try:
        params.add('qs1',value=init_guess['qs1'],min=0,max=12.5)
    except:
        params.add('qs1',value=q_s1_guess,min=0,max=12.5)
    try:
        params.add('qs2',value=init_guess['qs2'],min=0,max=12.5)
    except:
        params.add('qs2',value=q_s2_guess,min=0,max=12.5)
    try:
        params.add('b02',value=init_guess['b02'],min=1e-9,max=1e-2,vary=True)
        params.add('con2',value=init_guess['b01']-init_guess['b02'],min=0,vary=True)
    except:
        params.add('b02',value=b02_guess,min=1e-9,max=1e-2)
        params.add('con2',value=b01_guess-b02_guess,min=0,vary=True)
    try:
        params.add('deltaH2',value=init_guess['deltaH2'],min=15000,max=45000,vary=True)
        params.add('con1',value=init_guess['deltaH1']-init_guess['deltaH2'],min=0, vary=True)
    except:
        params.add('deltaH2',value=H2_guess,min=15000,max=45000,vary=True)
        params.add('con1',value=H1_guess-H2_guess,min=0, vary=True)
    params.add('deltaH1',expr='deltaH2+con1',min=15000,max=45000)
    params.add('b01',expr='b02+con2',min=1e-9,max=1e-2)
    
    # **************************************** FUNCTIONS TO CURVE FIT ***************************************************
    def residual(params,P,T,data):
        qs1     = params['qs1']
        qs2     = params['qs2']
        b01     = params['b01']
        b02     = params['b02']
        deltaH1 = params['deltaH1']
        deltaH2 = params['deltaH2']
        site_1 = qs1 * (P*b01 * np.exp(deltaH1/(8.314*T)))/(1+ (P*b01*np.exp(deltaH1/(8.314*T))))
        site_2 = qs2 * (P*b02 * np.exp(deltaH2/(8.314*T)))/(1+ (P*b02*np.exp(deltaH2/(8.314*T))))
        
        model = site_1 + site_2
        return (data-model)    
    
    # perform minimisation
    result = minimize(residual,
                      params,
                      args=(df[P_key],df[T_key],df[L_key]))

    
    # ********************************************* GET PARAMETERS FROM REGRESSION **********************************
    try_opt_count = 0  
    while not result.success:       # if regression not successful, check the reason and try again 
        
        try_opt_count += 1; 
    
        # randomly change the initial guess
        op_mappings = {"+":operator.add,
                       "-":operator.sub}
        op = random.choice(["+", "-"])
        
        for param in range(4):
            # Only change the b0 and deltaH parameters
            param_dict = result.params.valuesdict()
            guess = [param_dict['b01'],param_dict['b02'],param_dict['deltaH1'],param_dict['deltaH2']]
            guess[param] = op_mappings[op](guess[param], guess[param]*random.uniform(0,0.1))
        
        params.add('b02',value=guess[1],min=1e-9,max=1e-2)
        params.add('deltaH2',value=guess[3],min=15000,max=45000)
        params.add('con1',value=guess[2]-guess[3],min=0, vary=True)
        params.add('con2',value=guess[0]-guess[2],min=0,vary=True)   

        result = minimize(residual,
                          params,
                          args=(df[P_key],df[T_key],df[L_key]))
        print('Number of attempts at fitting model to data: {}'.format(try_opt_count))        
        if try_opt_count > 5:
            break


    # ***************************************** RETURN OPTIMISED INITIAL PARAMETERS *********************************
    return result


def CO2_DSL_regress_BACKUP(df,tol=1e-12,norm_factor=1,maxiter=100,disp=False,init_guess=[]):
    """
    Function to determine the initial parameters for a DSL model based on the experimental or GCMC data passed in.
    It is assumed that the maximum loading value is close to saturation. 
    This is a backup function which uses scipy's optimise function rather than lmfit to perform the regression 
    in case lmfit fails.
    
    Inputs:
        - df          = dataframe containing the data of an individual researcher at a single temperature
        - tol         = regression precision tolerance
        - norm_factor = normalisation factor which determines how accurately the RSS minimisation is (lower value 
            means a less accurate regression). Convergence is more difficult with higher values, and so this factor can be 
            reduced to improve any convergence issues at the expense of regression accuracy. 
        - maxiter     = number of function evaluations for regression. 
        - init_guess  = initial guess for model parameters
    Outputs: 
        - opt_res = optimised DSL model parameters for initial guess estimation
    """
    col_headers = df.columns.values.tolist()
    P_key = col_headers[0]
    L_key = col_headers[1]
    T_key = col_headers[2]
    ID_key = col_headers[3]
    error_key =  col_headers[4]
    name_key = col_headers[5]

    # ************************************* BOUNDS FOR OPTIMISATION PROBLEM *********************************************
    # define approximate bounds for params
    qs1_b = (0, 12.5)
    qs2_b = (0, 12.5)
    b01_b = (1e-9, 1e-2)
    b02_b = (1e-9, 1e-2)
    deltaH1_b = (15000, 45000)
    deltaH2_b = (15000, 45000)
    
    bnds = (qs1_b,qs2_b,b01_b,b02_b,deltaH1_b,deltaH2_b)  
    cons = ({'type': 'ineq', 'fun': lambda x: x[2] - x[3]},    # inequality const b01 - b02 >= 0   (i.e. b01 >= b02)
            {'type': 'ineq', 'fun': lambda x: x[4] - x[5]})    # inequality const deltaH1 - deltaH2 >= 0 (i.e. deltaH1 >= deltaH2)
    
    # ************************************* INITIAL GUESS FOR PARAMS ****************************************************
    if init_guess == []:
        lowT_df = df[df[T_key]==273]
        total_saturation = lowT_df[L_key].max()*1.05       # assume saturation value is the pressure point at 10bar 
        
        ## saturation capacity parameters
        q_s1_guess = total_saturation*0.5        # initial guess for q_s1 and q_s2 
        q_s2_guess = q_s1_guess 
        
        ## Heats of adsorption of each site
        # First determine K_H from the linear portion of the isotherm
        nonzero = lowT_df[lowT_df[L_key] != 0]                       # get all non-zero loading values
        delta_P = nonzero[P_key].nsmallest(2).values       # smallest values for P and L (to determine an approximation of K_H)
        delta_L = nonzero[L_key].nsmallest(2).values   
        
        if len(delta_P) < 2 or len(delta_L) < 2:
            K_H = 0
        else:
            K_H = (delta_L[1]-delta_L[0])/(delta_P[1]-delta_P[0])   # approximation to K_H based on smallest two loading values

        # Create an initial guess for the heat of adsorption based on the henry coefficient, defining appropriate 
        # upper and lower bounds for K_H to use for interpolating values of H1 
        UpperK_H = 400; logUpperK_H = np.log(UpperK_H)
        LowerK_H = 0.0001; logLowerK_H = np.log(LowerK_H) 
        interpK_H = min(max(LowerK_H,K_H),UpperK_H); loginterpK_H = np.log(interpK_H)
        H1_guess = 17500 +((interpK_H-LowerK_H)*((43500-17500)/(UpperK_H-LowerK_H)))        
        H2_guess = 0.9*H1_guess
        # H1_guess = np.random.uniform(17000,45000)
        # H2_guess = 0.9*H1_guess                     # create H2 = H1

        ## pre-exponential factor
        # assuming K_H1 = 0.7*K_H and K_H2 = 0.3*K_H (first site adsorbs stronger) we can determine b01 and b02 using the 
        # formula K_H1 = q_s1*b01*exp(H1_guess/R*T)
       
        b01_guess = K_H*0.5/(np.exp(H1_guess/(8.314*lowT_df[T_key].iloc[0]))*q_s1_guess)
        b02_guess = K_H*0.5/(np.exp(H2_guess/(8.314*lowT_df[T_key].iloc[0]))*q_s2_guess)
        guess = np.array([q_s1_guess,q_s2_guess,b01_guess,b02_guess,H1_guess,H2_guess])           # initial guess array 
    
    elif len(init_guess) == 6:
        guess  = init_guess
    else:
        raise Exception('''You have not supplied an initial guess vector for the CO2 DSL model which is incompatible with the 
                            the current function. Please omit submitting a parameter vector (in which case an estimate 
                            of the model parameters will be made), or submit a parameter vector that contains an initial guess 
                            for [q_s1,q_s2,b01,b02,deltaH1,deltaH2] ''')    

    
    # **************************************** FUNCTIONS TO CURVE FIT ***************************************************
    def DSL_loading(P,T):
        '''
        Function to determine the DSL loading as a function of the DSL model parameters guessed at pressure P.
        - guess[0] = saturation loading of site 1, q_s1
        - guess[1] = saturation loading of site 2, q_s2
        - guess[2] = pre-exponential factor for site 1, b01
        - guess[3] = pre-exponential factor for site 2, b02
        - guess[4] = Heat of adsorption for site 1, deltaH1
        - guess[5] = Heat of adsorption for site 2, deltaH2
        '''
        site_1 = guess[0] * (P*guess[2] * np.exp(guess[4]/(8.314*T)))/(1+ (P*guess[2]*np.exp(guess[4]/(8.314*T))))
        site_2 = guess[1] * (P*guess[3] * np.exp(guess[5]/(8.314*T)))/(1+ (P*guess[3]*np.exp(guess[5]/(8.314*T))))
        return site_1 + site_2

    def curve_fit(): 
        ''' 
        Function to perform nonlinear regression against the reference data that is passed into the parent function
        initial_params_guess
        '''
        
        def RSS(x):
            ''' 
            residual sum of squares to minimise the loss function for curve fitting
            '''
            for i in range(len(guess)):
                guess[i] = x[i]
                residual = np.sum((df[L_key].values-DSL_loading(df[P_key],df[T_key]).values)**2)/norm_factor                
            return residual 

        
        # minimise the RSS
        opt_res = opt.minimize(RSS,guess,
                               method='SLSQP',
                               bounds = bnds,
                               constraints=cons,
                               options={'maxiter': maxiter,
                                        'ftol': tol,
                                        'disp': disp})

        return opt_res
    
    # ********************************************* GET PARAMETERS FROM REGRESSION **********************************
    
    opt_res = curve_fit()            # perform the regression
    try_opt_count = 0  
    while not opt_res.success:       # if regression not successful, check the reason and try again 
        
        try_opt_count += 1; 
        failure_message = opt_res.message
        if failure_message == 'Positive directional derivative for linesearch':
            tol = min(1e-10,tol*10)  # tolerance is >= 1e-9, initially we start with 1e-12, but we can increase the 
        #                             # tolerance if linedirection not found. Don't go below 1e-9 though. 
            
            print(failure_message)
        
        if failure_message == 'Maximum number of function evaluations has been exceeded.':
            max_iter += 25        # Increase the number of function evaluations
    
        # randomly change the initial guess
        op_mappings = {"+":operator.add,
                       "-":operator.sub}
        op = random.choice(["+", "-"])
        
        for param in range(4):
            # Only change the b0 and deltaH parameters
            guess[param+2] = op_mappings[op](guess[param+2], guess[param+2]*random.uniform(0,0.1))
        
        opt_res = curve_fit()            # perform the regression again
        print('Number of attempts at fitting model to data: {}'.format(try_opt_count))
        failure_message = opt_res.message
        
        if try_opt_count > 10:
            break
    
#     if not opt_res.success:
#         print(opt_res.message)
#         raise Exception("""Minimization of RSS failed.""")

    # ***************************************** RETURN OPTIMISED INITIAL PARAMETERS *********************************
    return opt_res



# determine the initial model parameters for the DSL model based on characteristics of the data passed in
def N2_DSL_regress(df,CO2_params,N2_nonlin_tol=0.85,tol=1e-15,norm_factor=1,maxiter=100,disp=False,init_guess=[]):

    """
    Function to determine the initial parameters for a DSL model based on the experimental or GCMC data passed in.
    It is assumed that the maximum loading value is close to saturation.
    
    Inputs:
        - df            = dataframe containing the data of an individual researcher at a single temperature
        - CO2_params    = CO2_DSL model params which are used to determine the saturation capacity of N2
        - N2_nonlin_tol = R2 (coefficient of determination) tolerance to determine if N2 is flagged as linear or nonlinear. 
                          For example, 0.985 means that 98.5% of the variance in the N2 isotherm data is captured by a linear 
                          model, therefore any fit with R2 > 0.985 will be considered 'linear' for the N2 regression
        - tol           = regression precision tolerance
        - norm_factor   = normalisation factor which determines how accurately the RSS minimisation is (lower value 
            means a less accurate regression). Convergence is more difficult with higher values, and so this factor can be 
            reduced to improve any convergence issues at the expense of regression accuracy. 
        - maxiter       = number of function evaluations for regression. 
        - init_guess    = initial guess for model parameters
    Outputs: 
        - opt_res.x = optimised DSL model parameters for initial guess estimation
    """
    
    ## Keys
    col_headers = df.columns.values.tolist()
    P_key = col_headers[0]
    L_key = col_headers[1]
    T_key = col_headers[2]
    ID_key = col_headers[3]
    error_key =  col_headers[4]
    name_key = col_headers[5]
    
    # ************************************* BOUNDS FOR OPTIMISATION PROBLEM *********************************************
    # define approximate bounds for params
    b01_b = (1e-9, 1e-2)
    b02_b = (1e-9, 1e-2)
    deltaH1_b = (5000, 25000)
    deltaH2_b = (5000, 25000)
    bnds = (b01_b,b02_b,deltaH1_b,deltaH2_b) 
    
    def check_nonlinearity():
        # Function to check if the N2 isotherm exhibits significant nonlinearity
        lowT = pd.unique(df[T_key]).min()
        low_iso = df[df[T_key]==lowT]
        m, c, r_value, p_value, std_err = linregress(low_iso[P_key], low_iso[L_key])
        return r_value 
    r2 = check_nonlinearity()
     
        
    if r2 > N2_nonlin_tol:
        N2_nonlinear = 0
    else:
        N2_nonlinear = 1

    if N2_nonlinear == 1:
        cons = ({'type': 'ineq', 'fun': lambda x: x[0] - x[1]},  # inequality const b01 - b02 >= 0   (i.e. b01 >= b02)
                {'type': 'ineq', 'fun': lambda x: x[2] - x[3]})  # inequality const deltaH1 - deltaH2 = 0 (i.e. deltaH1 >= deltaH2)
    elif N2_nonlinear == 0:
        cons = ({'type': 'eq', 'fun': lambda x: x[0] - x[1]},    # equality const b01 - b02 = 0   (i.e. b01 = b02)
                {'type': 'eq', 'fun': lambda x: x[2] - x[3]})    # equality const deltaH1 - deltaH2 = 0 (i.e. deltaH1 = deltaH2)
    else: 
        raise Exception('''Nonlinear flag options are 1 (for nonlinear N2 adsorption isotherms) or 0 (for linear N2 adsorption
        isotherms).''')

    # ************************************* INITIAL GUESS FOR PARAMS ***************************************************

    if init_guess == []: 
        q_s1 = CO2_params['q_s1'].iloc[0]
        q_s2 = CO2_params['q_s2'].iloc[0]
        ## Heats of adsorption of each site
        # First determine K_H from the linear portion of the isotherm
        nonzero = df[df[L_key] != 0]                       # get all non-zero loading values       
        delta_P = nonzero[P_key].nsmallest(2).values       # smallest values for P and L (to determine an approximation of K_H)
        delta_L = nonzero[L_key].nsmallest(2).values   

        if len(delta_P) < 2 or len(delta_L) < 2: 
            K_H  = 0
        else:       
            K_H = (delta_L[1]-delta_L[0])/(delta_P[1]-delta_P[0])   # approximation to K_H based on smallest two loading values
            
        # Create an initial guess for the heat of adsorption based on the henry coefficient, defining appropriate 
        # upper and lower bounds for K_H to use for interpolating values of H1 
        UpperK_H = 25; logUpperK_H = np.log(UpperK_H)
        LowerK_H = 0.00001; logLowerK_H = np.log(LowerK_H) 
        interpK_H = min(max(LowerK_H,K_H),UpperK_H); loginterpK_H = np.log(interpK_H)
        H1_guess = 6000 +((interpK_H-LowerK_H)*((24000-6000)/(UpperK_H-LowerK_H)))        
        H2_guess = H1_guess
        
        # H1_guess = np.random.uniform(5000,25000)
        # H2_guess = H1_guess                     # create H2 = H1

        b01_guess = K_H*0.5/(np.exp(H1_guess/(8.314*df[T_key].iloc[0]))*q_s1)
        b02_guess = b01_guess
        guess = np.array([b01_guess,b02_guess,H1_guess,H2_guess])           # initial guess array 
        
    else: 
        if len(init_guess) == 4: # if only supplying b01,b02,deltaH1,deltaH2, then we need to define q_s1 and q_s2 from CO2_params
            guess = init_guess
            q_s1 = CO2_params['q_s1'].iloc[0]
            q_s2 = CO2_params['q_s2'].iloc[0]
        elif len(init_guess) == 6: 
            guess = init_guess[2:]
            q_s1 = init_guess[0]
            q_s2 = init_guess[1]
        else:
            raise Exception('''You have supplied an initial guess vector for the N2 DSL model which is incompatible with the 
                            the current function. Please omit submitting a parameter vector (in which case an estimate 
                            of the model parameters will be made), submit a parameter vector that contains an initial 
                            guess for [b01,b02,deltaH1,deltaH2], or submit a parameter vector that contains an initial guess 
                            for [q_s1,q_s2,b01,b02,deltaH1,deltaH2] ''')    

        
    # **************************************** FUNCTIONS TO CURVE FIT ***************************************************
    def DSL_loading(P,T):
        '''
        Function to determine the DSL loading as a function of the DSL model parameters guessed at pressure P.
        - guess[0] = pre-exponential factor for site 1, b01
        - guess[1] = pre-exponential factor for site 2, b02
        - guess[2] = Heat of adsorption for site 1, deltaH1
        - guess[3] = Heat of adsorption for site 2, deltaH2
        '''
        site_1 = q_s1 * (P*guess[0] * np.exp(guess[2]/(8.314*T)))/(1+ (P*guess[0]*np.exp(guess[2]/(8.314*T))))
        site_2 = q_s2 * (P*guess[1] * np.exp(guess[3]/(8.314*T)))/(1+ (P*guess[1]*np.exp(guess[3]/(8.314*T))))
        return site_1 + site_2
    
    def curve_fit(): 
        ''' 
        Function to perform nonlinear regression against the reference data that is passed into the parent function
        initial_params_guess()
        '''
        def RSS(x):
            ''' residual sum of squares to minimise the loss function for curve fitting
            - fun : callable; The objective function to be minimized.
            - fun(x, *args) -> float`` where ``x`` is an 1-D array with shape (n,) and ``args``
              is a tuple of the fixed parameters needed to completely
              specify the function.
            '''
            for i in range(len(guess)):
                guess[i] = x[i]
                residual = np.sum((df[L_key].values-DSL_loading(df[P_key],df[T_key]).values)**2)/norm_factor # divide by norm_factor to make opt easier
            return residual 

        # minimise the RSS
        opt_res = opt.minimize(RSS,guess,
                               method='SLSQP',
                               bounds = bnds,
                               constraints=cons,
                               options={'maxiter': maxiter,
                                        'ftol': tol,
                                        'disp': disp})
        return opt_res
    
    # ********************************************* GET PARAMETERS FROM REGRESSION **********************************
    
    opt_res = curve_fit()            # perform the regression
    try_opt_count = 0  
    while not opt_res.success:       # if regression not successful, check the reason and try again 
        
        try_opt_count += 1; 
        failure_message = opt_res.message
        
        if failure_message == 'Positive directional derivative for linesearch':
#             print('old tol = {}'.format(tol))
            tol = min(1e-12,tol*10)  # tolerance is >= 1e-12, initially we start with 1e-12, but we can increase the 
                                    # tolerance if linedirection not found. Don't go below 1e-9 though. 
            
#             print('new tol = {}'.format(tol))
            print(failure_message)
            
        if failure_message == 'Maximum number of function evaluations has been exceeded.':
            max_iter += 25        # Increase the number of function evaluations
    
        # randomly change the initial guess
        op_mappings = {"+":operator.add,
                       "-":operator.sub}
        op = random.choice(["+", "-"])        
        
        for param in range(len(guess)):
            guess[param] = op_mappings[op](guess[param], guess[param]*random.uniform(0,0.1))
        
       
        opt_res = curve_fit()            # perform the regression again
        print('Number of attempts at fitting model to data: {}'.format(try_opt_count))
        failure_message = opt_res.message
        
        if try_opt_count > 10:
            break
    
#     if not opt_res.success:
#         print(opt_res.message)
#         raise Exception("""Minimization of RSS failed.""")

    # ***************************************** RETURN OPTIMISED INITIAL PARAMETERS *********************************
    return opt_res
