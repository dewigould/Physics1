#Set to True or False to calculate results displayed in report
see_results = False
validation_routine = True

import copy
import numpy as np

#Function: To determine the machine precision for a given floating point scheme
#Input: Number of bits of floating point scheme (32,64,128)
#Output: Calculated value for machine precision of that variable type
def precision(num_bits):
    
    #Map argument of function (number of bits) to the corresponding argument required for numpy.dtype
    mapping = {32:'Float32', 64:'Float64', 128:'Float128'}
    
    if num_bits not in mapping.keys():
        print "ERROR: Function not able to cope with this number of bits"
        return 0
    #Store values 1.0 and 2.0 in the specific floating point scheme required
    information = np.array([1.0,2.0], dtype=mapping[num_bits])
    
    #avoid any issues with losing digits/ accumulating floating point errors by copying memory location
    one = copy.deepcopy(information[0])
    two = copy.deepcopy(information[1])
    
    #epislon is given an initial guess: 1.0
    #epsilon will never use many significant bits - so precision is not a problem
    epsilon = copy.deepcopy(information[0]) 
    
    #Algorithm subtracts epsilon guess from one and checks if the result is distinct from one
    #Epislon is repeatedly updated with epislon/2.0 until convergence.
    while one-epsilon != one:
        epsilon/=two
        
    print "Floating Point Type: ",num_bits," Floating Point Precision: ", epsilon
 
    
if see_results == True:
    
    #Calculate single, double and extended precision accuracies
    precision(32)
    precision(64)
    precision(128)

    #Calculate machine accuracy (default precision if 64-bit)
    precision(64)
    
if validation_routine == True:
    
    #Try input invalid number of bits - should raise error
    precision(10)
    
    #Validated with default machine precision 64-bit
    precision(64)
