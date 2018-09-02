see_results = False
validation_routine = False

if see_results == True:
    #File Location - .csv file containing data to be interpolated
    file_name = "Interpolation_data_csv.csv"

    #Define cubic spline increments (number of plotted points for each interval in x data)
    cubic_spline_increments = 100

import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt

#Import functions to do matrix methods from part2 script
import part2 as matrix_methods


#Function: Takes csv file of N points of 2 variables (one independent and one dependent) and creates Nx2 array called "data"
#Input: .csv file extension
#Output: array of x - independent data points, 
#        array of y - dependent data points
def process_data(file_extension):

    #Create pandas dataframe from information contained in .csv file
    df = pd.read_csv(file_extension)
    
    #Number of data points
    N_data = len(df)
    
    #Total number of variables in problem (dependent or independent)
    N_dimensions = len(df.columns)
    
    #generate empty array to store data
    data = np.zeros((N_data,N_dimensions))

    #populate data array with information in pandas data frame
    for i in df.itertuples():
        for j in range(N_dimensions):
            data[i[0]][j] = i[j+1]
    
    #Split into independent and dependent variables
    x = [i[0] for i in data]
    y = [i[1] for i in data]
    
    return x,y
            
#Function: Perform Linear interpolation on input data
#Input: x_set - array of 'x' data
#       y_set - array of 'y' data
#       increment - user defined number of increments over which interpolation performed
#Output: x_to_plot,y_to_plot: arrays containing x,y points to be plotted (forming the linear interpolation line)
def linear_interpolate(x_set,y_set):
    
    #Number of data points
    N_data = len(x_set)
    
    #Create empty arrays which will be populated by interpolant by algorithm
    x_to_plot = []
    y_to_plot = []

    #Loop over all data points and create interpolant for each. 
    #Not including the final data point in this loop as it has no point after it to "connect it to"
    for i in range(N_data-1):   
        
        #Define interpolant function to "join the dots"
        def f(x):
            result = ((x_set[i+1]-x)*y_set[i] + (x - x_set[i])*y_set[i+1])/ (x_set[i+1]-x_set[i])    
            return result

        #Split up x range in-between two points in questions, in even intervals
        x_range = np.linspace(x_set[i],x_set[i+1],3)
        
        #Calculate f(x) for each of these intermediate values
        y_values = [f(x) for x in x_range]
        
        #Append interpolant data to arrays to be plotted
        for j in range(3):
            x_to_plot.append(x_range[j])
            y_to_plot.append(y_values[j])
               
    return x_to_plot,y_to_plot
        

#Function: Perform Cubic Spline operation to yield fit curve
#Input: x_set - array of 'x' data
#       y_set - array of 'y' data
#       increment - user defined number of increments over which interpolation performed
#Output: x_to_plot,
#        y_to_plot - arrays containing x,y points to be plotted (forming the cubic spline)
def cubic_spline(x_set,y_set,increment):
    
    #Number of data points
    N_data = len(x_set)
    
    #Create empty arrays which will be populated by interpolant by algorithm
    x_to_plot = []
    y_to_plot = []
    
    #Create empty (N-2)x(N-2) array and vector of N-2 dimensions, A,b to be filled by algorithm.
    A = np.zeros((N_data-2,N_data-2))
    b = [0]*(N_data-2)
        
    #Fill the matrix A with coefficients determined by cubic spline boundary conditions 
    for i, j in itertools.product(np.arange(2,N_data),repeat=2):
        #Populate diagonal entries
        if i == j:
            A[i-2][i-2] = 2*(x_set[i]-x_set[i-2])
        #Populate Diagonal +1 entries
        if j == i+1:
            A[i-2][j-2] = x_set[i]-x_set[i-1]
        #Populate Diagonal -1 entries
        if j == i-1:
            A[i-2][j-2] = x_set[i-1]-x_set[i-2]
        #Populate entries in vector b
        if b[i-2] == 0.0:
            b[i-2] = 6* (((y_set[i]-y_set[i-1])/(x_set[i]-x_set[i-1])) - ((y_set[i-1]-y_set[i-2])/(x_set[i-1]-x_set[i-2])))

    #Use LU-Decompsition to split matrix of coefficients A into L and U
    L,U,R = matrix_methods.LU_Decomp(A)
    
    #Solve matrix equation Ax=b, where x is vector of second derivatives.
    #Manually add the two natural spline boundary conditions as 0.0 entries at beginning and end of vector f_double_primes
    f_double_primes = [0.0] + matrix_methods.solve_LU_equation(L,U,b)
    f_double_primes.append(0.0)
    
    #Find cubic function f(x) from data points and second derivatives calculated above
    for i in range(N_data-1):
        def A(x):
            result = (x_set[i+1]-x)/(x_set[i+1]-x_set[i])
            return result
        def B(x):
            result = 1.0 - A(x)
            return result
        def C(x):
            result = (1.0/6.0)*(((A(x))**3)-A(x))*((x_set[i+1]-x_set[i])**2)
            return result
        def D(x):
            result = (1.0/6.0)*(((B(x))**3)-B(x))*((x_set[i+1]-x_set[i])**2)
            return result
        #Function to Calculate Cubic Spline at any x value within range defined by data set.
        def f(x):
            result = (A(x)*y_set[i]) + (B(x)*y_set[i+1]) + (C(x)*f_double_primes[i]) + (D(x)*f_double_primes[i+1])
            return result
        
        #Generate evenly spaced x_range within data points in questions - for which the cubic spline is evaluated
        x_range = np.linspace(x_set[i],x_set[i+1],increment)
        y_values = [f(x) for x in x_range]
        
        #Append these x,y values of cubic spline to arrays for plotting.
        for j in range(increment):
            x_to_plot.append(x_range[j])
            y_to_plot.append(y_values[j])
            
    return x_to_plot,y_to_plot


if see_results == True:
    #Process data from .csv to arrays
    x,y = process_data(file_name)
    #Perform linear interpolation
    interpolate_x,interpolate_y = linear_interpolate(x,y)
    #Perform Cubic Spline Operation
    cubic_spline_x,cubic_spline_y = cubic_spline(x,y,cubic_spline_increments)

    #Plot Results
    plt.title("Interpolation Routine")
    plt.scatter(x,y,c="red",marker="x",label = "Data Points")
    plt.plot(interpolate_x,interpolate_y,c = "blue", label = "Linear Interpolation")
    plt.plot(cubic_spline_x,cubic_spline_y,c="green",label= "Cubic Spline")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()
    
if validation_routine == True:
    
    #Validate for simple linear function
    
    interpolate_x,interpolate_y = linear_interpolate(np.linspace(0,1000,10),[n+5 for n in np.linspace(0,1000,10)])
    cubic_spline_x,cubic_spline_y = cubic_spline(np.linspace(0,1000,10),[n+5 for n in np.linspace(0,1000,10)],1000)
    

    #Plot Results
    plt.title("Interpolation Routine")
    plt.scatter(np.linspace(0,1000,10),[n+5 for n in np.linspace(0,1000,10)],c="red",marker="x",label = "Data Points")
    plt.plot(interpolate_x,interpolate_y,c = "blue", label = "Linear Interpolation")
    plt.plot(cubic_spline_x,cubic_spline_y,c="green",label= "Cubic Spline")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()


