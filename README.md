# Physics1
Various computational methods for physics 

#——————————
- part1.py
- part2.py
- part3.py
- part4.py
- part5.py


In each case, mini validation routines have been programmed into the script - just set validation_routine = True at the top to see these run

#——————————
- part1.py:

To obtain results in report, set variable:

```
see_results = True
```
The script is capable of calculating precisions accuracy for 32, 64 and 128 bit floating-point schemes.

#——————————

- part2.py
To obtain results in report, set variable:
```
see_results = True
```
This script performs all elements asked in Part (b).

A user can input an NxN matrix A of their choice, and an N dimensional vector b of their choice, at the top of the script.

The script is pre-programmed to print:
	- The result of the LU Decomposition
	- A matrix containing all elements of U, and all non-diagonal elements on L
	- The solution to the equation Ax=b
	- The determinant of the matrix A
	- The inverse of the matrix A
	
#——————————

- part3.py
To obtain results, set variable:
```
see_results = True
```
This part takes the tabulated data in the form of a .csv file - this could easily be modified is the data is of another form in the future. 
The user can input what .csv file they want interpolated, and set the number of increments over which the interpolation is performed.

The script is pre-programmed to show:
	- Scatter Plot of Data Set
	- Linear Interpolating Function 
	- Cubic Spline Interpolating Function with number of increments set by user (100 found to be fine enough).

#——————————

- part4.py
To obtain results, set variable:
```
see_results = True
```
The user can define all of the quantities they desire at the start of the script: Period,Number of samples, delta_t…
The script is also designed to work for any input functions, they just have to be changed by the user.

The script is pre-programmed to plot the required results, plus some extra graphs for interest:
	- The two functions to be convolved
	- The result of the convolution
	- The DFT of both functions individually
	- The result of the convolution compared to scipy.convolve for robustness

#——————————

- part5.py
To obtain results required in To obtain results in report, set variable:
```
see_results = True,
numerical_solution = False
```
—— this variable get be set = True, to see numerical solutions to rejection method without analytically calculated CDF inverse (takes longer to run)

At the top of the script, the user can define all the parameters required - number of samples, bins for plotting..etc.
The user can also change what pdf’s they want to transform the uniform deviate into.

There is also a section where the comparison function (for rejection method), must be input
If the user requires the results to be found quickest, it is also necessary for the inverse of the CDF of the comparison function to be entered
If not, the program is set up to perform this calculation numerically, as long as the correct arguments are passed to the get_result_rej function.

The script is pre-programmed to plot several things:

- The uniform deviate of 100,000 samples as required in part (a)
- A graph showing the two pdf’s, a constant comparison function, and a better-conforming comparison function
	
- The result of a transformation method to 0.5sin(x)
- The result of a (numerical and slow) transformation method to (2/pi)sin^2(x)

- The result of a rejection method with constant comparison function for (2/pi)sin^2(x)
- The result of a rejection method wth better-conforming comparison function for (2/pi)sin^2(x)


The script is already set-up to print out all the times for each process, and the number of operations each process takes for analysis in the report.





