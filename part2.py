see_results = False
validation_routine = True

if see_results == True:
    #input required matrice and vector
    A = [[2,1,0,0,0],[3,8,4,0,0],[0,9,20,10,0],[0,0,22,51,-25],[0,0,0,-55,60]]
    b = [2,5,-4,8,9]
        
import numpy as np

#Function: Perform LU Decomposition for any arbritrary NxN Matrix, A
#Input: NxN Matrix 
#OutPut: L - lower diagonal matrix
#        U - upper diagonal matrix
#        new_matrix - (matrix of all elements of U and non-diagonal elements of L)
def LU_Decomp(A):

    #Number of rows and columns
    N = len(A)   
    
    #set-up three empty NxN matrices to be filled by algorithm
    L = np.zeros((N,N))
    U = np.zeros((N,N))
    new_matrix = np.zeros((N,N))    
    
    #assign arbritrary values (=1) to N parameters to reduce problem to N^2 equations in N^2 unknowns 
    for i in range(N):
        L[i][i] = 1.0               

    #Implement Crout's Algorithm, values of U[i][j] and L[i][j] are calculated one step before they are needed
    for j in range(N):
        for i in range(N):    
        
            if i <= j:
                summation_one = 0
                for k in range(i):
                    summation_one += L[i][k]*U[k][j]
                U[i][j] = A[i][j] - summation_one
                new_matrix[i][j] = U[i][j] 
                
            if i > j:
                summation_two = 0
                for k in np.arange(j):
                    summation_two += L[i][k]*U[k][j]
                L[i][j] = (1.0/(U[j][j])) * ( A[i][j] - summation_two)
                new_matrix[i][j] = L[i][j]  
   
    return L,U,new_matrix



#Function: Calculate Determinant of Matrix A, given upper diagonal matrix U obtained from LU Decomposition
#Input: Upper diagonal matrix U of A
#Output: Determinant of A
def LU_determinant(U):
    
    #Give determinant start value of 1.0 to apply multiplication to
    det_A = 1.0
    #Number of rows or columns
    N = len(U)          
    
    #Calculate by multiplying diagonal elements of U
    for i in range(N):
        det_A *= U[i][i]
        
    return det_A



#Function: Solve matrix equation Ax=b for x, in form LUx=b
#Input: L,U and b, two matrices formed by LU Decomposition and vector b
#Output: solution to equation - vector x
def solve_LU_equation(L,U,b):
    
    #Number of rows or columns
    N = len(U)
    
    #Create empty vectors to be filled by algorithm
    x = [0]*N
    y = [0]*N
    
    #Perform Forward Substitution to populate y
    y[0] = b[0]/L[0][0]
    for i in np.arange(1,N):
        summation_one = 0.0
        for k in range(i):
            summation_one += L[i][k]*y[k]
        y[i] = (1.0/L[i][i])*(b[i] - summation_one)
        
    #Perform Backward Subsitution to populate solution x, using vector y
    x[N-1] = y[N-1]/ (U[N-1][N-1])
    for i in range(N)[::-1]:
        summation_two = 0.0
        for k in np.arange(i+1,N):
            summation_two += U[i][k]*x[k]    
        x[i] = (1.0/U[i][i])*(y[i] - summation_two)
        
    return x


#Function: Find Inverse of matrix A, using LU Decomposition techniques
#Input: L, U - obtained by LU Decomposition of matrix A
#Output: Inverse of A (in matrix form)
def LU_method_inverse(L,U):
    
    #Number of rows or columns
    N = len(L)
    
    #Create NxN identity matrix
    I = np.zeros((N,N)) 
    for i in range(N):
        I[i][i] = 1
    
    #Set-up empty NxN matrix to be filled by algorithm
    inv_A = np.zeros((N,N))
    
    #Perform algorithm for calculating inverse using ith column of identity matrix (ie. ith basis vector in R^n)
    for i in range(N):
        e_i = [j[i] for j in I]             
        j_col = solve_LU_equation(L,U,e_i)
        
        #Transpose rows and columns to obtained ith row of Inverse of A
        for element in j_col:    
            inv_A[j_col.index(element)][i] = element    
            
    return inv_A


if see_results == True:
    #Decompose matrix A into L and U
    L,U,R = LU_Decomp(A)

    #Calculate determinant of A
    det = LU_determinant(U)

    #Solve matrix equation Ax=b for x
    x = solve_LU_equation(L,U,b)

    #Print out Results
    print "A: ", A
    print "A decomposed by LU Decomposition into: "
    print "L: ", L
    print "U: ", U
    print "Matrix containing U and non-diagonal L: ", R
    print "Solution to equation Ax=b: ", x
    print "Determinant of A: ", det
    print "Inverse of A: ", LU_method_inverse(L,U)
    

if validation_routine == True:
    
    A = [[1,0,0],[0,1,0],[0,0,1]]
    b = [1,1,1]
    
    L,U,R = LU_Decomp(A)
    if np.array_equal(L,A) != True or LU_determinant(U) != 1.0:
        print "Function no longer working with trivial example of identity matrix"
    else:
        print "Validated for identity matrix"



