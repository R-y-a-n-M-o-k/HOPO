import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_coeffs_powers(filename):
    """Obtain the linear coefficients of DMS
       obtains associated powers, torsion argument for given coefficients
       filename is the name of the res file containing fit data"""
    
    x_params, y_params, z_params = [], [], []
    x_powers, y_powers, z_powers = [], [], []

    with open(filename, 'r') as file:
        #skip unnecessary lines 
        line = file.readline()
        while line[0] != ' ':
            line = file.readline()
    
        line_split = line.split()
        while line_split[0][0] == 'd':
            
            #get powers, convert to float 
            powers = line_split[2:-2]
            for j in range(len(powers)):
                powers[j] = float(powers[j])
            
            param = float(line_split[-1])
            
            #append into correct dipole component 
            if line_split[0][1] == 'x':
                x_powers.append(powers)
                x_params.append(param)
            elif line_split[0][1] == 'y':
                y_powers.append(powers)
                y_params.append(param)
            else:
                z_powers.append(powers)
                z_params.append(param)
                
            line = file.readline()
            line_split = line.split()

    return (np.array(x_params), np.array(y_params), np.array(z_params),
            np.array(x_powers), np.array(y_powers), np.array(z_powers))

def get_non_linear_p(filename):
    """Obtain the non-linear coefficients of DMS fit from 
       given residuals file, filename."""
    
    non_linear_params = []
    with open(filename, 'r') as file:
        #skip unnecessary lines 
        line = file.readline()
        while line[0] == '*':
            line = file.readline()
    
        #get nonlin_params 
        while line[0] != ' ':
            line_split = line.split()
            non_lin_param = float(line_split[-1])
            non_linear_params.append(non_lin_param)
            
            line = file.readline()
    return non_linear_params

def get_dms(linear_params, nonlin_eq, powers, coordinate):
    """computes the dms function for a given component direction,
       using input parameters (linear params, nonlin_eq) and powers
       coordinate a string eg. 'x', 'y', 'z'
    """
    def dms(geometries):
        """Fitted dms function (R^6 -> R)
           - Assumes geometries is a n by 6 ndarray of input values, 
           - Assumes angles are given in degrees, (function converts to radians)
           Outputs value of electric dipole moment for those given geometries.
           """
        dipoles = np.zeros(np.shape(geometries)[0])
        
        #set up variable arg to be modified to calculate dipole moment values
        #start with input geometries
        arg = geometries.copy()
        arg[:,:5] = arg[:,:5] - nonlin_eq #geometry relative to equilbrium
        arg[:,4] = -arg[:,4]
        arg[:,3:] = arg[:,3:] * np.pi/180 #convert deg to rad 
        
        #compute one sum-product term per loop
        for j in range(len(powers[:,0])):
            #initialise product parts
            term = np.zeros(np.shape(arg)[0])
            term.fill(linear_params[j])
            
            products = arg[:].copy()
            #compute polynomial, fourier step
            products[:,:5] = products[:,:5]**powers[j,:-1]
            if coordinate == 'y':       
                products[:,5] = (np.sin(products[:,5])*
                                 np.cos(powers[j,-1]*products[:,5]) )
            else:
                products[:,5] = np.cos(powers[j,-1]*products[:,5])
            
            #compute term
            for k in range(6):
                term *= products[:,k]
            
            dipoles += term
            
        return dipoles
            
    return dms