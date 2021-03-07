# Computes the volume of a 10-dimensional sphere using midpoint integration 
# and the hit and miss monte carlo method

import math as math
import numpy as np
from time import process_time
from matplotlib import pyplot as plt
import random

def discretization(bound, n):
    """Returns a list of n equidistant points between -bound and bound"""
    return np.linspace(-float(bound) + float(bound) / n, float(bound) + float(bound) / n, n, False)


def recursiveIntegral(radius2, volume2, dim, n):
    """Recursively computes an integral of a dim-dimensonal sphere
    with radius sqrt(radius). Pass volume2 = 1 at start.
    n is the number of points used for the midpoint method."""
    volume2 *= radius2 * 2 * 2  / (n * n)
    if (dim > 1):
        partIntegral = 0
        for x in discretization(math.sqrt(radius2), n):
            if (dim == 10):
                print(x)
            partIntegral += recursiveIntegral(radius2 - x * x, volume2, dim - 1, n)
    else:
        partIntegral = math.sqrt(volume2) * n

    return partIntegral


def compute_volume(numPointsPerDim=4, numDims=10):
    """Computes the volume of a 10-dimensional sphere using midpoint integration.
    numPointsPerDim is the number of points in the midpoint method along one dimension.
    numDims is the number of dimensins of the sphere""" 

    t = process_time()
    integral = recursiveIntegral(1, 1, numDims, numPointsPerDim)
    t = process_time() - t
    
    return t, integral


def plot_error(numDims=10):
    """Plots the error between the analytical solution and the midpoint method"""
    analytical = math.pi**(numDims / 2) / math.factorial(numDims / 2)
    time = []
    errors = []
    
    for numPointsPerDim in range(1,6):
        t, integral = compute_volume()
        error = integral - analytical
        errors.append(error)
        time.append(t)
        
    plt.figure()
    plt.title('Midpoint method')
    plt.xlabel('Error')
    plt.ylabel('Computational time')
    plt.plot(errors, time)


def montecarlo(dim, iterations):
    """Monte Carlo approximation of the volume of a dim-dimensonal unit sphere"""
    
    count_in_sphere = 0
    for count_loops in range(iterations):
        point = np.random.uniform(-1.0, 1.0, dim)
        distance = np.linalg.norm(point)
        if distance < 1.0:
            count_in_sphere += 1

    return np.power(2.0, dim) * (count_in_sphere / iterations)

