#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 14:43:32 2020

@author: viggomoro

"""
import numpy as np
import random
from matplotlib import pyplot as plt


def next_move():
    """Generate a new move out of the four possibilities in a random walk"""
    move = int(4 * random.random())
    if move == 0:
        return [1, 0]
    elif move == 1:
        return [-1, 0]            
    elif move == 2:
        return [0, 1]    
    else:
        return [0, -1]


def Greens_function_approxRW(n=10, nwalks=200, start_position=[5,5]):
    """Approximates Greens function for start_position by a ranom walk approach. 
    nwalks is the number of walsk used for this approximation and n is the grid
    size."""
   
    G_func = np.zeros((n+1, n+1))    # Matrix to indicate how often the walk end att each boundary position
    
    for i in range(nwalks):
        position = start_position[:]
        # Perform one random walk and store the boundary position it reaches
        while (position[0] != 0 and position[0] != n and position[1] != 0 and position[1] != n):
            new_move = next_move()
            position[0] = position[0] + new_move[0]
            position[1] = position[1] + new_move[1]
        G_func[position[1], position[0]] += 1
    
    return G_func
    

def potential_from_Greens(boundary_potential, n=10, G=Greens_function_approxRW(n=10, nwalks=200, start_position=[5,5]), nwalks=200):
    """Computes the potential at a the position start_position. G is a Greens 
    function computed by Greens_function_approxRW. The potential at 
    start_position is returned."""
    
    v = 0           # potential
    G /= nwalks     # Scale Greens function to get probabilities
    
    # Iterate over all boundary positions and add their contributions to the potential
    for i in range(1,n):
        v1 = boundary_potential[0,i] * G[0,i]
        v2 = boundary_potential[n,i] * G[n,i]
        v3 = boundary_potential[i,0] * G[i,0]
        v4 = boundary_potential[i,n] * G[i,n]
        v += v1 + v2 + v3 + v4

    return v 
    
def find_potential(n=10,nwalks=200):
    """Calculates the potential from the Greens function at every position and 
    produces a plot of the potential. """
    
    boundary = np.zeros((n+1, n+1))   # Store boundary potentials
    v = np.zeros((n+1, n+1))          # Store potential for all positions
    
    # Set the boundary conditions
    for i in range(1,n):
        boundary[0,i] = 10
        boundary[n,i] = 10
        boundary[i,0] = 5
        boundary[i,n] = 5
    
    boundary[3,0] = boundary[4,0] = boundary[5,0] = boundary[6,0] = boundary[7,0] = 20   # Set the the boundary position that maximizes the potential at [3, 5] to 20
    v = np.copy(boundary)      # Store potential for all positions

    # Compute Greens function for each position
    for x in range(1,n):
        for y in range(1,n):
            position = [x, y]       # Position to compute Greens function for
            Greens_func = Greens_function_approxRW(n=n, nwalks=nwalks, start_position=position)  # The Greens function
            
            # Find potential at current position
            v_pos = potential_from_Greens(boundary, n=n, G=Greens_func, nwalks=nwalks)
            v[position[1], position[0]] = v_pos
    
    # v is now computed for all locations and can be plotted
    fig = plt.figure()
    plt.title('Potential for all positions computed using the Greens function for\n every position. Boundary values that maximizes the potentials\n at position [3,5] are used.')
    ax = fig.add_subplot(111)
    im = ax.imshow(v, cmap=None, interpolation='nearest')
    fig.colorbar(im)
    
    
def plot_Greens_function(n=10, nwalks=500, position=[5,5]):
    """Plots the Greens function for the position position."""

    G = Greens_function_approxRW(n=n, nwalks=nwalks, start_position=position)
    G /= nwalks
       
    fig = plt.figure()
    plt.title('Greens function for position ' + str(position) + '. Every position indicates the\n fraction of the random walks that ended up in that position.')
    ax = fig.add_subplot(111)
    im = ax.imshow(G, cmap=None, interpolation='nearest')
    fig.colorbar(im)
    

# 10.18 a) ===================================================================

# plot_Greens_function(n=10, nwalks=100000, position=[5,5])
# plot_Greens_function(n=10, nwalks=100000, position=[3,5])
# plot_Greens_function(n=10, nwalks=100000, position=[5,3])
  
# 10.18 b) ===================================================================  
# Ändra de positioner med potential 20 för att maximera potentialen för olika positioner. 
# och ändra titeln i plotten till den positione
# find_potential(n=10, nwalks=1000)
