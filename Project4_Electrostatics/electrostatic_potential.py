#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 07:18:48 2020

@author: viggomoro
"""

import numpy as np
import matplotlib.pyplot as plt
import random

import matplotlib
matplotlib.rcParams.update({'font.size': 12})

def laplace_approx(n=10, nsteps=10, checker=1, initial_grid=0, sequential=False):
    """Determines the potential in a square region with the relaxation method. 
    The grid is n+1 points along x and y, including boundary points 0 and n. 
    nsteps is the number of iterations. checker indicates order of update, 
    checker=1 is regular and 2 is like a chess board. Initial_grid is the 
    guess of the potential in the interior and is the value which is set for 
    all interior positions. sequential indicates if the sequential relaxation
    methid should be used. """
    
    # Initialize potentiqal of the grid to initial_grid
    v = np.zeros((n+1, n+1))    # Matrix with potentials
    vnew = np.zeros((n+1, n+1))
    v += initial_grid
    vnew += initial_grid
    
    # Set the boundary conditions
    for i in range(n+1):
        v[0,i] = 10
        v[n,i] = 10
        v[i,0] = 10
        v[i,n] = 10
    
    # Successive approximations for all positions
    for i in range(nsteps):
        if not sequential:
            v, vnew = relax(n, v, vnew, checker)
        else:
            v, vnew = relax_sequentially(n, v, vnew, checker)
        
    return v


def relax(n, v, vnew, checker):
    """Perform one step of relaxation. The relaxation method determines the 
    potential as the mean of the four nearby potentials"""
    
    # Update potentials in vnew
    for check in range(0,checker):  
        for x in range(1,n):
            for y in range(1,n):
                if (x*(n+1) + y) % checker == check:
                    vnew[x,y] = (v[x-1][y] + v[x+1][y] + v[x][y-1] + v[x][y+1])*0.25

        # Copy back the new values to v
        for x in range(1,n):
            for y in range(1,n):
                if (x*(n+1) + y) % checker == check:
                    v[x,y] = vnew[x,y]
    return v, vnew


def relax_sequentially(n, v, vnew, checker):
    """Perform one step of sequentialrelaxation. The difference between this 
    and relax is that the an updated potential is directly stored in v and thus
    used to compute sequential potentials. The relaxation method determines the 
    potential as the mean of the four nearby potentials"""
    
    # Update potentials in v
    for check in range(0,checker):  
        for x in range(1,n):
            for y in range(1,n):
                if (x*(n+1) + y) % checker == check:
                    v[x,y] = (v[x-1][y] + v[x+1][y] + v[x][y-1] + v[x][y+1])*0.25

    return v, vnew
    


def accuracy(n=10, checker=1, initial_grid=9, plot=True):
    """Plots accuracy versus number of steps for laplace_approx. If plot=false 
    the smallets number of steps that produces an accuract < 0.01 is returned"""
    
    
    nsteps = np.arange(1, 130, 1)  # different number of steps
    accuracy_list = []             # store accuracy
    guess = initial_grid
    grid_size = n
    
    # Determine accuracy for att different steps
    for step in nsteps:
        v = laplace_approx(n=grid_size, nsteps=step, checker=1, initial_grid=guess)
        minimum = np.amin(v)
        error = 10 - minimum
        accuracy = error / 10
        accuracy_list.append(accuracy)
        
        if not plot and accuracy < 0.01:
            return step
        
    plt.figure()
    plt.title('Accuracy versus number of steps for grid size n = ' + str(n))
    plt.xlabel('Number of steps')
    plt.ylabel('Accuracy')
    plt.plot(nsteps, accuracy_list)
    plt.plot([nsteps[0], nsteps[-1]], [0.01, 0.01])
    plt.legend(['Accuracy', '0.01'])


def countour_plot(n=10, nsteps=100, checker=1, initial_grid=5, title=''):
    """Plots the equipotential surfaces of the solution of laplace_approx"""
    grid_size = n
    guess = initial_grid
    steps = nsteps
    
    v = laplace_approx(n=grid_size,nsteps=steps, checker=checker, initial_grid=guess)
    x = np.arange(n+1)
    y = np.arange(n+1)
    
    plt.figure()
    plt.title('Equipotential surfaces for boundaries ' + title)
    plt.contour(x,y,v)    
    
 
    
def compare_accuracy(n=10, checker=1, initial_grid=9, plot=True):
    """Plots accuracy versus number of steps for laplace_approx for both the 
    sequential and non-sequential relaxation method."""
    
    nsteps = np.arange(1, 100, 1)
    accuracy_list = []
    accuracy_list_sequential = []
    guess = initial_grid
    grid_size = n
    
    for step in nsteps:
        v = laplace_approx(n=grid_size, nsteps=step, checker=1, initial_grid=guess)
        minimum = np.amin(v)
        error = 10 - minimum
        accuracy = error / 10
        accuracy_list.append(accuracy)
        
        v_sequential = laplace_approx(n=grid_size, nsteps=step, checker=1, initial_grid=guess, sequential=True)
        minimum_sequential = np.amin(v_sequential)
        error_sequential = 10 - minimum_sequential
        accuracy_sequential = error_sequential / 10
        accuracy_list_sequential.append(accuracy_sequential)
        
    plt.figure()
    plt.title('Comparison of accuracy between sequential and \nnon-sequential relaxation for grid size n = ' + str(n))
    plt.xlabel('Number of steps')
    plt.ylabel('Accuracy')
    plt.plot(nsteps, accuracy_list)
    plt.plot(nsteps, accuracy_list_sequential)
    plt.plot([nsteps[0], nsteps[-1]], [0.01, 0.01])
    plt.legend(['Non-sequential relaxation', 'Sequential relaxation', '0.01'])
    
    
def compare_checkering(n=10, checker=1, initial_grid=9, plot=True):
    """Plots accuracy versus number of steps for the two different orders in 
    which the potential can be computed, i.e value of checker."""
    
    nsteps = np.arange(1, 100, 1)
    accuracy_list = []
    accuracy_list_checkured = []
    guess = initial_grid
    grid_size = n
    
    for step in nsteps:
        v = laplace_approx(n=grid_size, nsteps=step, checker=1, initial_grid=guess)
        minimum = np.amin(v)
        error = 10 - minimum
        accuracy = error / 10
        accuracy_list.append(accuracy)
        
        v_checkured = laplace_approx(n=grid_size, nsteps=step, checker=2, initial_grid=guess)
        minimum_checkured = np.amin(v_checkured)
        error_checkured = 10 - minimum_checkured
        accuracy_checkured = error_checkured / 10
        accuracy_list_checkured.append(accuracy_checkured)
        
    plt.figure()
    plt.title('Comparison of accuracy between the two different orders of \nupdating the potenital for grid size n = ' + str(n))
    plt.xlabel('Number of steps')
    plt.ylabel('Accuracy')
    plt.plot(nsteps, accuracy_list)
    plt.plot(nsteps, accuracy_list_checkured)
    plt.plot([nsteps[0], nsteps[-1]], [0.01, 0.01])
    plt.legend(['checker = 1', 'checker = 2', '0.01'])
    
    
def random_walk_approx(n=10, nwalks=100, start_position=[5,5]):
    """Computers the potential at positions start_position with a random walk 
    approach. nwalks is the number of walks that are performed to determine the potential."""
   
    v = np.zeros((n+1, n+1))    # Matrix with boundary potentials
    
    # Set the boundary conditions
    for i in range(n+1):
        v[0,i] = 10
        v[n,i] = 10
        v[i,0] = 5
        v[i,n] = 5
        
    # Determine potential for each walk and take the average
    potential = 0 
    for i in range(nwalks):
        position = start_position[:]
        # Perform one random walk and store the potential
        while (position[0] != 0 and position[0] != n and position[1] != 0 and position[1] != n):
            new_move = next_move()
            position[0] = position[0] + new_move[0]
            position[1] = position[1] + new_move[1]
        potential += v[position[1], position[0]]
    
    potential /= nwalks
    return potential


def next_move():
    """Generate a new move out of the four possibilities"""
    move = int(4 * random.random())
    if move == 0:
        return [1, 0]
    elif move == 1:
        return [-1, 0]            
    elif move == 2:
        return [0, 1]    
    else:
        return [0, -1]
    
    
def comparison_rw_and_relaxation(n=10, initial_grid=7, start_position=[5,5]):
    """Comparison between accuracy of the relaxation method and the random wlk method for 
    the position start position and with the initial guess for relaxation equal
    to initial_grid."""
    
    nsteps = np.arange(10, 103, 3)

    exact_value = laplace_approx(n=n, nsteps=10000, checker=1, initial_grid=initial_grid)[start_position[1], start_position[0]]  # approxiate exact value with value for many steps 
    # List to store results
    accuracy_relaxation_list = []
    accuracy_rw_list = []
    
    for step in nsteps:
        # Accuracy for relaxation method
        v_relaxation = laplace_approx(n=n, nsteps=step, checker=1, initial_grid=initial_grid)
        error = exact_value - v_relaxation[start_position[1], start_position[0]]
        accuracy_relaxation = error / exact_value
        accuracy_relaxation_list.append(accuracy_relaxation)
        
        # Accuracy for random walk method, mean of 1000 results are used
        accuracy_rw = 0
        iters = 1000 
        for i in range(iters):
            v_rw = random_walk_approx(n=n, nwalks=step, start_position=start_position)
            error = abs(exact_value - v_rw)
            accuracy_rw += error / exact_value
        accuracy_rw /= iters      
        accuracy_rw_list.append(accuracy_rw)
        
    plt.figure()
    plt.title('Comparison between acuracy of relaxation method with an initial \nguess of 7 everywhere and the random walk method. Each \nrandom walk data point is the mean of ' + str(iters) + ' random walks')
    plt.xlabel('Number of steps')
    plt.ylabel('Accuracy at position ' + str(start_position))
    plt.plot(nsteps, accuracy_relaxation_list)
    plt.plot(nsteps, accuracy_rw_list)
    plt.legend(['Relaxation method', 'Random walk method'])


def comparison_rw_relaxtion_fix_step(n=10, nsteps=100, initial_grid=7.5, start_position=[5,5]):
    """Comparison between the relaxation method and the random wlk method for 
    the position start position and with the initial guess for relaxation equal
    to initial_grid. The difference compared to comparison_rw_and_relaxation is 
    that the comparison now is for a fixed step."""
    
    number_simulations = np.arange(1,101)   # Number of simulations
    rw_values = []                          # List to store different values for the potential in each iteration
    
    for simulation in number_simulations:
        v_rw = random_walk_approx(n=n, nwalks=nsteps, start_position=start_position)                                               # Potenital according to random walk method
        rw_values.append(v_rw)
    
    laplace_value = laplace_approx(n=n, nsteps=nsteps, checker=1, initial_grid=initial_grid)[start_position[1], start_position[0]] # Potential according to relaxation method
    # exact_value = laplace_approx(n=n, nsteps=10000, checker=1, initial_grid=initial_grid)[start_position[1], start_position[0]]    # Potenital for many steps are considered to be the exact potential
        
    plt.figure()
    plt.title(str(len(number_simulations)) + ' comparisons between the relaxation method with an initial guess of\n ' + str(initial_grid) + ' and the random walk method. The number of steps/walks are ' + str(nsteps))
    plt.xlabel('Simulation')
    plt.ylabel('Value of potential at position ' + str(start_position))
    plt.plot([number_simulations[0], number_simulations[-1]], [laplace_value, laplace_value])
    # plt.plot([number_simulations[0], number_simulations[-1]], [exact_value, exact_value])
    plt.scatter(number_simulations, rw_values, c='orange')
    plt.legend(['Relaxation method', 'Random walk method'])

 
# 10.10 a) ===================================================================

# accuracy(n=10, checker=1, initial_grid=9, plot=True)
# num_steps = accuracy(n=10, checker=1, initial_grid=9, plot=False)    
# accuracy(n=5, checker=1, initial_grid=9, plot=True)    
# num_steps = accuracy(n=5, checker=1, initial_grid=9, plot=False)   
    
# 10.10 b) ===================================================================
    
# Ändra dessutom potentialen när den innitieras i laplace_approx så att den är 
# 4 i centrum. Kolla också animeringen för att kolla på skillnader i utvecklingen
# accuracy(n=10, checker=1, initial_grid=0, plot=True)
# num_steps = accuracy(n=10, checker=1, initial_grid=0, plot=False)    
# accuracy(n=5, checker=1, initial_grid=0, plot=True)    
# num_steps = accuracy(n=5, checker=1, initial_grid=0, plot=False)   
    
# 10.10 c) ===================================================================

# Ändra hur v innitieras i laplace_approx enligt instruktion 
# countour_plot(n=10, nsteps=100000, checker=1, initial_grid=7.5, title='5, 10, 5 and 10')
# countour_plot(n=10, nsteps=100000, checker=1, initial_grid=7.5, title='\n10 on three sides and 0 on the fourth')

# 10.11 a) ===================================================================
    
# compare_accuracy(n=10, checker=1, initial_grid=9, plot=True)   

# 10.11 b) ===================================================================
    
# compare_checkering(n=10, checker=1, initial_grid=9, plot=True)  

# 10.17 a) ===================================================================  
# Kom ihåg att ändra i laplace_approx till rätta boundary value

# comparison_rw_relaxtion_fix_step(n=10, nsteps=100, initial_grid=7.5, start_position=[5,5])
# comparison_rw_relaxtion_fix_step(n=10, nsteps=1000, initial_grid=7.5, start_position=[5,5])

# 10.17 b) ===================================================================  
# Kom ihåg att ändra i laplace_approx till rätta boundary value

# comparison_rw_relaxtion_fix_step(n=10, nsteps=100, initial_grid=7.5, start_position=[1,1])
# comparison_rw_relaxtion_fix_step(n=10, nsteps=1000, initial_grid=7.5, start_position=[1,1])

# comparison_rw_relaxtion_fix_step(n=10, nsteps=100, initial_grid=7.5, start_position=[3,3])
# comparison_rw_relaxtion_fix_step(n=10, nsteps=1000, initial_grid=7.5, start_position=[3,3])

# comparison_rw_and_relaxation(n=10, initial_grid=7, start_position=[5,5])
# comparison_rw_and_relaxation(n=10, initial_grid=7, start_position=[1,1])
# comparison_rw_and_relaxation(n=10, initial_grid=7, start_position=[3,3])
