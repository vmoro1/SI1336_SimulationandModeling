#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 10:51:12 2020

@author: viggomoro
"""
import numpy as np
import random
from matplotlib import pyplot as plt
import math


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
    

def random_walk(steps, plot=True):
    """Perform a random walk. It's possible to move up, down, left and right,
    i.e move in the x.y plane along the axis. steps is the number of steps for 
    the random walk and plot is a boolean indicating wether to plot or not"""

    positions = np.zeros((2, steps+1))       # Array to store positions
    for step in range(1, steps+1):
        new_move = next_move()                                  # Generate a new move
        positions[:, step] = positions[:, step - 1] + new_move  # Update positions
    
    if plot:
        # TODO Centrera plotten kring origo
        plt.figure()
        plt.title('A random walk for ' + str(steps) + ' steps')
        plt.xlabel('x')
        plt.ylabel('y')   
        plt.plot(0,marker='o', color='k')
        plt.plot(positions[0, :], positions[1, :])
        plt.plot(positions[0,-1], positions[1,-1], marker='o', color='k')
    else:
        return positions[:, -1]
    
    
class random_number_generator:
    """A random generator"""
    
    def __init__(self, r_0=1, a=3, c=4, m=128):
        self.r_nprev = r_0
        self.a = a
        self.c = c
        self.m = m
        
    def random_number(self):
        """Generate a random number bases on the previous generates number. 
        Generates a rendom number in the range 0,1,..., m-1"""
        r_n = (self.a * self.r_nprev + self.c) % self.m
        self.r_nprev = r_n
        return r_n
         
        
# Kommer fungera bättre då m kan delas jämnt med 4, annars kommer det vara någon som blir vald oftare
def random_walk_rng(steps, r_0=1, a=3, c=4, m=128):
    """Performs a random walk in the same way as random_walk() but using the 
    random generator """
    
    positions = np.zeros((2, steps+1))                # Array to store positions
    number = random_number_generator(r_0, a, c, m)    # Innitiate a random_number_generator object  
    
    for step in range(1, steps+1):
        move = number.random_number()
        move = move // (number.m // 4)  # Convert to int between 0 and 3
        
        if move == 0:
            new_move = [1, 0]
        elif move == 1:
            new_move = [-1, 0]            
        elif move == 2:
            new_move = [0, 1]    
        else:
            new_move = [0, -1]
            
        positions[:, step] = positions[:, step - 1] + new_move  # Update position

    # TODO Centrera plotten kring origo
    plt.figure()
    plt.title('A random walk with the random number generator for ' + str(steps) + '\n steps with r = ' + str(r_0) + ', a = ' + str(a) + ', c = ' + str(c) + ' and m = ' + str(m))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(0,marker='o', color='k')
    plt.plot(positions[0, :], positions[1, :])
    plt.plot(positions[0,-1], positions[1,-1], marker='o', color='k')
    
    
def distance_travelled(number_walks, steps):
    """performs many random walks of a given size. number_wals is the number of 
    walks that are completed. steps is the number of steps the random walk 
    performs. Returns a tuple of expected value of distance and distance squared."""
    
    R = 0       # Used for expected value of distance
    R_2 = 0     # Used for expected value of squared distance
    for walk in range(number_walks):
        final_position = random_walk(steps, False)      # Final position of random walk
        distance = math.sqrt(final_position[0] ** 2 + final_position[1] ** 2)  # Distance travelled
        R += distance
        R_2 += distance ** 2
        
    R = R / number_walks   # Expected value(i.e mean) of the end to end distance
    R_2 = R_2 / number_walks # Expected value(i.e mean) of the squared end to end distance
    
    return R, R_2

def distance_plots(different_steps=np.arange(1,1000, 50).tolist(), plot=True):
    """Plots the below the plots(described bt axis labeling) by repeadly 
    calling distance_travelled. Different steps is a list with the number of 
    steps for the random walk to perform"""
    number_walks = 50                                    # Number of walks for each step   
    
    # Lists to store values
    R_2_list = []
    fluctuation_list = []
    standard_error_list = []
    
    # Update lists for all different step sizes in different_steps
    for steps in different_steps:
        R, R_2 = distance_travelled(number_walks, steps)
        fluctuation = math.sqrt(((R_2 - R ** 2) * number_walks) / (number_walks - 1)) # root-mean-square fluctuation estimate
        standard_error = math.sqrt((R_2 - R ** 2) / (number_walks - 1))               # Standard error estimate  
        
        R_2_list.append(R_2)
        fluctuation_list.append(fluctuation)
        standard_error_list.append(standard_error)
    
    R_2_list = np.array(R_2_list)
    sqrt_R_2_list = np.sqrt(R_2_list).tolist()
    
    if plot:
        plt.figure()
        # plt.title('root mean squared end to end distance versus number of steps')
        plt.ylabel('root mean squared end to end distance')
        plt.xlabel('steps performed by the random walk')
        plt.plot(different_steps, sqrt_R_2_list)
        
        plt.figure()
        # plt.title('root-mean-square fluctuation estimate')
        plt.ylabel('root-mean-square fluctuation estimate')
        plt.xlabel('steps performed by the random walk')
        plt.plot(different_steps, fluctuation_list)
        
        plt.figure()
        plt.ylabel('standard error estimate')
        plt.xlabel('steps performed by the random walk')
        plt.plot(different_steps, standard_error_list) 
        
        plt.show()
    else:
        return sqrt_R_2_list


def self_avoiding_random_walk(steps, returned='plot'):
    """This is a self avoiding random walk. The difference between this and 
    random_walk is that a walk is discarded if it crosses itself. If a walk 
    crosses it self False is returned. Everything else works as the random_walk 
    function with the exception that you also have the possibility to return True 
    if a walk is successful by setting returned='successful"""
    positions = np.zeros((2, steps+1))      # Store previous positions
    
    for step in range(1, steps+1):
        new_move = next_move()                                  # Generate a new move
        positions[:, step] = positions[:, step - 1] + new_move  # Update positions     
        
        # Check if we have already been to the new position
        current_pos = positions[:,step]
        previous_pos = positions[:,:step]
        already_visited = (previous_pos == current_pos.reshape(2,1)).all(axis=0).any() # If we have already been at the position, i.e the walk crosses itself   
        if already_visited:
            successful = False
            return successful
    
    if returned == 'plot':
        # TODO Centrera plotten kring origo
        plt.figure()
        plt.title('A self avoiding random walk for ' + str(steps) + ' steps')
        plt.xlabel('x')
        plt.ylabel('y')   
        plt.plot(positions[0, :], positions[1, :])
    elif returned == 'final_position':
        return positions[:, -1]
    elif returned == 'successful':
        successful = True
        return successful
    
    
def successful_walks(type_walk,title):
    """Plots the fraction of successful random walks versus the number of steps 
    the random wlak takes. A successful walk is a walk that never revisits an 
    old position, i.e crosses itself. type_walk is eitehr self_avoiding_random_walk 
    or self_avoiding_random_walk_improved and title indicates this in the plot."""
    
    steps = np.arange(0, 60, 1).tolist()    # Steps of the random walk
    attempts = 100                          # Attempts to complete a random walk for each stepsize
    fraction_successful_walks = []          # List to store fraction of successful walks for different steps
    
    # Loop over different steps. For each step this will be done the number of attempts 
    for step in steps:
        successful_attempts = 0
        for attempt in range(attempts):
            success = type_walk(step, 'successful')
            if success == True:
                successful_attempts += 1
               
        fraction_successful = successful_attempts / attempts
        fraction_successful_walks.append(fraction_successful)
        
    plt.figure()
    plt.title('Fraction of successfull random walks for the ' + title + ' \nself avoiding random walk')
    plt.xlabel('Steps')
    plt.ylabel('Fraction of successfull walks')
    plt.plot(steps, fraction_successful_walks)        
    

def self_avoiding_random_walk_improved(steps, returned='plot'):
    """This is almost the same function as self_avoiding_random_walk with the 
    exception that it doesn't allow you to perform a step that return you to 
    your previous position. Thus it's better att self avoiding."""
    
    positions = np.zeros((2, steps+1))  # Array with previous positions
    prev_move = [None, None]            # Previous move, initially there is no previous position
    
    for step in range(1, steps+1):
        while True:
            new_move = next_move()      # Generate new move       
            if (-1 * np.array(new_move) == np.array(prev_move)).all(): # If we return to the previous position, stay in the loop and generate a new position
                continue
            else:
                break
        
        # update position and prev_move
        positions[:, step] = positions[:, step - 1] + new_move
        prev_move = new_move
        
        # Check if we have already been to the new position. 
        current_pos = positions[:,step]
        previous_pos = positions[:,:step]
        already_visited = (previous_pos == current_pos.reshape(2,1)).all(axis=0).any() # If we have already been at the position   
        if already_visited:         
            successful = False
            return successful
    
    if returned == 'plot':
        # TODO Centrera plotten kring origo
        plt.figure()
        plt.xlabel('x')
        plt.ylabel('y')   
        plt.plot(positions[0, :], positions[1, :])
    elif returned == 'final_position':
        return positions[:, -1]
    elif returned == 'successful':
        successful = True
        return successful
    


def distance_travelled_self_avoiding(type_walk, number_walks, steps):
    """performs many random walks of a given size. number_wals is the number of 
    walks that are completed. steps is the number of steps the random walk 
    performs. type_walk is eitehr self_avoiding_random_walk or 
    self_avoiding_random_walk_improved. Returns the expected value of the 
    distance squared."""
    
    R_2 = 0                 # Used for expected value of squared distance
    successful_walks = 0    # Count of successful walks
    
    for walk in range(number_walks):
        final_position = type_walk(steps, returned='final_position')    # Final position of random walk      
        
        # If walk is successful
        if type(final_position) != bool:
            distance = math.sqrt(final_position[0] ** 2 + final_position[1] ** 2)  # Distance travelled
            R_2 += distance ** 2
            successful_walks += 1
        
    R_2 = R_2 / successful_walks  # Expected value(i.e mean) of the squared end to end distance  
    return R_2


def distance_plot_self_avoiding(type_walk, different_steps, title, plot=True):
    """x"""
    """Plots the below the plots(described bt axis labeling) by repeadly 
    calling distance_travelled_self_avoiding. type_walk is eitehr self_avoiding_random_walk 
    or self_avoiding_random_walk_improved and title indicates this in the plot. 
    different_steps is a list with the different number of steps the random walk
    will perform."""
    number_walks = 100       # Number of walks for each step   
    R_2_list = []            # Lists to store values
    
    # Update lists for all different step sizes in different_steps
    for steps in different_steps:
        R_2 = distance_travelled_self_avoiding(type_walk, number_walks, steps)      
        R_2_list.append(R_2)
    
    R_2_list = np.array(R_2_list)
    sqrt_R_2_list = np.sqrt(R_2_list).tolist()
    
    if plot:
        plt.figure()
        plt.title('Root mean squared end to end distance for the ' + title + ' \nself avoiding random walk')
        plt.ylabel('Root mean squared end to end distance')
        plt.xlabel('Steps')
        plt.plot(different_steps, sqrt_R_2_list)
    else:
        return sqrt_R_2_list
    
              
# 2.1 a) =====================================================================
# random_walk(10, plot=True)
# random_walk(100, plot=True)
# random_walk(1000, plot=True)

# 2.1 b) =====================================================================
# random_walk_rng(100, r_0=1, a=3, c=4, m=128)
# random_walk_rng(100, r_0=1, a=3, c=4, m=129)
# random_walk_rng(100, r_0=1, a=3, c=4, m=130)

# random_walk_rng(100, r_0=5, a=10, c=15, m=128)
# random_walk_rng(100, r_0=1, a=1, c=2, m=128)
# random_walk_rng(100, r_0=20, a=2, c=5, m=128)
        
# 2.1 c) =====================================================================
# distance_plots()

# 2.1 d) =====================================================================

# Kanske ej ska ha med plotts för själva random walks ty de lyckas nästan aldrig
# self_avoiding_random_walk(2, returned='plot')
# self_avoiding_random_walk(100, returned='plot')
# successful_walks(self_avoiding_random_walk, 'standard')

# # Improved version
# self_avoiding_random_walk_improved(10, returned='plot')
# self_avoiding_random_walk_improved(50, returned='plot')
# successful_walks(self_avoiding_random_walk_improved, 'improved')

# 2.1 e) =====================================================================

# standard self avoiding random walk
# different_steps = np.arange(12)
# distance_plot_self_avoiding(self_avoiding_random_walk, different_steps, 'standard')
    
# sqrt_R_2 = distance_plots(different_steps, plot=False)
# sqrt_R_2_self_avoiding = distance_plot_self_avoiding(self_avoiding_random_walk, different_steps, title='', plot=False)
# plt.figure()
# plt.title('Log-log plot of root mean square end-to-end distances')
# plt.xlabel('Normal random walk')
# plt.ylabel('Standard self avoiding random walk')
# plt.loglog(sqrt_R_2, sqrt_R_2_self_avoiding)

# Improved self avoiding random walk
# different_steps_improved = np.arange(35)
# distance_plot_self_avoiding(self_avoiding_random_walk_improved, different_steps_improved, 'improved')   

# sqrt_R_2 = distance_plots(different_steps_improved, plot=False)
# sqrt_R_2_self_avoiding_improved = distance_plot_self_avoiding(self_avoiding_random_walk_improved, different_steps_improved, title='', plot=False)
# plt.figure()
# plt.title('Log-log plot of root mean square end-to-end distances')
# plt.xlabel('Normal random walk')
# plt.ylabel('Improved self avoiding random walk')
# plt.loglog(sqrt_R_2, sqrt_R_2_self_avoiding_improved)