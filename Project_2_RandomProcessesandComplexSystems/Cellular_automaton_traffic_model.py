#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 13:10:55 2020

@author: viggomoro
"""

import random
import matplotlib.pyplot as plt
import numpy as np

        
def traffic_simulation(road_length=100, number_cars=10, timesteps=100, v_max=2, p=0.5):
    """Simulates traffic according to the cellular automation traffic modell. 
    Two matrices with entry (i,j) corresponding to the position/velocity of car j at 
    timestep i are returned. road_length is length of road, numbers_cars is the 
    number of cars in the road and timesteps is the number of steps the simulation 
    runs. v_max is maximum velocity and p is probability of braking."""
    # v_max = 2   # Maximum velocity
    # p = 0.5     # Probability of braking
    
    # Matrices to store position and velocities for each timestep, entry (i,j) is at timestep i for car j
    positions = np.zeros((timesteps + 1, number_cars))
    velocities = np.zeros((timesteps + 1, number_cars))
    
    # Initial conditions
    current_pos = []
    current_vel = []
    for i in range(number_cars):
        current_pos.append(int(i % road_length))   # Append position and make apply periodicity  
        current_vel.append(1)
    
    positions[0,:] = current_pos[:]
    velocities[0,:] = current_vel[:]
    
    # Simulate one timestep
    for t in range(1, timesteps + 1):       
        # accelerate on step towards max speed
        for i in range(number_cars):
            if current_vel[i] < v_max:
                current_vel[i] += 1
        # brake if close to next car
        for i in range(number_cars):
            distance = (current_pos[(i+1)%number_cars] - current_pos[i]) % road_length  # Distance to car in front
            if current_vel[i] >= distance:
                current_vel[i] = distance - 1
        # random braking
        for i in range(number_cars):
            if random.random() <= p and current_vel[i] >= 1 :
                current_vel[i] -= 1
        # move cars forwards, the road is a periodic i.e a circle
        for i in range(number_cars):
            current_pos[i] += current_vel[i]
            current_pos[i] %= road_length     # Creates periodicity
            
        # Store position and velocity
        positions[t,:] = current_pos[:]
        velocities[t,:] = current_vel[:]
    
    return positions, velocities   
    

def plot_simulation(road_length=100, number_cars=10, timesteps=100, title=''):
    """Plots time versus the position of the cars as a scatter plot where each 
    color is a individual car."""
    time = list(range(timesteps + 1))
    position, vel = traffic_simulation(road_length, number_cars, timesteps)
    
    plt.figure()
    plt.title(title)
    for car in range(number_cars):
        plt.xlabel('Position')
        plt.ylabel('Time')
        plt.scatter(position[:,car], time, s=3)
    
    
def plot_flow_rate(road_length=100, timesteps=100, number_simulations=10, different_number_cars=np.arange(2,20), title=''):
    """Plots flow rate versus car density. The flow rate is defined as the sum 
    of the cars verlocities at the end at equilibrium. The car density is 
    defined as the number of cars divided by road length."""
    flow_rates = []
    average_flow_rates = []
    density = []
    
    # Calculate flow rate and density
    for number_cars in different_number_cars:
        # Do it multiple times for and take average
        for simulatoin in range(number_simulations):
            positions, velocities = traffic_simulation(road_length, number_cars, timesteps)
            final_velocities = velocities[timesteps,:]    
            flow_rate = sum(final_velocities) / road_length                    
            flow_rates.append(flow_rate)
        
        average_flow_rates.append(np.mean(flow_rates))
        density.append(number_cars / road_length)
        flow_rates = []
        
    plt.figure()
    plt.title(title)
    plt.xlabel('Density')
    plt.ylabel('Flow rate')
    plt.plot(density, average_flow_rates)
    
    
def standard_error_flow_rate(number_simulations, road_length=50, number_cars=25, timesteps=100):
    """Computes the standard error of the flow rate. number simulation is the 
    number of sample of the flow rate that is generated."""
    flow_rates = []
    for i in range(number_simulations):
        positions, velocities = traffic_simulation(road_length, number_cars, timesteps)
        final_velocities = velocities[timesteps,:]    
        flow_rate = sum(final_velocities) / road_length
        flow_rates.append(flow_rate)
     
    #Calculate the standard error of the flow rate    
    variance = np.var(flow_rates)
    standard_error = np.sqrt(variance / (number_simulations - 1))
    return standard_error


def plot_standard_error_flow_rate():
    """Plots the standard error of the flowrite versus the size of the sample, 
    i.e number of simulation."""
    
    different_number_simulations = np.arange(100, 3100, 100)
    standard_errors = []
    for number_simulations in different_number_simulations:
        standard_error = standard_error_flow_rate(number_simulations)
        standard_errors.append(standard_error)
    
    plt.figure()
    plt.title('Standard error for 25 cars and a road length of 50 meters.')
    plt.xlabel('Number of simulations')
    plt.ylabel('Standard error')
    plt.plot(different_number_simulations, standard_errors)


def plot_flow_rates_different_lengths(timesteps=100, number_simulations=10):
    """Plots flow rate versus car density for different road lengths in the 
    same plot. The flow rate is defined as the sum of the cars verlocities at 
    the end at equilibrium. The car density is defined as the number of cars 
    divided by road length."""
    road_lengths = [18, 25, 50, 70, 90]
    different_number_cars_list = [list(range(2,16)), list(range(2,22)), list(range(2,45,2)), list(range(2,60, 3)), list(range(2,80, 3))]
    
    plt.figure()
    plt.title('')
    plt.xlabel('Car density')
    plt.ylabel('Flow rate')
    
    # Determine plot for each road length
    for i in range(len(road_lengths)):
        road_length = road_lengths[i]
        different_number_cars = different_number_cars_list[i]
        
        # Store values
        flow_rates = []
        average_flow_rates = []
        density = []
       
        # Calculate flow rate and density
        for number_cars in different_number_cars:
            # Do it multiple times for and take average
            for simulatoin in range(number_simulations):
                positions, velocities = traffic_simulation(road_length, number_cars, timesteps)
                final_velocities = velocities[timesteps,:]    
                flow_rate = sum(final_velocities) / road_length                    
                flow_rates.append(flow_rate)
            
            average_flow_rates.append(np.mean(flow_rates))
            density.append(number_cars / road_length)
            flow_rates = []
        plt.plot(density, average_flow_rates)
        flow_rates = []
        density = []
    for i in range(len(road_lengths)):
        road_lengths[i] = 'L = ' + str(road_lengths[i])
    plt.legend(road_lengths)
    
def plot_flow_rates_different_vmax(timesteps=100, number_simulations=10):
    """Plots flow rate versus car density for different v_max in the 
    same plot. The flow rate is defined as the sum of the cars verlocities at 
    the end at equilibrium. The car density is defined as the number of cars 
    divided by road length."""
    different_vmax = [1, 2, 5]
    different_number_cars = np.arange(2,45, 2)
    road_length = 50
    
    plt.figure()
    plt.title('')
    plt.xlabel('Car density')
    plt.ylabel('Flow rate')
    
    # Determine plot for v_max
    for vmax in different_vmax:
        
        # Store values
        flow_rates = []
        average_flow_rates = []
        density = []
       
        # Calculate flow rate and density
        for number_cars in different_number_cars:
            # Do it multiple times for and take average
            for simulatoin in range(number_simulations):
                positions, velocities = traffic_simulation(road_length, number_cars, timesteps, v_max=vmax)
                final_velocities = velocities[timesteps,:]    
                flow_rate = sum(final_velocities) / road_length                    
                flow_rates.append(flow_rate)
            
            average_flow_rates.append(np.mean(flow_rates))
            density.append(number_cars / road_length)
            flow_rates = []
        plt.plot(density, average_flow_rates)
        flow_rates = []
        density = []
    for i in range(len(different_vmax)):
        different_vmax[i] = 'v_max = ' + str(different_vmax[i])
    plt.legend(different_vmax)


# 2.2 a) =====================================================================
# plot_simulation(road_length=100, number_cars=10, timesteps=100, title='Plot of time versus position for 10 cars and roadlength 100 meters')
# plot_simulation(road_length=50, number_cars=10, timesteps=100, title='Plot of time versus position for 10 cars and roadlength 50 meters')

# plot_flow_rate(road_length=50, timesteps=100, number_simulations=1000, different_number_cars=np.arange(1,50, 1), title='The fundamental diagram for roadlength 50 meters')
# plot_flow_rate(road_length=100, timesteps=100, number_simulations=1000, different_number_cars=np.arange(2,90, 2), title='The fundamental diagram for roadlength 100 meters')


# 2.1 b) =====================================================================

# print(standard_error_flow_rate(2500))   # this results in a standard error of approximatley 0.001
# plot_standard_error_flow_rate()
        
# 2.2 c) =====================================================================

# plot_flow_rates_different_lengths(timesteps=200, number_simulations=1000)

# 2.2 d) =====================================================================

# plot_simulation(road_length=50, number_cars=10, timesteps=100, title='Plot of time versus position for 10 cars and roadlength\n 50 meters when v_max = 5')
# plot_flow_rate(road_length=50, timesteps=100, number_simulations=1000, different_number_cars=np.arange(2,20, 1), title='The fundamental diagram for roadlength 50 meters when v_max = 5')
# plot_flow_rates_different_vmax(timesteps=100, number_simulations=1000)

# 2.2 e) =====================================================================
    
# plot_simulation(road_length=50, number_cars=10, timesteps=100, title='Plot of time versus position for 10 cars and roadlength\n 50 meters when p = 0.8')
# plot_flow_rate(road_length=50, timesteps=100, number_simulations=1000, different_number_cars=np.arange(2,20, 1), title='The fundamental diagram for roadlength 50 meters when p = 0.8')

