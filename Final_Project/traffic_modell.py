#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 08:04:55 2020

@author: viggomoro
"""

import random
import matplotlib.pyplot as plt
import numpy as np

class Car:
    """Class representing the cars in the simulation. It holds all necessary 
    data for every car."""
    
    def __init__(self, position, velocity, v_max, p, lane):
        self.position = position    # Position of the car
        self.velocity = velocity    # Velocity of the car
        self.v_max = v_max          # Maximum speed for the car
        self.p = p                  # Probability of speed reduction
        self.lane = lane            # Current lane of the car
        
        
def compute_distance_ahead(cars, current_index, desired_lane, road_length):
    """Computes the distance to the car closest ahead in the lane desired_lane. 
    cars is a list of Car objects, current_index is the index of the car which 
    velocity is currently being updated in traffic_simulation and road_length 
    is the road_length in traffic_simulation"""
    
    cars_rearanged = cars[current_index:] + cars [0:current_index]  # Arange cars starting with the current one
    car = cars_rearanged[0]
    
    # Calculate distance to car in front. If there is only one car in the lane the distance is set to road_length
    distance = road_length
    for current_car in cars_rearanged[1:]:
        if current_car.lane == desired_lane:
            distance = (current_car.position - car.position) % road_length
            break
    return distance
    

def compute_distance_behind(cars, current_index, desired_lane, road_length):
    """Identical to compute_distance_ahead with the exception of calculating 
    the distance to the car closest behind."""
    
    cars_rearanged = cars[current_index+1:] + cars[0:current_index+1] # Arange cars ending with the current one
    car = cars_rearanged[-1]
    
    # Calculate distance to car behind. If there is only one car in the lane the distance is set to road_length
    distance = road_length
    for current_car in cars_rearanged[::-1][1:]:
        if current_car.lane == desired_lane:
            distance = (car.position - current_car.position) % road_length
            break
    return distance


def switch_lane_left(cars, current_index, road_length, number_lanes):
    """Returns a boolean indicating if it's possible to change lane to the left i.e 
    overtake or not. cars is a list of Car objects, current_index is the index 
    of the car which velocity is currently being updated in traffic_simulation 
    and road_length is the road_length in traffic_simulation. number_lanes is 
    the number of lanes in the simulation."""
    
    car = cars[current_index]  # Car currently being updated in traffic_simulation
    # If the current lane is number_lanes there is no lane to overtake in
    if car.lane == number_lanes:
        return False
    
    # Compute distances to the closest cars in front and behind in the desired lanes
    dist_ahead = compute_distance_ahead(cars, current_index, car.lane+1, road_length)
    dist_behind = compute_distance_behind(cars, current_index, car.lane+1, road_length)
    
    # For an overtake to be possible it shouldn't risk causing collitions
    if car.velocity < dist_ahead and cars[current_index-1].velocity < dist_behind:
        return True
    else:
        return False
    
    
def switch_lane_right(cars, current_index, road_length):
    """Identical to switch_lane_left with the exception that it now checks if 
    it's possible to change lane to the right. For example after completing an 
    overtake."""
    car = cars[current_index]  # Car currently being updated in traffic_simulation
    # If the current lane is the most to the right you can't switch further to the right
    if car.lane == 1:
        return False
    
    # Compute distances to the closest cars in fron and behind
    dist_ahead = compute_distance_ahead(cars, current_index, car.lane-1, road_length)
    dist_behind = compute_distance_behind(cars, current_index, car.lane-1, road_length)
    
    # For an overtake to be possible it shouldn't risk causing collitions
    if car.velocity < dist_ahead and cars[current_index-1].velocity < dist_behind:
        return True
    else:
        return False
    
    
def simulate_one_timestep(cars, number_cars, lanes, road_length):
    """Simulates one time step of the graffic simulation. cars is a list of Car
    objects, number_cars is the number of cars in the simulation, lanes is the 
    number of lanes in the simulation and road_length is the length of the 
    periodic road in the simulation."""
    
    # Update velocies
    for i in range(number_cars):
        car = cars[i]     # the current car
        
        # accelerate on step towards the maximum speed of the current car
        if car.velocity < car.v_max:
            car.velocity += 1
            
        distance = compute_distance_ahead(cars, i, car.lane, road_length) # Distance to the car ahead
        # breake or overtake if to close to car in front
        if car.velocity >= distance:    
            overtrake = switch_lane_left(cars, i, road_length, lanes)  # boolean to indicate if it's possible to overtake
            if overtrake:
                car.lane += 1
            else:
                car.velocity = distance - 1
                
        # random braking, reduce the velocity of a moving car by one unit
        if random.random() <= car.p and car.velocity >= 1 :
            car.velocity -= 1
        
    # move cars forwards, the road is a periodic i.e  a circle
    for i in range(number_cars):
        cars[i].position += cars[i].velocity 
        cars[i].position %= road_length     # Create periodicity
        
    # Sort cars based on positions
    cars = sorted(cars, key=lambda car: car.position)
    
    # The cars change one lane to the right if it won't risk causing a collisoin
    for i in range(number_cars):
        change_lane = switch_lane_right(cars, i, road_length)
        if change_lane:
            cars[i].lane -= 1
    
    return cars
    
        
def traffic_simulation(road_length=100, number_cars=10, timesteps=500, lanes=1):
    """Simulates traffic according to the cellular automation traffic modell 
    with overtaking and different maximum velocities and speed reduction 
    probabilities. number_cars is the number of cars in the simulation, lanes 
    is the number of different lanes in the simulation, road_length is the length of the 
    periodic road in the simulation and timesteps is the number of steps the 
    simulation runs."""

    cars = []  # List to store Car objects
    
    # Create list for diffent v_max
    different_v_max = [2,3]
    num = np.ceil(number_cars / 2)
    different_v_max = int(num) * different_v_max
    different_v_max = different_v_max[:number_cars]

    # Initial conditions
    for i in range(number_cars):
        pos = int(road_length * random.random())
        vel = int(3 * random.random()) + 1
        v_max = 2#different_v_max[i] #int(3 * random.random()) + 1
        p = 0.5
        lane = 1
        new_car = Car(pos, vel, v_max, p, lane)
        cars.append(new_car)
        
    # Sort cars based on positions
    cars = sorted(cars, key=lambda car: car.position)
    
    # Test
    for i in range(number_cars):
        cars[i].v_max = different_v_max[i]
        
    # Run simulation
    for t in range(timesteps):
        # Simulate one time step
        cars = simulate_one_timestep(cars, number_cars, lanes, road_length) # Kan testa att bara kalla på funktionen utan att retunera något
        
    return cars
    
    
def plot_flow_rate(road_length=100, timesteps=100, number_simulations=10, different_number_cars=np.arange(1,50), lanes=1):
    """Plots flow rate versus car density, i.e the fundamental diagram. The flow
    rate is defined as the sum of the cars verlocities at equilibrium divided 
    by the road length. The car density is defined as the number of cars 
    divided by road length."""
    flow_rates = []
    average_flow_rates = []
    density = []
    standard_errors = []
    
    # Calculate flow rate and density
    for number_cars in different_number_cars:
        # Do it multiple times for and take average
        for simulation in range(number_simulations):
            cars = traffic_simulation(road_length, number_cars, timesteps, lanes)
            sum_velocities = sum(car.velocity for car in cars)
            flow_rate = sum_velocities / road_length
            flow_rates.append(flow_rate)
        
        # Standard error
        standard_error = np.sqrt(np.var(flow_rates) / (number_simulations- 1))
        standard_errors.append(standard_error)
        
        average_flow_rates.append(np.mean(flow_rates))
        density.append(number_cars / road_length)
        flow_rates = []
        
    plt.figure()
    plt.title('The fundamental diagram for ' + str(lanes) + ' lanes and a road length of ' + str(road_length) + '. The\n different v_max are with equal proportions from [2]')
    plt.xlabel('Density')
    plt.ylabel('Flow rate')
    plt.errorbar(density, average_flow_rates, yerr=np.array(standard_errors))
    
    return density, average_flow_rates, standard_errors


def standard_error_flow_rate(number_simulations, road_length=50, number_cars=25, timesteps=100, lanes=1, flow_rates=[]):
    """Computes the standard error of the flow rate. number_simulations is the 
    number of samples of the flow rate that is generated. flow_rates is the 
    list where the flow rate is stored for each simulation. It's included as 
    an argument to utilize past results when plotting the standard error."""

    for i in range(number_simulations):
        cars = traffic_simulation(road_length, number_cars, timesteps, lanes)
        sum_velocities = sum(car.velocity for car in cars)
        flow_rate = sum_velocities / road_length
        flow_rates.append(flow_rate)
     
    # Calculate the standard error of the flow rate    
    variance = np.var(flow_rates)  # Total simulations up to this point
    total_number_simulations = len(flow_rates)
    standard_error = np.sqrt(variance / (total_number_simulations - 1))
    return standard_error


def plot_standard_error_flow_rate(road_length, number_cars, timesteps, lanes):
    """Plots the standard error of the flow rate versus the size of the sample, 
    i.e number of simulations."""
    
    # Intervall between number of simulations in list
    intervall = 25
    different_number_simulations = np.arange(intervall, 3101, intervall)
    
    # Stor in lists
    standard_errors = []
    flow_rates = []  # This list is an argument to standard_error_flow_rate to utilize past results
    
    for number_simulations in different_number_simulations:
        standard_error = standard_error_flow_rate(intervall, road_length, number_cars, timesteps, lanes, flow_rates)
        standard_errors.append(standard_error)
    
    plt.figure()
    plt.title('Standard error for ' + str(lanes) + ' lanes, ' + str(number_cars) + ' cars and a road length of ' + str(road_length) + '. The\n different v_max are with equal proportions from [1, 2, 3]')
    plt.xlabel('Number of simulations')
    plt.ylabel('Standard error')
    plt.plot(different_number_simulations, standard_errors)
    return list(different_number_simulations), list(standard_errors)


# =============================================================================
# Comparison of flow rate

# from time import process_time 
# t1_start = process_time() 
# list1 = list(range(1,50,2))
# den1, flow1, yerr1 = plot_flow_rate(50,100,1000, list1, lanes=1)

# list2 = list(range(1,25, 3)) + list(range(26,75, 2)) + list(range(76,101, 3))
# den2, flow2, yerr2 = plot_flow_rate(50,150,500, list2, lanes=2)

# list3 = list(range(1,50, 3)) + list(range(51,103, 2)) + list(range(103,151, 3))
# den3, flow3, yerr3 = plot_flow_rate(50,250,250, list3, lanes=3)    # större mellanrum speciellt i början och i slutet
# t1_stop = process_time() 
# time = t1_stop - t1_start

# =============================================================================
# Standard error estimates

# sim1, stan_err1 = plot_standard_error_flow_rate(50, 25, 100, 1)
# sim2, stan_err2 = plot_standard_error_flow_rate(50, 40, 150, 2)
# sim3, stan_err3 = plot_standard_error_flow_rate(50, 52, 250, 3)


