import numpy as np
import matplotlib.pyplot as plt

# Defining constants

sim_time = 30.0
dt = 0.01
num_iterations = int(sim_time / dt)

mass_car = 350
mass_wheel = 23
r_wheel = 0.329

spring_k = 60000

r_bump = 0.7

g = 9.81

# Defining functions
