import numpy as np
import matplotlib.pyplot as plt

### DEFINING CONSTANTS ###

### ALL THE CONSTANTS ARE DEFINED IN SI UNITS ###

WHEEL_R = 0.329
WHEEL_MASS = 23
CAR_MASS = 350
SPRING_K = 60000
SPRING_MAX = 1
SPRING_MIN = 0.2
BUMP_R = 0.7
BUMP_CENTER_1 = 3
BUMP_CENTER_2 = 6
BUMP_CENTER_3 = 9
BUMP_DISTANCE_TO_SURFACE = 0.5
SIM_DISTANCE = 65

### SIMULATION PARAMETERS ###

v_kmh = 1
dx = 0.001

### FUNCTIONS ###


def velocity_in_ms(v_kmh):
    """
    Converts velocity from km/h to m/s.

    :param v_kmh: car velocity in km/h
    :return: car velocity in m/s
    """
    v_ms = v_kmh / 3.6
    return v_ms


def sim_time(v_ms, dx):
    """
    Calculates simulation time and time between interations.

    :param v_ms: car velocity in m/s
    :param dx: lenght traveled with each iteration
    :return: simulation time and time between interations
    """
    dt = dx / v_ms
    sim_time = v_ms / SIM_DISTANCE
    return dt, sim_time
