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
SIM_DISTANCE = 15

### SIMULATION PARAMETERS ###

v_kmh = 1
dx = 0.001
n_bumps = 0

### INITIAL CALCULATIONS ###


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
    sim_time = SIM_DISTANCE / v_ms
    return dt, sim_time


### DEFINING ROAD SURFACE ###


def create_road(dx, n_bumps):
    """
    Creates the road surface based on the given bump parameters.

    :param dx: lenght traveled with each iteration
    :param bump_center: center of the bump to be used in the current iteration
    :return: array containing the coordenates of the road surface
    """
    v_ms = velocity_in_ms(v_kmh)
    dt, sim_t = sim_time(v_ms, dx)
    road_surface = np.zeros((int(sim_t / dt) + 1, 2))
    t = 0
    x = 0
    list_position = 0
    if n_bumps == 0:
        while t <= sim_t:
            road_surface[list_position, 0] = t
            x += dx
            t += dt
            list_position += 1
    else:
        while t <= sim_t:
            road_surface[list_position, 0] = t
            bump_center = used_bump(x)
            road_surface[list_position, 1] = solve_position_equation(x, bump_center)
            x += dx
            t += dt
            list_position += 1
    return road_surface


def used_bump(x):
    """Defines which bump is used for wheel to road contact calculation.

    :param x: position of the current iteration
    :return: center of the bump to be used in the current iteration
    """
    if x < 4.5:
        bump_center = BUMP_CENTER_1
    elif x >= 4.5 and x < 7.5 and n_bumps > 1:
        bump_center = BUMP_CENTER_2
    elif n_bumps > 2:
        bump_center = BUMP_CENTER_3
    else:
        bump_center = 0
    return bump_center


def solve_position_equation(x, bump_center):
    """
    Calculates y position of road surface for bump geometry.

    :param x: position of the current iteration
    :param bump_center: center of the bump to be used in the current iteration
    :return: y position
    """
    eq_discriminant = 0.49 - (x - bump_center) ** 2
    if eq_discriminant >= 0:
        if np.sqrt(eq_discriminant) - 0.5 > 0:
            return np.sqrt(eq_discriminant) - 0.5
        else:
            return 0
    else:
        return 0


### PLOTTING DATA ###


def plot_road_surface(road_surface, wheel_position):
    """
    Plots wheel position and road surface relative to travel time.

    :param wheel_position: numpy array with the position of the wheel center
    :param road_surface: numpy array with the position of the road surface
    """
    v_ms = velocity_in_ms(v_kmh)
    dt, sim_t = sim_time(v_ms, dx)
    time = np.arange(0, sim_t + dt, dt)
    # plt.plot(road_surface[:, 0], road_surface[:, 1], "r--", label="Road surface xy")
    # plt.plot(wheel_position[:, 0], wheel_position[:, 1], "b", label="Wheel position xy")
    plt.plot(time, road_surface[:, 1], "r--", label="Road surface")
    plt.plot(time, wheel_position[:, 1], "b", label="Wheel position")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical position (m)")
    plt.title("Road surface and Wheel position")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    road_surface = create_road(dx, n_bumps)
    wheel_position = np.copy(road_surface)
    wheel_position[:, 1] = road_surface[:, 1] + WHEEL_R
    plot_road_surface(road_surface, wheel_position)
