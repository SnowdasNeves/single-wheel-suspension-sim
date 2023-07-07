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


def used_bump(x):
    """Defines which bump is used for wheel to road contact calculation.

    :param x: gives the position of the current iteration
    :return: center of the bump to be used in the current iteration
    """
    if x < 4.5:
        bump_center = BUMP_CENTER_1
    elif x >= 4.5 and x < 7.5:
        bump_center = BUMP_CENTER_2
    else:
        bump_center = BUMP_CENTER_3
    return bump_center


def bump_contact_height(x, bump_center, hypotenuse):
    """
    Calculates the height of the contact point between wheel and road surface.

    :param x: gives the position of the current iteration
    :param bump_center: center of the bump being used for the current iteration
    :param hypotenuse: distance between centers of the wheel and bump
    :return: height of the contact point
    """
    bump_contact_height = 0
    while hypotenuse < WHEEL_R + BUMP_R:
        bump_contact_height += 0.05
        hypotenuse = np.sqrt(
            (np.abs(bump_center - x)) ** 2
            + (BUMP_DISTANCE_TO_SURFACE + WHEEL_R + bump_contact_height) ** 2
        )
    return bump_contact_height


def check_wheel_position(dx):
    """
    Checks if the wheel has reached a bump.

    :param dx: lenght traveled with each iteration
    :return: numpy array with the position of the wheel center
    """
    wheel_position_array = np.zeros((int(SIM_DISTANCE / dx) + 1, 2))
    list_position = 0
    x = 0
    while x <= SIM_DISTANCE:
        bump_center = used_bump(x)
        hypotenuse = np.sqrt(
            (np.abs(bump_center - x)) ** 2 + (BUMP_DISTANCE_TO_SURFACE + WHEEL_R) ** 2
        )
        if hypotenuse >= WHEEL_R + BUMP_R:
            wheel_position_array[list_position] = [x, WHEEL_R]
        elif hypotenuse < WHEEL_R + BUMP_R:
            bump_contact = bump_contact_height(x, bump_center, hypotenuse)
            wheel_position_array[list_position] = [x, (bump_contact + WHEEL_R)]
        list_position += 1
        x += dx
    return wheel_position_array


def plot_wheel_position(wheel_position_array, road_surface_array):
    """
    Plots wheel position and road surface relative to travel time.

    :param wheel_position_array: numpy array with the position of the wheel center
    :param road_surface_array: numpy array with the position of the road surface
    """
    v_ms = velocity_in_ms(v_kmh)
    time = np.arange(0, len(wheel_position_array) * dx / v_ms, dx / v_ms)
    plt.plot(
        time,
        wheel_position_array[:, 1],
        "b-",
        label="Wheel position",
    )
    plt.plot(time, road_surface_array[:, 1], "r--", label="Road surface")
    plt.xlabel("Time (s)")
    plt.ylabel("Position in yy (m)")
    plt.title("Wheel position")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    wheel_position_array = check_wheel_position(dx)
    road_surface_array = wheel_position_array - WHEEL_R
    plot_wheel_position(wheel_position_array, road_surface_array)
