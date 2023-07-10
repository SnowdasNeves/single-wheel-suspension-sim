import numpy as np
import matplotlib.pyplot as plt

### DEFINING CONSTANTS ###

### ALL THE CONSTANTS ARE DEFINED IN SI UNITS ###

WHEEL_R = 0.329
WHEEL_MASS = 23
CAR_MASS = 350
SPRING_K = 60000
SPRING_DAMPING = 3929.73  # For critical damping
SPRING_MAX = 1
SPRING_MIN = 0.2
SPRING_EQ = 0.442
BUMP_R = 0.7
BUMP_CENTER_1 = 3
BUMP_CENTER_2 = 6
BUMP_CENTER_3 = 9
BUMP_DISTANCE_TO_SURFACE = 0.5
SIM_DISTANCE = 65

GRAVITY_ACC = -9.81

### ROAD SIMULATION PARAMETERS ###

v_kmh = 20
dx = 0.001
n_bumps = 1

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
    eq_discriminant = BUMP_R**2 - (x - bump_center) ** 2
    if eq_discriminant >= 0:
        if np.sqrt(eq_discriminant) - 0.5 > 0:
            return np.sqrt(eq_discriminant) - 0.5
        else:
            return 0
    else:
        return 0


### EULER-CROMER SIMULATION PARAMETERS ###

v_ms = velocity_in_ms(v_kmh)
dt, sim_t = sim_time(v_ms, dx)
time = np.arange(0, sim_t + dt, dt)

y_car0 = 0.0
v_car0 = 0.0

### SIMULATION FUNCTIONS ###


def simulate(road_surface, spring_damping):
    """
    Uses the Euler-Cromer method to simulate the displacement of the car given
    the wheel displacement and system characteristics.

    :param road_surface: array containing the coordenates of the road surface
    :return: vertical position of the center of the car
    """
    y_wheel = road_surface[:, 1]
    v_wheel = np.zeros_like(y_wheel)
    a_wheel = np.zeros_like(y_wheel)
    f_wheel_spring = np.zeros_like(y_wheel)
    f_ground = np.zeros_like(y_wheel)
    y_car = np.zeros_like(y_wheel)
    v_car = np.zeros_like(y_wheel)

    gravity = (CAR_MASS + WHEEL_MASS) * GRAVITY_ACC

    v_wheel[0] = (y_wheel[1] - y_wheel[0]) / dt
    y_car[0] = y_car0
    v_car[0] = v_car0
    f_ground[0] = gravity

    # Iterates through dt and calculates the displacement of the car

    for i in range(1, int(sim_t / dt) + 1):
        v_wheel[i] = (y_wheel[i] - y_wheel[i - 1]) / dt

        a_wheel[i] = (v_wheel[i] - v_wheel[i - 1]) / dt

        f_wheel_spring[i] = -SPRING_K * y_wheel[i]
        f_ground[i] = gravity + f_wheel_spring[i]

        a_car = -SPRING_K / CAR_MASS * (
            y_car[i - 1] - y_wheel[i - 1]
        ) - spring_damping / CAR_MASS * (v_car[i - 1] - v_wheel[i - 1])

        v_car[i] = v_car[i - 1] + a_car * dt

        y_car[i] = y_car[i - 1] + v_car[i] * dt

        # Restricts string to max and min lengths

        if y_car[i] > y_wheel[i] + (SPRING_MAX - SPRING_EQ):
            y_car[i] = y_wheel[i] + (SPRING_MAX - SPRING_EQ)
        if y_car[i] < y_wheel[i] - (SPRING_EQ - SPRING_MIN):
            y_car[i] = y_wheel[i] - (SPRING_EQ - SPRING_MIN)
    return y_car, f_ground


def spring_compression(y_car, wheel_position):
    """
    Calculates the length of the spring at each time step.

    :param y_car: vertical position of the center of the car
    :param wheel_position: numpy array with the position of the wheel center
    :return: spring length
    """
    wheel_center = wheel_position[:, 1]
    spring_compression = y_car - wheel_center
    return spring_compression


def stabilization(spring_compression):
    """
    Calculates the time it takes to stabilize an oscillation after passing one
    bump.

    :param spring_compression: spring length
    :return: time it takes to stabilize the oscillation
    """
    amp = np.abs(spring_compression - SPRING_EQ)
    max_amp = np.max(np.abs(amp))
    stability_condition = 0.05 * max_amp

    stabilization_time = 0
    for i in range(0, len(amp)):
        if amp[i] > stability_condition:
            stabilization_time = (i + 1) * dt
    return stabilization_time


### PLOTTING DATA ###


def plot_displacements(road_surface, wheel_position, y_car):
    """
    Plots wheel position, road surface and car displacement relative to travel time.

    :param wheel_position: numpy array with the position of the wheel center
    :param road_surface: numpy array with the position of the road surface
    :param y_car: vertical position of the center of the car
    :return: graph of road, wheel and car positions over time
    """
    plt.plot(time, road_surface[:, 1], "r--", label="Road surface")
    plt.plot(time, wheel_position[:, 1], "b", label="Wheel displacement")
    plt.plot(time, y_car, "g", label="Car displacement")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical position (m)")
    plt.title("Road surface and displacements")
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_spring_compression(spring_compression):
    """
    Plots spring compression. The lower the value in yy, the more compressed the
    spring is.

    :param spring_compression: spring length
    :return: graph of spring length over time
    """
    plt.plot(time, spring_compression, "m", label="Spring compression")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical position (m)")
    plt.title("Spring compression, eq. around 0.442 m")
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_force(f_ground):
    """
    Plots the force applied by the wheel on the ground.

    :param force: force applied to the ground on each iteration
    :return: graph of force applied over time
    """
    plt.plot(time, f_ground, "c", label="Force applied")
    plt.xlabel("Time (s)")
    plt.ylabel("Force (N)")
    plt.title("Force applied on the ground")
    plt.legend()
    plt.grid(True)
    plt.show()


### PROGRAM INITIALIZATION ###


if __name__ == "__main__":
    # Plots data using the best damping coefficient and the selected number of bumps

    print("Best damping coefficient:", SPRING_DAMPING, "Ns/m")
    road_surface = create_road(dx, n_bumps)
    wheel_position = np.copy(road_surface)
    wheel_position[:, 1] = road_surface[:, 1] + WHEEL_R
    y_car, f_ground = simulate(road_surface, SPRING_DAMPING)
    y_car += WHEEL_R + SPRING_EQ
    spring_compression = spring_compression(y_car, wheel_position)
    plot_displacements(road_surface, wheel_position, y_car)
    plot_spring_compression(spring_compression)
    plot_force(f_ground)
