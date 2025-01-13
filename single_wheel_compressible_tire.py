import numpy as np
import matplotlib.pyplot as plt

WHEEL_RADIUS = 0.329
WHEEL_MASS = 23
CAR_MASS = 350
TIRE_K = 191000  # K is Elastic Coefficient
SPRING_K = 60000
SPRING_DAMPING = 8000  # Value for critical damping
SPRING_MAX_LEN = 1
SPRING_MIN_LEN = 0.2
SPRING_EQ_LEN = 0.442
BUMP_RADIUS = 0.7
BUMP_CENTER_1 = 3  # Position of successive bumps on the road
BUMP_CENTER_2 = 6
BUMP_CENTER_3 = 9
BUMP_DISTANCE_TO_SURFACE = 0.5
SIM_DISTANCE = 65  # Distance traveled by the car
GRAVITY_ACCEL = -9.81

# Defining simulation conditions
velocity_kmh = 20
dx = 0.001  # Distance traveled between iterations
n_bumps = 3

velocity_ms = velocity_kmh / 3.6


def simulation_time():
    dt = dx / velocity_ms
    sim_time = SIM_DISTANCE / velocity_ms
    return dt, sim_time


def create_road(dt, sim_t):
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
    if x < 4.5:
        bump_center = BUMP_CENTER_1
    elif 4.5 <= x < 7.5 and n_bumps > 1:
        bump_center = BUMP_CENTER_2
    elif n_bumps > 2:
        bump_center = BUMP_CENTER_3
    else:
        bump_center = 0
    return bump_center


def solve_position_equation(x, bump_center):
    eq_discriminant = BUMP_RADIUS**2 - (x - bump_center) ** 2
    if eq_discriminant >= 0:
        if np.sqrt(eq_discriminant) - 0.5 > 0:
            return np.sqrt(eq_discriminant) - 0.5
        else:
            return 0
    else:
        return 0


def simulate(road_surface, spring_damping, sim_time):
    y_road = road_surface[:, 1]
    y_wheel = np.zeros_like(y_road)
    v_wheel = np.zeros_like(y_road)
    y_car = np.zeros_like(y_road)
    v_car = np.zeros_like(y_road)
    f_wheel_spring = np.zeros_like(y_road)
    f_ground = np.zeros_like(y_road)

    gravity = (CAR_MASS + WHEEL_MASS) * GRAVITY_ACCEL

    y_wheel[0] = y_wheel0
    v_wheel[0] = v_wheel0
    y_car[0] = y_car0
    v_car[0] = v_car0
    f_ground[0] = gravity

    # Iterates through dt and calculates the displacement of the car
    for i in range(1, int(sim_time / dt) + 1):
        f_wheel_spring[i] = -SPRING_K * y_wheel[i]
        f_ground[i] = gravity + f_wheel_spring[i]

        a_car = -SPRING_K / CAR_MASS * (
            y_car[i - 1] - y_wheel[i - 1]
        ) - spring_damping / CAR_MASS * (v_car[i - 1] - v_wheel[i - 1])

        a_wheel = (
            spring_damping / WHEEL_MASS * (v_car[i - 1] - v_wheel[i - 1])
            + SPRING_K / WHEEL_MASS * (y_car[i - 1] - y_wheel[i - 1])
            - TIRE_K / WHEEL_MASS * (y_wheel[i - 1] - y_road[i - 1])
        )

        v_car[i] = v_car[i - 1] + a_car * dt
        v_wheel[i] = v_wheel[i - 1] + a_wheel * dt

        y_car[i] = y_car[i - 1] + v_car[i] * dt
        y_wheel[i] = y_wheel[i - 1] + v_wheel[i] * dt

        # Restricts string to max and min lengths

        if y_car[i] > y_wheel[i] + (SPRING_MAX_LEN - SPRING_EQ_LEN):
            y_car[i] = y_wheel[i] + (SPRING_MAX_LEN - SPRING_EQ_LEN)
        if y_car[i] < y_wheel[i] - (SPRING_EQ_LEN - SPRING_MIN_LEN):
            y_car[i] = y_wheel[i] - (SPRING_EQ_LEN - SPRING_MIN_LEN)
    return y_wheel, y_car, f_ground


def spring_compression(y_car, y_wheel):
    spring_compression = y_car - y_wheel
    return spring_compression


def stabilization(spring_compression):
    amp = np.abs(spring_compression - SPRING_EQ_LEN)
    max_amp = np.max(np.abs(amp))
    stability_condition = 0.05 * max_amp
    stabilization_time = 0
    for i in range(0, len(amp)):
        if amp[i] > stability_condition:
            stabilization_time = (i + 1) * dt
    return stabilization_time


def plot_displacements(road_surface, y_wheel, y_car):
    plt.plot(time, road_surface[:, 1], "r--", label="Road surface")
    plt.plot(time, y_wheel, "b", label="Wheel displacement")
    plt.plot(time, y_car, "g", label="Car displacement")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical position (m)")
    plt.title("Road surface and displacements")
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_spring_compression(spring_compression):
    plt.plot(time, spring_compression, "m", label="Spring compression")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical position (m)")
    plt.title("Spring compression, eq. around 0.442 m")
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_force(f_ground):
    plt.plot(time, f_ground, "c", label="Force applied")
    plt.xlabel("Time (s)")
    plt.ylabel("Force (N)")
    plt.title("Force applied on the ground")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    dt, sim_time = simulation_time()
    time = np.arange(0, sim_time + dt, dt)

    y_wheel0 = 0.0
    v_wheel0 = 0.0
    y_car0 = 0.0
    v_car0 = 0.0

    road_surface = create_road(dt, sim_time)

    spring_damping = np.arange(100.0, 20000.0, 100.0)
    best_damping = None
    low_stabilization_time = float("inf")
    for i in spring_damping:
        y_wheel, y_car, f_ground = simulate(road_surface, i, sim_time)
        y_wheel += WHEEL_RADIUS
        y_car += WHEEL_RADIUS + SPRING_EQ_LEN
        spring_comp = spring_compression(y_car, y_wheel)

        stabilization_time = stabilization(spring_comp)

        if best_damping is None or stabilization_time < low_stabilization_time:
            best_damping = i
            low_stabilization_time = stabilization_time

    print("Best damping coefficient:", best_damping, "Ns/m")
    road_surface = create_road(dt, sim_time)
    y_wheel, y_car, f_ground = simulate(road_surface, best_damping, sim_time)
    y_wheel += WHEEL_RADIUS
    y_car += WHEEL_RADIUS + SPRING_EQ_LEN
    spring_compression = spring_compression(y_car, y_wheel)
    plot_displacements(road_surface, y_wheel, y_car)
    plot_spring_compression(spring_compression)
    plot_force(f_ground)
