import numpy as np
import matplotlib.pyplot as plt

# Defining constants

sim_time = 30.0
dt = 1
num_iterations = int(sim_time / dt)

mass_car = 350
mass_wheel = 23
r_wheel = 0.329
v_car = 0.2778

spring_k = 60000

bump_height = 0.7
bump_start = 2.3
bump_end = 3.7

g = 9.81

time_frame = np.linspace(0.0, sim_time, num_iterations)
displacement = np.zeros(num_iterations, dtype=float)

# Defining functions


def bump_effect():
    pass


def sim_bump(num_iterations, spring_k, displacement, v_car, dt):
    for i in range(1, num_iterations):
        F_spring = -spring_k * displacement[i - 1]

        acceleration = F_spring / mass_car
        displacement[i] = displacement[i - 1] + v_car * dt
        return displacement


def plot_displacement(time_frame, displacement):
    plt.plot(time_frame, displacement)
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.title("Displacement as a function of time")
    plt.show()


if __name__ == "__main__":
    displacement = sim_bump(num_iterations, spring_k, displacement, v_car, dt)
    print(displacement)
