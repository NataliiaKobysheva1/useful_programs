from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


V0 = 37.0
ANGLE_DEG = 75.0
G = 9.81
MASS = 0.145
RADIUS = 0.037
RHO = 1.225
CD = 0.47
DT = 0.005
OUTPUT_DIR = Path(__file__).resolve().parent
AREA = np.pi * RADIUS**2
DRAG_COEFFICIENT = RHO * CD * AREA / (2 * MASS)
ANGLE = np.deg2rad(ANGLE_DEG)
VX0 = V0 * np.cos(ANGLE)
VY0 = V0 * np.sin(ANGLE)


def derivatives(state):
    x, y, vx, vy = state
    speed = np.hypot(vx, vy)
    return np.array([
        vx,
        vy,
        -DRAG_COEFFICIENT * speed * vx,
        -G - DRAG_COEFFICIENT * speed * vy,
    ])


def rk4_step(state, dt):
    k1 = derivatives(state)
    k2 = derivatives(state + 0.5 * dt * k1)
    k3 = derivatives(state + 0.5 * dt * k2)
    k4 = derivatives(state + dt * k3)
    return state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate_with_drag():
    times = [0.0]
    states = [np.array([0.0, 0.0, VX0, VY0])]
    while states[-1][1] >= 0.0 or len(states) == 1:
        next_state = rk4_step(states[-1], DT)
        next_time = times[-1] + DT
        if next_state[1] < 0.0:
            fraction = states[-1][1] / (states[-1][1] - next_state[1])
            next_state = states[-1] + fraction * (next_state - states[-1])
            next_time = times[-1] + fraction * DT
        times.append(next_time)
        states.append(next_state)
        if next_state[1] == 0.0:
            break
    return np.array(times), np.vstack(states)


def main():
    flight_time = 2 * VY0 / G
    time_vacuum = np.linspace(0.0, flight_time, 1200)
    x_vacuum = VX0 * time_vacuum
    y_vacuum = VY0 * time_vacuum - 0.5 * G * time_vacuum**2
    vx_vacuum = np.full_like(time_vacuum, VX0)
    vy_vacuum = VY0 - G * time_vacuum
    time_drag, state_drag = simulate_with_drag()
    x_drag, y_drag, vx_drag, vy_drag = state_drag.T

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax.plot(x_vacuum, y_vacuum, label="Bez oporu", linewidth=2.2)
    ax.plot(x_drag, y_drag, label="Z oporem", linewidth=2.2)
    ax.scatter([x_vacuum.max() / 2], [y_vacuum.max()], s=35, zorder=3)
    ax.scatter([x_drag[y_drag.argmax()]], [y_drag.max()], s=35, zorder=3)
    ax.set(title="Tor rzutu ukośnego", xlabel="x [m]", ylabel="y [m]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(OUTPUT_DIR / "trajectory.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True,
                             constrained_layout=True)
    for ax, vacuum, drag, label in zip(
        axes,
        (vx_vacuum, vy_vacuum, y_vacuum),
        (vx_drag, vy_drag, y_drag),
        ("vx [m/s]", "vy [m/s]", "y [m]"),
    ):
        ax.plot(time_vacuum, vacuum, label="Bez oporu", linewidth=2)
        ax.plot(time_drag, drag, label="Z oporem", linewidth=2)
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
    axes[0].legend()
    axes[-1].set_xlabel("t [s]")
    fig.savefig(OUTPUT_DIR / "time_series.png", dpi=180)
    plt.close(fig)

    vacuum_kinetic = 0.5 * MASS * (vx_vacuum**2 + vy_vacuum**2)
    vacuum_potential = MASS * G * y_vacuum
    drag_kinetic = 0.5 * MASS * (vx_drag**2 + vy_drag**2)
    drag_potential = MASS * G * y_drag
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax.plot(time_vacuum, vacuum_kinetic, label="Ek, bez oporu")
    ax.plot(time_vacuum, vacuum_potential, label="Ep, bez oporu")
    ax.plot(time_vacuum, vacuum_kinetic + vacuum_potential,
            label="E, bez oporu", linewidth=2.4)
    ax.plot(time_drag, drag_kinetic + drag_potential,
            label="E, z oporem", linewidth=2.4)
    ax.set(title="Energia mechaniczna", xlabel="t [s]", ylabel="E [J]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(OUTPUT_DIR / "energy.png", dpi=180)
    plt.close(fig)

    print(f"Pole przekroju: {AREA:.6f} m^2")
    print(f"Współczynnik oporu c: {DRAG_COEFFICIENT:.6f} 1/m")
    print(f"Bez oporu: czas = {flight_time:.3f} s, "
          f"wysokość = {y_vacuum.max():.3f} m, zasięg = {x_vacuum[-1]:.3f} m")
    print(f"Z oporem: czas = {time_drag[-1]:.3f} s, "
          f"wysokość = {y_drag.max():.3f} m, zasięg = {x_drag[-1]:.3f} m")


if __name__ == "__main__":
    main()