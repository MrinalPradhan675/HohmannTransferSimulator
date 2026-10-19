# Hohmann Transfer Orbit Simulator

A MATLAB simulator that finds the two burns of a Hohmann transfer between user-defined circular orbits by integrating two-body motion with ODE45, and checks the results against the closed-form delta-v, transfer time, and conservation of energy and angular momentum.

## Why This Is Not Just The Formula

The Hohmann equations give delta-v directly. This project does not use them to find the answer. It uses them to check the answer.

1. A shooting method (`fzero`) searches for the departure speed whose transfer orbit reaches the target radius at apoapsis. Each guess is tested by integrating the two-body equations of motion with ODE45.
2. The arrival burn is the difference between the target circular speed and the speed ODE45 actually reaches at apoapsis.
3. The numerical delta-v values and transfer time are compared with the analytic Hohmann solution.
4. Specific energy and specific angular momentum are tracked along the path to check the integration.

Agreement between steps 1 to 3 and the formulas means the numerical model is correct, since the two methods are independent.

## Results (Default Case: 400 km LEO to GEO)

| Quantity | Numerical | Analytic | Error [%] |
|---|---|---|---|
| Departure Dv [m/s] | 2397.4703 | 2397.4703 | 1.26e-10 |
| Arrival Dv [m/s] | 1456.4868 | 1456.4868 | 4.98e-10 |
| Total Dv [m/s] | 3853.9571 | 3853.9571 | 2.67e-10 |
| Transfer Time [s] | 19048.48 | 19048.48 | 6.39e-10 |

- Transfer time: 5.291 hr
- Max relative energy drift: 2.56e-12
- Max relative angular momentum drift: 8.43e-13

## Usage

1. Open `hohmann_transfer_simulator.m` in MATLAB.
2. Run the script.
3. Enter values at each prompt, or press Enter to keep the default.

| Input | Default | Notes |
|---|---|---|
| Gravitational parameter [m^3/s^2] | 3.986004418e14 | Earth. Moon: 4.9048695e12. Mars: 4.282837e13 |
| Body radius [m] | 6378137 | Moon: 1737400. Mars: 3389500 |
| Initial altitude [km] | 400 | Circular orbit |
| Target altitude [km] | 35785.863 | GEO for Earth |

The script prints a results table and opens a figure with the transfer geometry, speed along the transfer, and the conservation check.

## Requirements

- MATLAB R2019b or later (uses `tiledlayout`)
- No additional toolboxes

## Method

**Dynamics.** Two-body motion in a plane, state `[x; y; vx; vy]`:

```
r_ddot = -mu * r_vec / |r|^3
```

**Integration.** `ode45` with `RelTol = 1e-12` and `AbsTol = 1e-9`. An event function stops the integration at apoapsis, where the radial velocity crosses zero from positive to negative.

**Shooting.** The initial state is periapsis at `(r1, 0)` moving in `+y` with speed `vp`. `fzero` solves for `vp` so that the apoapsis radius equals `r2`. The bracket runs from just above circular speed to just below escape speed.

**Analytic reference.**

```
a_t    = (r1 + r2) / 2
dv1    = sqrt(mu * (2/r1 - 1/a_t)) - sqrt(mu / r1)
dv2    = sqrt(mu / r2) - sqrt(mu * (2/r2 - 1/a_t))
t_hohm = pi * sqrt(a_t^3 / mu)
```

## Assumptions And Limitations

- Impulsive burns (instantaneous velocity changes)
- Coplanar, circular initial and target orbits
- Two-body gravity only. No J2, drag, or third-body effects
- Ascending transfers only (target above initial orbit)
- The solver needs the target radius below about 49 times the initial radius. The script checks this and exits with a message if it is exceeded
- Central body presets are for convenience. Large altitudes around small bodies are not physically meaningful because the model ignores other gravitating bodies

## Possible Extensions

- Bi-elliptic transfer and the Hohmann crossover near a radius ratio of 11.94
- Plane change burns
- Finite burn durations and gravity losses
- J2 perturbation

## Files

| File | Description |
|---|---|
| `hohmann_transfer_simulator.m` | Simulator, validation, and plots |
