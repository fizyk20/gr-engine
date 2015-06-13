# gr-engine
A simple physics engine capable of simulating General Relativity.

## Features
- Kerr and Schwarzschild spacetimes with separate coordinate systems for near-pole regions increasing accuracy
- Propagation of point particles and entities with orientation
- Integration of the equation of motion with either Runge-Kutta 4 or Dormand-Prince integrators

## Tests/Examples
Folder "tests" contains some examples:
- A sine wave integrated with RK4/DP
- Shapiro delay calculator - by propagating a photon near the sun and reading the round-trip time

## Documentation
Some documentation of the available classes is provided at http://fizyk20.github.io/gr-engine
