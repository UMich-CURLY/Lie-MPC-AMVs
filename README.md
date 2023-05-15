# Lie MPC for AMVs

### Convex Geometric Trajectory Tracking using Lie Algebraic MPC for Autonomous Marine Vehicles

<img src="https://github.com/UMich-CURLY/Lie-MPC-AMVs/blob/main/figures/framework3.jpg" width="600">
We propose a Convex geometric error-state MPC for marine vehicles. The proposed algorithm is compared with Nonlinear MPC.

## Used packages
#### OSQP : QP solver 
(https://osqp.org/docs/index.html)

`osqp`

Installation is required to run the code.

#### CasADi : NMPC solver 
(https://web.casadi.org/get/)

`casadi-linux-matlabR2014b-v3.5.5`

#### MPC-tools for CasADi 
(https://bitbucket.org/rawlings-group/octave-mpctools)

`octave-mpctools`

Tool for easy implementation of NMPC with CasADi by autonomously transforming the given dynamics model into the proper format.

#### USV Otter model from MSS 
(https://github.com/cybergalactic/MSS)

`GNC, HYDRO`

Control input of the otter model is changed from motor speed to thrust force.

#### Original Error-state MPC code 
(https://github.com/UMich-CURLY/Error-State-MPC)

`utils`

## Demonstration
run `main.m` in MATLAB.

### Parameters
* Reference trajectory, ocean currents, prediction horizon, simulation parameters are defined in `main.m`
* MPC cost weights (P, Q, R) are defined in `sim.m`
* For calculating dynamics of USV Otter, NMPC uses `otter.m` and NMPC-simple uses `otter-simple.m`. They regard ocean current as 0 m/s. True simulation is processed using `otter_true.m` with ocean currents.

