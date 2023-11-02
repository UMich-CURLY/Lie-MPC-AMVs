# Lie MPC for AMVs

### Convex Geometric Trajectory Tracking using Lie Algebraic MPC for Autonomous Marine Vehicles
IEEE Robotics & Automation Letters
paper link: https://ieeexplore.ieee.org/document/10301632

<img src="https://github.com/UMich-CURLY/Lie-MPC-AMVs/blob/main/figures/framework3.jpg" width="600">

We propose a computationally-efficient Convex geometric error-state MPC for marine vehicles. The proposed algorithm is compared with Nonlinear MPC.

## Used packages
#### 1) [OSQP](https://osqp.org/docs/index.html) : QP solver 

Folder : `osqp`

Installation is required to run the code.


#### 2) [CasADi](https://web.casadi.org/get/) : NMPC solver 

Folder : `casadi-linux-matlabR2014b-v3.5.5`


#### 3) [MPC-tools](https://bitbucket.org/rawlings-group/octave-mpctools) for CasADi 

Folder : `octave-mpctools`

Tool for easy implementation of NMPC with CasADi by autonomously transforming the given dynamics model into the proper format.


#### 4) USV Otter model from [MSS](https://github.com/cybergalactic/MSS)

Folder : `GNC`, `HYDRO`

Control input of the otter model is changed from motor speed to thrust force.



## Demonstration
Run `main.m` in MATLAB


### Parameters
* Reference trajectory, ocean currents, prediction horizon, simulation parameters are defined in `main.m`
* MPC cost weights (P, Q, R) are defined in `sim.m`
* For calculating dynamics of USV Otter, NMPC uses `otter.m` and NMPC-simple uses `otter-simple.m`. They regard ocean current as 0 m/s. True simulation is processed using `otter_true.m` with ocean currents.


## Results
### Trajectory tracking results
<p float="left">
<img src="https://github.com/UMich-CURLY/Lie-MPC-AMVs/blob/main/figures/result_1.jpg" width="300">
<img src="https://github.com/UMich-CURLY/Lie-MPC-AMVs/blob/main/figures/result_2.jpg" width="300">
</p>


### Computation time (msec) for single optimization
| Ocean Current         | Proposed MPC | NMPC |  NMPC-simple |
|-----------------|:--------:|:--------:|:-----:|
| 0 m/s      |   49 |   764  | 478 |
| 0.5 m/s    |   50  |   953  | 460 |

