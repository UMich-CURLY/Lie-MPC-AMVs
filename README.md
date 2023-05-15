# Geometric trajectory tracking MPC for Autonomous Marine Vehicles
![framework](figures/framework3.jpg?raw=true "Title" | width=100)

We propose a Convex geometric error-state MPC for marine vehicles. The proposed algorithm is compared with Nonlinear MPC.

## Used packages
#### OSQP : QP solver (https://osqp.org/docs/index.html)
'osqp' folder
installation is required to run the code.

#### CasADi : NMPC solver (https://web.casadi.org/get/)
'casadi-linux-matlabR2014b-v3.5.5' folder

#### MPC-tools for CasADi (https://bitbucket.org/rawlings-group/octave-mpctools)
'octave-mpctools'

#### USV Otter model from MSS (https://github.com/cybergalactic/MSS)
'GNC, HYDRO'

#### Original Error-state MPC code (https://github.com/UMich-CURLY/Error-State-MPC)
'utils, controllers'

## Demonstration
run 'main.m' in MATLAB.

### Parameters
* Reference trajectory, ocean currents, prediction horizon, simulation parameters are defined in 'main.m'
* MPC cost weights (P, Q, R) are defined in 'sim.m'
* For calculating dynamics of USV Otter, NMPC uses 'otter.m' and NMPC-simple uses 'otter-simple.m'. They regard ocean current is 0 m/s. True simulation is processed using 'otter_true.m' with ocean currents.

