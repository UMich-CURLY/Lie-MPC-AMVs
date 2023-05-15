# MPCTools: Nonlinear Model Predictive Control Tools for CasADi (Octave Interface) #

Copyright (C) 2016-2017

Michael J. Risbeck and James B. Rawlings.

MPCTools is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 3, or (at your option) any later
version.

MPCTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file
`COPYING.txt` for more details.

## Installation ##

MPCTools is written for [Octave](https://octave.org) but also supports
[Matlab](https://www.mathworks.com/products/matlab.html). You will also need
to download [CasADi](https://casadi.org). Version requirements are as follows:

* Octave >= 4.0
* Matlab >= R2015b
* CasADi >= 3.1 ([installation instructions here](https://github.com/casadi/casadi/wiki/InstallationInstructions))

Older versions of Octave are not compatible (as they lack `classdef` support).
Older versions of Matlab may be compatible but have not been tested.

For MPCTools, users should download `mpctools.zip` from the
[Downloads](https://bitbucket.org/rawlings-group/octave-mpctools/downloads/)
page and unzip to a convenient location. This zip includes the documentation
pdfs, example scripts, and MPCTools itself. If you wish to clone the repository
and build the documentation yourself, you will need recent versions of GNU Make,
Pandoc, and Python 3, as well as a LaTeX installation.

To use MPCTools, the file `import_mpctools.m` and the directory `+mpctools`
should be placed in the *same* directory, and that directory should be added to
the Octave/Matlab path.

For example, if you unzip `mpctools.zip` to a folder called
`/home/user/octave`, then `/home/user/octave/mpctools` should now contain the
file `import_mpctools.m` and the directory `+mpctools` (among other files).
In Octave (or Matlab), you should run

    addpath('/home/user/octave/mpctools');

To add that folder to the path. You should then be able to run

    mpc = import_mpctools();

See `doc/install.pdf` for more information.

To check that CasADi and MPCTools have both been installed correctly, change
to the appropriate examples directory (either `mpctools/examples-octave` or
`mpctools/examples-matlab`) and run `runall`, which will run all of the
example scripts distributed with MPCTools; plots will appear after the script
finishes (about 2 minutes on standard hardware).

## Documentation ##

Documentation for MPCTools is included in each function and also in the file
`doc/documentation.pdf`. See sample files in the `examples-octave` and
`examples-matlab` folders for complete example scripts.

## Citing MPCTools ##

Because MPCTools is primarily an interface to CasADi, you should cite CasADi as
described on its [website](https://github.com/casadi/casadi/wiki/Publications).
In addition, you can cite MPCTools as

- Risbeck, M.J., Rawlings, J.B., 2016. MPCTools: Nonlinear model predictive
  control tools for CasADi (Octave interface).
  `https://bitbucket.org/rawlings-group/octave-mpctools`.

## Bugs ##

Questions, comments, bug reports, and contributions should be sent to
pratyushkumar@ucsb.edu.

Pratyush Kumar 
<pratyushkumar@ucsb.edu>  
University of California - Santa Barbara  
Department of Chemical Engineering
