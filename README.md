## Documentation customization from the HENT Team

The Head Empty, No Thoughts™ team uses three different 3D N-body integrator to integrate toy-model stellar clusters based on a piecewise power-law IMF and a Plummer model. The N-body integrators are the following:
- direct Leap-Frog integrator for small clusters (~100 particles).
- Barnes-Hut algorithm for large clusters (~1000 particles).
- an implicit RK15 from the Rebound Python package for higher accuracy for small clusters (~100 particles).

Additionally, the presentation slides are available [here](https://docs.google.com/presentation/d/1UEkMNe21g3MoQCzrra4jgIt7j6J2d666Kll3J-W7jzA/edit?usp=sharing). 

## Documentation

The code themselves are contained in the `src/` folder. There are 3 main "modules" or "components":

1. `initialConditionsBuilder.py`: used to generate realistic cluster initial conditions for the other integrators to use. This script contains a class with methods to build arbitrary piecewise power-law IMFs and initial conditions in phase space that are either drawn from a virialized Plummer sphere or a uniform spherical distribution in position and velocity spaces. It also has methods to compute the kinetic and potential energy to check the virial ratio of the created cluster.

2. `nbody_python.py`: loads (with some edition; see `nbody_python_demo.ipynb` [here](https://github.com/ua-2025q3-astr501-513/p1-head-empty-no-thought/blob/main/nbody_python_demo.ipynb) for a more in-detail explanation of what's in there) the calculated initial condition, then integrate the system onwards, using either leapfrog or RK4 method. This is not done in the format of a package; instead users can download the script and run it themselves, out of the consideration for a lighter installation. The dependency in this module is only `numpy`. This part of the code referenced [5 steps to N-body simulation](https://alvinng4.github.io/grav_sim/5_steps_to_n_body_simulation/). Additionally, `make_gif.py` is the python script for generating gifs used for presentation, dependent on `matplotlib`. It should be used in conjunction with `bb_colors.txt` (from [m15y](https://temperature.m15y.com/)) for maximum ability - mainly, this adds the ability to make blackbody temperature colored stars in a gif, but I think it's cool.) 

3. (c code readme here...)
astrophysics problem.
We look forward to your results!
