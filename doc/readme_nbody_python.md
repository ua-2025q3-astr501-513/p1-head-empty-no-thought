## Nbody_python.py

The source code for the python integrators. 

It includes a class `nbdsys` for an n-body system to be integrated, and 3 more routines: `rk4`, `leapfrog` and `energy`.

- `nbdsys(nparticles, mass_list = masses, init_pos = init_pos, init_vel = init_vel)`: the initialization method for `nbdsys` class.
  - `nparticles`: numbers of particles in this simulation.
  - `mass_list`: list of masses for each particle. If not specified, this will be set to 1 $M_\odot$ for all the particles.
  - `init_pos`: list of initial position for each particle. If not specified, this will randomly draw n
    samples from a roughly n $\textrm{AU}^3$ space. This will be very tight and is not a good initial condition, so it motivates the user to get better initial conditions.
  - `init_vel`: list of initial velocity for each particle. If not specified, this will be set to 0 for all of them. 
- `energy(system, both = False)`: calculates the energy of the system in its current state.
  - `system`: the system for which the current energy is calculated.
  - `both`: if set to `True`, this method will return an array `[KE, PE]` instead of just a number for total energy (= KE + PE).
- `leapfrog(system, dt)`: using the specified dt, integrate the system for one timestep using the leapfrog method.
  - `system`: the system to be integrated for one timestep.
  - `dt`: the timestep. 
- `rk4(system, dt)`: using the specified dt, integrate the system for one timestep using the 4th order Runge-Kutta method.
  - `system`: the system to be integrated for one timestep.
  - `dt`: the timestep.
- `integrate(system, dt, integrator = 'leapfrog')`: integrates the system using the specified integrator.
  - `system`: the system to be integrated for one timestep.
  - `dt`: the timestep.
  - `integrator`: which integrator to use. This is either `'leapfrog'` or `'rk4'`. 
