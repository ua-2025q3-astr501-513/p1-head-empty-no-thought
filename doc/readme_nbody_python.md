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
  - `system`: the `nbdsys` for which the current energy is calculated.
  - `both`: if set to `True`, this method will return an array `[KE, PE]` instead of just a number for total energy (= KE + PE).
- `leapfrog(system, dt)`: using the specified dt, integrate the system for one timestep using the leapfrog method.
  - `system`: the `nbdsys` to be integrated for one timestep.
  - `dt`: the timestep. 
- `rk4(system, dt)`: using the specified dt, integrate the system for one timestep using the 4th order Runge-Kutta method.
  - `system`: the `nbdsys` to be integrated for one timestep.
  - `dt`: the timestep.
- `integrate(system, dt, integrator = 'leapfrog')`: integrates the system using the specified integrator.
  - `system`: the `nbdsys` to be integrated for one timestep.
  - `dt`: the timestep.
  - `integrator`: which integrator to use. This is either `'leapfrog'` or `'rk4'`.

## make_gif.py

Please note that this file has to be present with `bb_colors.txt` to be used properly, or you will get a file not exist error!

This is the collection of gif-making routine. It is also very slow so be warned. 

- `make_gif(sol_x, filename = "nbody", frames = 100, xlim = None, ylim = None)`: makes a gif with whatever the current default plotting scheme is.
  - `sol_x`: array of positions.
  - `filename` (string): file name of the output file. Code will attach `.gif` for you so no need to keep that.
  - `frames`: number of frames to save.
  - `xlim`, `ylim`: for plotting uses, symmetric with respect to (0, 0, 0). If `None`, the code will use the maximum value in `sol_x` to determine a limit that may or may not be good. 
- `make_gif_bb(sol_x, mass = None, temp = None, filename = "nbody", frames = 100, xlim = None, ylim = None)`: makes a gif but colors the point mass using a very, very crude method of nearest neighbors and may-not-be-valid scaling laws. This is mostly for pretty and for fun.
  - `sol_x`: array of positions.
  - `mass`: list of mass to calculate an approximate blackbody temperature from. If provided with `temp` and `mass`, `temp` will take precedence.
  - `temp`: list of blackbody temperatures for the stars in the simulation. If provided with `temp` and `mass`, `temp` will take precedence.
  - `filename` (string): file name of the output file. Code will attach `.gif` for you so no need to keep that.
  - `frames`: number of frames to save.
  - `xlim`, `ylim`: for plotting uses, symmetric with respect to (0, 0, 0). If `None`, the code will use the maximum value in `sol_x` to determine a limit that may or may not be good. 
