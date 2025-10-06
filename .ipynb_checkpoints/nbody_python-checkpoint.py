import numpy as np
import scipy as sp

# global constant definition
msun = 1.989e30 # kg
pc   = 3.086e13 # km

G    = sp.constants.G * (1000)**3 * (1/msun) # in km^3 Msun^-1 s^-2

# Class implementation for the nbody system
class nbdsys:

    def __init__(self, nparticles, mass_list = None, init_vel = None, init_pos = None, test = False) -> None:
        """
        Initializes an object of the class nbdsys. 
        Note that nparticles has to be defined at the time of object generation; for others, the program provides a uniform distribution
        of 1 solar mass objects in n**(1/3) pc**3 with 0 initial velocity and uniformly random initial positions. 

        Note that after the list was put in, the stars will be recentered so that CoM is always at (0, 0, 0).
        """
        
        # check if the lists are of appropriate size. 
        # ***********
        # ** TO DO **
        # ***********

        self.nparticles = nparticles

        # attach the attributes.
        if mass_list is None:
            mass_list = np.ones(nparticles)
            print("Mass list is not provided; initiated as 1 Msun.")
        else: 
            mass_list = np.asarray(mass_list)
        self.mlist = mass_list

        if init_vel is None:
            init_vel = np.zeros((nparticles, 3))
            print("Velocity list is not provided; initiated as 0.")
        else: 
            init_vel = np.asarray(init_vel)
        self.vel = init_vel

        if init_pos is None:
            init_pos = np.random.uniform(low=-1*nparticles**(1/3), high=nparticles**(1/3), size=(nparticles, 3)) * pc
            print("Positions list is not provided; initiated as uniform random distribution.")
        else: 
            init_pos = np.asarray(init_pos)
        self.pos = init_pos

        # Now this class will center the system before beginning the integration.
        self.center()

    def center(self) -> None:
        # calculate the center of mass as-is. 
        
        M = np.sum(self.mlist) # total mass of the system
        
        x_cm = np.einsum("i,ij->j", self.mlist, self.pos) / M # einsum: sum(ith entry in mlist * ijth entry on pos) = ith entry on final answer.
        v_cm = np.einsum("i,ij->j", self.mlist, self.vel) / M

        self.pos -= x_cm
        self.vel -= v_cm

# methods related to this class?
def acceleration(system): 
    # Again use the einsum method detailed in the reference document.
    # Note that this returns km s^-2.

    # Initialize the acceleration array. 
    a = np.zeros((system.nparticles, 3))

    pos = system.pos
    m = system.mlist

    # Step 1: broadcast pos into (N, N, 3) to compute pairwise distance vectors. 
    r_ij = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]

    # Step 2: compute the norm.
    r_norm = np.linalg.norm(r_ij, axis=2)

    # Step 3: compute 1/R_norm**3. For the diagonal 0 elements, set them all to 0 directly. 
    with np.errstate(divide='ignore', invalid='ignore'):
        inv_r_cubed = 1.0 / (r_norm * r_norm * r_norm)
    np.fill_diagonal(inv_r_cubed, 0.0)

    # Step 4: Calculate the gravitational acceleration. 
    # Use einsum here: we would end up with a Nx3 matrix.
    # for r_ij, we need all 3 dimensions (this tells us the direction). -> ijk
    # for r_norm, we use two dimensions and broadcast to the 3rd (this tells us the distance between pairs) -> ij
    # Finally for m, we broadcast it from 1 to 3 dimensional -> i
    # For output, we want to sum along the 1st dimension (each point mass occupies 1 "column") -> jk.
    # Putting this into einsum:
    a[:] = G * np.einsum("ijk,ij,i->jk", r_ij, inv_r_cubed, m)

    # finally because I don't understand the reference document's class coding style, return the acceleration
    return a

# well. I forgot to push but. Here is the bare bones of one of the integration algorithms.
# \Leapfrog method/: it conserves energy in the long term, and we will probably be running pretty long term. 

def leapfrog(system, dt):

    # Initial kick: calculate velocity after half the timestep. 
    a = acceleration(system)
    system.vel += 0.5 * a * dt
    
    # Position drift: calculate position at the end with the half timestep velocity as an approximator. 
    system.pos += dt * system.vel
    
    # Final velocity kick for another half timestep:
    a = acceleration(system) # calculate the acceleration from the new system!
    system.vel += 0.5 * a * dt