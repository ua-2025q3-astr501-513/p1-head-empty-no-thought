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

def integrate(system, dt = None, tf = None, energy_conservation = True): 
    """
    Integrate the system. 
    dt  - timestep; defult: None. If dt = None, the system will switch to adaptive timestep (RKF4(5)) integration method. 
    tf  - final time; default: None. Used only if dt = None. 
    energy_conservation: if set to True, we will use leapfrog algorithm. If False, RK4 will be used. Used only if dt != None. 
    One of dt or tf must be specified, and dt takes priority over tf (if a timestep is specified, that timestep and
    non-adaptive time step methods will always be used).
    *** NOTE: adaptive (dt = None) is not implemented yet; please use a specified dt. 
    """
    if dt is None:
        # Check if the user gave a final timestep:
        if tf is None:
            print("If no timestep is specified, you must specify an end time tf.")
            return None
        # Implement RKF4(5) here.
        cal_cst = [25./216., 1408./2565., 2197./4104., 1./5.]
        err_cst = [16./135., 6656./12825., 28561./56430., 9./50., 2./55.]

        coeff = np.array((
        [1.0 / 4.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0],
        [1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0],
        [439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0],
        [-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0],
        ))

        # then here comes the timesteps...
        print("I'll implement this if we have time, but RK4 should be good enough along with leapfrog.")
        
    else: 
        # dt is not none: decide if we want energy conservation. 
        if energy_conservation:
            # Initial kick: calculate velocity after half the timestep. 
            a = acceleration(system)
            system.vel += 0.5 * a * dt
            
            # Position drift: calculate position at the end with the half timestep velocity as an approximator. 
            system.pos += dt * system.vel
            
            # Final velocity kick for another half timestep:
            a = acceleration(system) # calculate the acceleration from the new system!
            system.vel += 0.5 * a * dt
            
        else: 
            # RK4 algorithm here. 
            coeff = [0.5, 0.5, 1]
            cst = [1., 3., 3., 1.] / 6
            order = 4

            # Allocate memories:
            x0 = system.pos.copy()
            v0 = system.vel.copy()
            xk = np.empty(order, system.nparticles, 3)
            vk = np.empty(order, system.nparticles, 3)

            # Initial stage
            a = acceleration(system)
            xk[0] = v0
            vk[0] = a

            # Loop to calculate following xk, vk
            # Evaluation at each step is required, so we can't do vectorization here
            for stage in range(1, num_stages):
                # Compute acceleration:
                system.pos = x0 + dt * coeff[stage - 1] * xk[stage - 1]
                a = acceleration(system)
        
                # Then compute xk and vk:
                xk[stage] = v0 + dt * coeff[stage - 1] * vk[stage - 1]
                vk[stage] = a

            # Finally, step forward:
            dx = np.einsum("i,ijk->jk", weights, xk)
            dv = np.einsum("i,ijk->jk", weights, vk)

            # Update system:
            system.pos = precise_add(x0, dt * dx)
            system.vel = precise_add(v0, dt * dv)
            # Update in place so we don't need to return anything. 


def precise_add(x, dx):
    """
    Helper function to reduce roundoff error in the final step for integrators where we add together small numbers. 
    """
    # to reduce roundoff error, use the technique mentioned in 
    # https://alvinng4.github.io/grav_sim/docs/documentations/reducing_round_off_error/
    e = 0.0
    x_0 = x
    e = e + dx
    x = x_0 + e
    e = e + (x_0 - x)
    return (x_0 + dx + e)