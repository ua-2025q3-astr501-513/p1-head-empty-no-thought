import numpy as np
import astropy.units as u
from astropy import constants as const

def randomUniformSpherical(n,
                           nDim,
                           R):
    """
    This method samples n position in a uniform 2D or 3D solid sphere using the information from the Jacobian of
    either the 3D spherical transformation Jacobian or the 2D cylindrical one

    Parameters
    ----------
    n (int): number of arrays.
    nDim (int): dimensions of the space in which lives the sphere (only 2D and 3D available)
    R (float or array): -if float then radii will be sample from 0 to R. 
                        -if array these will be the values of radii used (warning must have length n)

    Returns
    -------
    data (array): 3xn (or 2xn) array with sampled positions.
    """
# For 3D
    if nDim == 3:
        # if needed, we sample the radiuses
        if isinstance(R, np.ndarray) == False:
            R     = R*(np.random.uniform(0,1,n))**(1/3)
        # we sample the angles
        theta = np.arccos(1-2*np.random.uniform(0,1,n)) # polar angle
        phi   = np.random.uniform(0,2*np.pi,n)          # azimuthal angle
        data  = np.array([R*np.cos(phi)*np.sin(theta),  # X
                          R*np.sin(phi)*np.sin(theta),  # Y
                          R*np.cos(theta)])             # Z

# For 2D
    if nDim == 2:
        # if needed, we sample the radiuses
        if isinstance(R, np.ndarray) == False:
            R     = R*(np.random.uniform(0,1,n))**(1/2)
        # we sample the angles
        phi   = np.random.uniform(0,2*np.pi,n)
        data  = np.array([R*np.cos(phi),  # X
                          R*np.sin(phi)]) # Y
    return data


class initialConditions:
    """Description

    Attributes:
        attribute (type): description

    """
    def __init__(self, N):
        """Initialize
        
        Args:
            N (int): number of particles to sample
        
        """
        self.N = N
        self.masses = None

    def sample_piecewise_powerlaw(self, alphas, mass_intervals):
        """
        This method samples N masses from a continuous piecewise power-law PDF.
    
        Parameters
        ----------
        N (int): number of masses to sample.
        alphas (list or array): power-law exponents for each interval.
        mass_intervals (list or array): mass breakpoints, length = len(alphas)+1.
    
        Returns
        -------
        masses (array): array of sampled masses.
        """
        # Step 1. Compute continuity coefficients for the breaks of the power-law
        C               = [1.0]
        for i in range(1, len(alphas)):
            prev_alpha  = alphas[i-1]
            curr_alpha  = alphas[i]
            m_break     = mass_intervals[i]
            C          += [C[-1]*m_break**(prev_alpha-curr_alpha)]
        C               = np.array(C)

        # Step 2. Normalize the coefficients so the whole distribution is normalized
        integrals       = []
        for i, alpha in enumerate(alphas):
            m1, m2      = mass_intervals[i], mass_intervals[i+1]
            if alpha == -1:
                val     = np.log(m2/m1)
            else:
                val     = (m2**(alpha+1)-m1**(alpha+1))/(alpha+1)
            integrals  += [C[i]*val]
        total           = sum(integrals)
        C              /= total

        # Step 3. Set the weights for the sample of intervals
        weights     = np.array(integrals)/total
        cdf_edges   = np.cumsum(weights)

        # Step 4. Sample
        u_all           = np.random.rand(self.N)
        samples         = []

        for i, alpha in enumerate(alphas):
            m1, m2      = mass_intervals[i], mass_intervals[i+1]
            # how many fall in this interval
            mask        = (u_all <= cdf_edges[i]) if i == 0 else (u_all > cdf_edges[i-1]) & (u_all <= cdf_edges[i])
            count       = np.sum(mask)
            if count == 0:
                continue
            u_local     = np.random.rand(count)

            if alpha == -1:
                masses  = m1*(m2/m1)**u_local
            else:
                A       = m1**(alpha+1)
                B       = m2**(alpha+1)
                masses  = (u_local*(B-A)+A)**(1/(alpha+1))
            samples     += [masses]

        self.masses = np.concatenate(samples)
    
    def build_phasespace(self,
                         nDim,
                         R0,
                         Vmax    = None,
                         model   = 'uniform'):
        """
        This method samples a distribution of N phase-space positions following a certain
        distribution.
    
        Parameters
        ----------
        nDim (int): number of dimensions of position space
        R0 (astropy quantity): characteristic physical length
        Vmax (astropy quantity): maximum velocity for the uniform model
        model (str): model for the phase-space DF.
                     available models:
                         - uniform: creates a uniform spherical unstable distribution
                         - plummer: this DF correspond to the Plummer sphere model that obeys
                                    both Poisson and Boltzmann's equations
    
        Returns
        -------
        posVels (array): array with the sampled phase-space positions
        """
        
        # This builds a (unstable) uniform spherical distribution
        if model == 'uniform':
            positions       = randomUniformSpherical(n    = self.N,
                                                     nDim = nDim,
                                                     R    = R0)
            velocities      = randomUniformSpherical(n    = self.N,
                                                     nDim = nDim,
                                                     R    = Vmax)
        
        # This builts a phase-space distribution that follows a Plummer sphere and
        # obeys both Poisson and Boltzmann's equations
        elif model == 'plummer':
            
            M   = (np.sum(self.masses)*u.Msun).to(u.kg)
            G   = const.G
            R0  = R0.to(u.m)
            
            # Radii from the Plummer mass CDF: mr = M(r)/M = r^3/(r^2+a^2)^(3/2)
            mr  = np.random.uniform(size = self.N)
            r   = R0/np.sqrt(mr**(-2/3)-1.0)      
            
            # We built the positions with random directions
            positions = randomUniformSpherical(n    = self.N,
                                               nDim = nDim,
                                               R    = r)
            
            # Local potential and escape speed at each r
            # Plummer potential: Phi(r) = -GM/sqrt(r^2+R0^2),  psi = -Phi
            psi     = G*M/np.sqrt(r**2+R0**2) # relative potential
            vesc    = np.sqrt(2.0*psi)
        
            # Speeds from DF at fixed r:
            # P(v|r) = Av^2f(E) with f(E) = B(psi-v^2/2)^{7/2} for 0 <= v <= vesc
            # With q = v/vesc in [0,1]; then P(q) = Cq^2(1-q^2)^{7/2}.
            # If y = q^2, then y ~ Beta(alpha = 3/2, beta = 9/2) where Beta is the Beta
            # distribution. Hence:
            y = np.random.beta(1.5, 4.5, size=self.N)
            q = np.sqrt(y)
            v = q*vesc
        
            # We built the velocities with random directions
            velocities = randomUniformSpherical(n    = self.N,
                                                nDim = nDim,
                                                R    = v)
        self.posVels    = np.concatenate([positions,velocities],axis=0)
        
        
        
#%%

from astropy import units as u

alphas = [2.03,0.3,-2.35]
mass_intervals = [0.01,0.1,1,50]
R0 = 10*u.au

ic = initialConditions(1000)
ic.sample_piecewise_powerlaw(alphas = alphas, mass_intervals = mass_intervals)
ic.build_phasespace(nDim = 3,R0 = R0, model = 'plummer')

posVels = ic.posVels
masses = ic.masses
particles = np.concatenate((posVels,masses.reshape(1, -1)))

plt.scatter(particles[0],particles[1],c=np.log10(particles[6]),s=1)
plt.show()


plt.scatter(particles[3]/1000,particles[4]/1000,c=np.log10(particles[6]),s=10)
plt.show()



#%%

particles
np.savetxt("initialPositions.txt", particles, fmt="%.8e")























        
        
        
        
        
        