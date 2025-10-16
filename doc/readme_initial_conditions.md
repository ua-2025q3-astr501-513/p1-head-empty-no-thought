# **cluster_dynamics/initial_conditions.py**

Generate and manage particle distributions in phase space for N-body simulations.

---

### **function: randomUniformSpherical**(n, nDim, R)

Sample *n* random positions inside a solid 2D or 3D sphere.

- **Parameters:**
  - `n` – Number of points to sample.
  - `nDim` – Dimensionality of the space (2 or 3).
  - `R` – Sphere radius or array of radii.  
    - If a scalar, radii are sampled randomly from 0 to *R*.  
    - If an array, each element is treated as a fixed radius.

- **Returns:**
  - `data` – 2×n or 3×n NumPy array of sampled positions.

- **Description:**
  - Uses the Jacobian of the coordinate transformation to ensure uniform density.  
  - In 3D, samples both polar and azimuthal angles.  
  - In 2D, samples angles uniformly over the circle.

---

### **function: kinetic_energy**(m, vx, vy, vz)

Compute individual kinetic energies for a set of particles.

- **Parameters:**
  - `m` – Particle masses (astropy quantity).  
  - `vx`, `vy`, `vz` – Velocity components (astropy quantities).

- **Returns:**
  - `K` – Array of kinetic energies in Joules.

- **Description:**
  - Converts all inputs to SI units.  
  - Subtracts mean velocity from each component to remove bulk motion.  
  - Computes kinetic energy for each particle.

---

### **function: potential_energy**(m, x, y, z)

Compute total gravitational potential energy of a set of particles.

- **Parameters:**
  - `m` – Particle masses (astropy quantity).  
  - `x`, `y`, `z` – Position components (astropy quantities).

- **Returns:**
  - `U` – Total potential energy in Joules.

- **Description:**
  - Converts all quantities to SI units.  
  - Computes pairwise distances between particles.  
  - Accumulates total potential energy from all pairs.  
  - Prints warnings if any particles overlap.

---

## **class: initialConditions**

Create, store, and sample initial phase-space configurations for an N-body system.

- **Attributes:**
  - `N` – Number of particles.  
  - `masses` – Array of sampled masses.  
  - `C` – Continuity coefficients used in piecewise normalization.  
  - `posVels` – Combined array of positions and velocities.

---

### **method: sample_piecewise_powerlaw**(alphas, mass_intervals)

Sample *N* particle masses from a continuous, piecewise power-law distribution.

- **Parameters:**
  - `alphas` – Power-law exponents for each mass interval.  
  - `mass_intervals` – List of mass breakpoints, length = len(alphas) + 1.

- **Returns:**
  - Updates `self.masses` with sampled values.

- **Description:**
  - Ensures continuity between adjacent power-law segments.  
  - Normalizes the probability distribution.  
  - Divides random samples among intervals based on relative weights.  
  - Generates random masses in each segment using inverse sampling.  
  - Handles the special case of flat logarithmic distributions.

---

### **method: build_phasespace**(nDim, R0, Vmax=None, model='uniform')

Build positions and velocities for the N-body system.

- **Parameters:**
  - `nDim` – Number of spatial dimensions (2 or 3).  
  - `R0` – Characteristic radius (astropy quantity).  
  - `Vmax` – Maximum velocity (for uniform distribution only).  
  - `model` – Choice of distribution: `'uniform'` or `'plummer'`.

- **Returns:**
  - Updates `self.posVels` with concatenated positions and velocities.

- **Description:**
  - For `'uniform'`:  
    - Samples positions uniformly inside a solid sphere of radius *R0*.  
    - Samples velocities uniformly within a sphere of radius *Vmax*.  
  - For `'plummer'`:  
    - Samples positions and velocities following the Plummer distribution.  
    - Computes random directions for isotropy.  
    - Stores escape velocities for reference in `self.vesc`.  
  - Final output contains both positions and velocities for each particle.

---

**Summary:**

This module provides:
- Uniform and Plummer sphere samplers for initial conditions.  
- Tools to generate particle masses following realistic mass functions.  
- Simple kinetic and potential energy estimators for gravitational systems.
