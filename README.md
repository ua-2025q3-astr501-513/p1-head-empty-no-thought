[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/nqfiwWTG)
[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-2e0aaae1b6195c2367325f4f02e2d04e9abb55f0b24a779b69b11b9e10269abc.svg)](https://classroom.github.com/online_ide?assignment_repo_id=20675659&assignment_repo_type=AssignmentRepo)
# ASTR 513 Mid-Term Project

Welcome to the repository of the **Mid-Term Project** for ASTR 513.
This project is worth **20 points** and gives you the opportunity to
apply the numerical techniques we have covered in the course so far.

You are recommended to work in teams of 3 to 6 students.
To form team, come up with a unique team name and put it in this
[GitHub Classroom link](https://classroom.github.com/a/nqfiwWTG).

## Documentation customization from the HENT Team

The Head Empty, No Thoughts™ team uses three different 3D N-body integrator to integrate toy-model stellar clusters based on a piecewise power-law IMF and a Plummer model. The N-body integrators are the following:
- direct Leap-Frog integrator for small clusters (~100 particles).
- Barnes-Hut algorithm for large clusters (~1000 particles).
- an implicit RK15 from the Rebound Python package for higher accuracy for small clusters (~100 particles).

Additionally, the presentation slides are available [here](https://docs.google.com/presentation/d/1UEkMNe21g3MoQCzrra4jgIt7j6J2d666Kll3J-W7jzA/edit?usp=sharing) (will be released after the presentation in-class is done). 

## Documentation for each file

The code themselves are contained in the `src/` folder. There are 3 main "modules" or "components":

1. `initialConditionsBuilder.py`: used to generate realistic cluster initial conditions for the other integrators to use. This script contains a class with methods to build arbitrary piecewise power-law IMFs and initial conditions in phase space that are either drawn from a virialized Plummer sphere or a uniform spherical distribution in position and velocity spaces. It also has methods to compute the kinetic and potential energy to check the virial ratio of the created cluster.

2. `nbody_python.py`: loads (with some edition; see `nbody_python_demo.ipynb` [here](https://github.com/ua-2025q3-astr501-513/p1-head-empty-no-thought/blob/main/nbody_python_demo.ipynb) for a more in-detail explanation of what's in there) the calculated initial condition, then integrate the system onwards, using either leapfrog or RK4 method. This is not done in the format of a package; instead users can download the script and run it themselves, out of the consideration for a lighter installation. The dependency in this module is only `numpy`. This part of the code referenced [5 steps to N-body simulation](https://alvinng4.github.io/grav_sim/5_steps_to_n_body_simulation/). Additionally, `make_gif.py` is the python script for generating gifs used for presentation, dependent on `matplotlib`. It should be used in conjunction with `bb_colors.txt` (from [m15y](https://temperature.m15y.com/)) for maximum ability - mainly, this adds the ability to make blackbody temperature colored stars in a gif, but I think it's cool.) 

3. (c code readme here...)

## Timeline & Deliverables

* Prsentation dates:
  October 15th.
* Submission deadline:
  By 11:59pm (Arizona tgime) on the day of your presentation
* Submission platform: GitHub Classroom

Your final submission should include:

* Project code (inside the `src/` directory of this git repository)
* Documentation (inside the `doc/` directory)
* Presentation materials (slides or Jupyter notebook, also version
  controlled with this git repository)

Only **one submission per team** is needed.

**Late submissions may not be accepted. Please plan ahead.**

## Project Ideas

The file `doc/ideas.yaml` contains a compilation of topics from
homework set \#1.
Please use this list to help you look for other students with similar
interests and form teams.

Example from a past project:
[Exoplanet Statistics](https://github.com/ua-2024q3-astr513/ASTRSTATS513_final).

## Requirements

### 1. Code

* Submit well-documented, runnable source code.
* Include docstrings and inline comments.
* Update this `README.md` file to explain:
  * How to install and run your project
  * Any dependencies or data required

### 2. Presentation

* Deliver a ~ 15 minute presentation on your project.
* You may use either:
  * Slides (traditional format), or
  * A Jupyter notebook (similar to our class style).

Your presentation should:
* clearly explain the problem you tackled;
* show the numerical techniques you applied;
* present results with relevant plots, tables; or figures
* highlight your findings and insights.

Each team member should be prepared to discuss their contributions.

## Grading (20 points total)

Projects will be graded based on the following criteria:
1. Originality & clarity of the idea
2. Quality of the solution (numerical methods, implementation,
   correctness)
3. Thoroughness of documentation (code comments, docstrings, README)
4. Effectiveness of presentation (clarity, structure, visualizations,
   teamwork)

## Collaboration & GitHub Use

Projects are managed through GitHub Classroom.
* Multiple students can share a single repository.
  You can join by putting your unique team name in this
  [GitHub Classroom link](https://classroom.github.com/a/nqfiwWTG).
* Use GitHub to track progress, manage code, and collaborate.
* Only one final submission per team is needed.

## Final Note

This project is your chance to be creative, apply what you have
learned so far, and work collaboratively on a meaningful computational
astrophysics problem.
We look forward to your results!
