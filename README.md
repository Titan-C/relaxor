Monte Carlo Simulation of a Relaxor Ferroelectric
=================================================

Contains the code for an Ising-like model of a relaxor
ferroelectric, where nodes represent polar nano regions
and are arranged in a cubic lattice. Interaction between
nodes follows a Gaussian distribution. Includes the
interaction with an external AC electric field.

Much of the code is written in Spanish, a translation
of comments, and other necessary stuff is still due, as
well as documentation.

To compile the code just run make, you will need GSL
Alternatively use cmake. Call it from a build directory
or use kdevelop.
The later method is prefered because it also builds a python
module. Which allows to call the simulation from python
scripts and use the main.py for simulation needs

Additional tools required are Python, SciPy, and pyeq2
