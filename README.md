READ ME

File Name: 2d_isng_4.cpp
Team: Omri, Max and Angela Chen


function name: Init_lattice
returns: void
Purpose: The purpose of this function is to generate a spin matrix with a random
         spin configuration.

function name: ising_MC
returns: delta_e which is the change in energy due to the accepted/rejected energy 
         change.
Purpose: The purpose of this function is to do a Metropolis procedure. We select
         a lattice position at random and compute the change in energy due to the
         new configuration. If it is accepted, we compute the running sum of the 
         expectation values.

function name: transient_effects
returns: void
Purpose: The purpose of this function is to equilibrate the lattice configration before
         applying the metropolis procedure onto it. 

function name: get_energy
returns: Energy which represents at a lattice position.
Purpose: A function to return the energy of the lattice due to its interaction with
         its neighboring spin, periodic boundary conditions where implemented into
         the algorithm.

function name: energy_total
returns: energy
Purpose: Function used to return the value of the total energy of the entire lattice.
         Nested loop in order to cycle through entire latitce. The functon get energy
         is called to return the energy at each postion. 

function  name: magnet_total()
returns: int m, representing total megnetization 
Purpose: Returns the magnetization of the entire lattice. It is computed by cycling
         through each element of the lattice and summing each spin.   
