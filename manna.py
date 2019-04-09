#### Chain CLG Model ####

import numpy as np
import random


##### Lattice Methods #####


def create_manna_lattice(n,L):
    ''' creates a manna lattice of L sites and N particles

    This function creates a manna lattice with particles randomly distributed
    among its sites. The distribution process is implement through randomly
    choosing a site and adding it a particle for a total number of N (number of
    particles) times.

    Args:
        L (int): number of the lattice sites
        N (int): number of particles


    Returns:
        numpy array of a manna lattice with N randomly distributed particles

    '''
    lattice = np.zeros(L)
    for particle in range(n):
        particle_loc = random.randint(0,L-1)
        lattice[particle_loc] += 1
    return lattice

def count_particles(lattice):
    ''' counts the number particles in a lattice

    this function counts the number of particles in a lattice by interating
    over the lattice and counting the number of particles at each site, summing
    them all and returning the answer

    Args:
        lattice (numpy array): manna Lattice

    Returns:
        int the number of particles in the lattice

    '''
    n = 0
    for site in range(len(lattice)):
        n += lattice[site]
    return int(n)


def find_active_sites(lattice,z = 0):
    ''' finds the active sites of a manna lattice

    This function will find and return the active sites of a manna lattice
    according to the threshold value Z.
    Args:
        l (numpy array): the manna lattice
        Z (int): threshold value for the activity of a site

    Returns:
        list of the active sites for the provided manna lattice

    '''
    if (z == 0):
        print("Z is the threshold value for activity and it cannot be less than 1")
        return

    L = len(lattice)
    active_sites = []
    for site in range(len(lattice)):
        if lattice[site] > z:
            active_sites.append(site)
    return active_sites


def parallel_manna_update(lattice,timesteps=1, z = 0):
    ''' updates all the lattice sites in parallel

    This function will update the lattice by displacing the active particles in
    each active site randomly between its nearest neighbors. The timesteps
    argument dictates how many updates will be performed on the lattice. It is
    important to note that the updates take place on the provided lattice and
    not on a copy of it
    Args:
        lattice (numpy array): the manna lattice
        timesteps (int): the number of updates to be performed
        Z (int) :  threshold value for the activity of a site

    Returns:
        None

    '''
    if (z == 0):
        print("Z is the threshold value for activity and it cannot be less than 1")
        return

    L = len(lattice)

    for t in range(timesteps):
        active_sites = find_active_sites(lattice, z)
        propagated_lattice = lattice

        if (len(active_sites) == 0):
            break
        for active_site in active_sites:
            particles_to_the_right = random.randint(0,lattice[active_site]- z)
            particles_to_the_left = lattice[active_site] - z - particles_to_the_right

            lattice[active_site] -= particles_to_the_right + particles_to_the_left
            lattice[(active_site+1)%L] += particles_to_the_right
            lattice[(active_site-1)%L] += particles_to_the_left

def manna_activity(lattice,Z):
    ''' returns the activity of the manna lattice as a fraction of its active
    sites.

    Args:
        lattice (numpy array) : the manna lattice
    Returns:
        the density of active sites - float
    '''
    return float(len(find_active_sites(lattice,Z)))/len(lattice)
