#### Chain CLG Model ####

import numpy as np
import random


##### Lattice Methods #####

#==============================================================================
# create_clg_lattice(n,L)
# Arguments - N is the number of particles
#             L is the number of sites
# returns an array with N particles randomly organized on L sites
#==============================================================================
def create_clg_lattice(n,L): #length of the dimension
    lattice = np.zeros(L)
    if n < L:
        for it in random.sample(range(L), n):
            lattice[it]+=1
    elif n==L:
        for it in range(L):
            lattice[it]+=1
    else:
        print("number of particles cant be bigger than the system's size")
    return lattice

#### Update Methods ####
#==============================================================================
# find_active_sites(lattice):
# a site is considered according to the following logic -
# it has to be occcupied
# it has to have atleast one empty neighbor and atleast one occupied neighbor
# returns an array of the active sites
#==============================================================================
def find_active_sites(lattice):
    L = len(lattice)
    active_sites = []
    for site in range(len(lattice)):
        if lattice[site]:
            if (lattice[(site+1)%L]==1.0) ^ (lattice[(site-1)%L]==1.0):
                active_sites.append(site)
    return active_sites

#==============================================================================
# find_empty_neighbor(lattice,site)
# returns the empty neighbor of the actice site
# assumes that the site is active
# if the empty neighbor is on the left of the active site return -1
# else, if the empty neighbor is on the right of the active site returns 1
#==============================================================================
def find_empty_neighbor(lattice,site):
    L = len(lattice)
    if lattice[(site+1)%L] == 1.0:
        return -1
    return 1

#==============================================================================
# fix_competition(active_particles, lattice)
# this function aims to resolve the case where two active particles are
# competiting for the same unoccupied neighbor. It will automatically resolve
# all such cases by randomly (with equal probabilities) deactivate one of the
# active competing sites
# although not utilized it will also return the resulting active sites
#==============================================================================
def fix_competition(active_particles,lattice):
    L_a = len(active_particles)
    L = len(lattice)
    active_particles_to_be_deactivated = set()
    for activity_index,particle_loc in enumerate(active_particles):
        if active_particles[(activity_index+1)%L_a] == particle_loc + 2:
            if (random.uniform(0,1) < 0.5):
                active_particles_to_be_deactivated.add(particle_loc)
            else:
                active_particles_to_be_deactivated.add(active_particles[(activity_index+1)%L_a])
    if (0 in active_particles and (L-2) in active_particles):
        if (random.uniform(0,1) < 0.5):
            active_particles_to_be_deactivated.add(0)
        else:
            active_particles_to_be_deactivated.add(L-2)
    if (1 in active_particles and L-1 in active_particles):
        if (random.uniform(0,1) < 0.5):
            active_particles_to_be_deactivated.add(1)
        else:
            active_particles_to_be_deactivated.add(L-1)
    for particle in active_particles_to_be_deactivated:
        active_particles.remove(particle)
    return active_particles

#==============================================================================
# parallel_update(lattice,timesteps)
# propages the lattice using parallel update for every active site
# if two active sites are competing over the same empty neighbor one will be
# randomly chosen
# will make a total number of timesteps updates to the lattice
#==============================================================================
def parallel_update(lattice,timesteps=1,randomize=False):
    L = len(lattice)
    for t in range(timesteps):
        active_sites = find_active_sites(lattice)

        if (len(active_sites) == 0):
            break
        if randomize:
            active_site = random.choice(active_sites)
            lattice[active_site] -= 1
            lattice[(active_site+(find_empty_neighbor(lattice,active_site)))%L] += 1

        else:
            fix_competition(active_sites,lattice)
            for active_site in active_sites:
                ## assumes competion is resolved -> displace active particles
                lattice[active_site] -= 1
                lattice[(active_site+(find_empty_neighbor(lattice,active_site)))%L] += 1

def clg_activity(lattice):
    return float(len(find_active_sites(lattice)))/len(lattice)
