from simulator import create_multiple_realizations
from data_analysis_tools import analyze_data,visualize_results

import pandas as pd
import numpy as np
import time

N = range(0,10001,500)
L = 10000
T = [150000]
cores = 8
model = 'clg'

def run_simulation(L,N,T,cores,model = 'clg'):
    start_time = time.clock()
    """ perform the simulation according to the given parameters

     This simulation creates a one dimensional lattice with periodic b.c i.e
     a ring of gridded sites. Each site can hold up to one particle and the
     entire setup will hold N different particles. Those particles will be
     propagated according to the clg rules. This simulation replicates the
     same external parameters for an ensemble of different realizations each
     with randomly choosen initial conditions

     Args:
         L (int): length of the lattice
         N (list): number of particles
         T (list): propagation times
         cores (int): the number of different realizations

     Returns:
         pandas_DataFrame: the results of the simulation where N are the
         columns, T dictates the rows (row i equals sum of first i elements in
         T) and each cell holds a list of the CID values for each realization
    """
    create_multiple_realizations(L,N,T,cores,model,'activity')

    data = pd.DataFrame(data = analyze_data(N,T,cores,True))

    visualize_results(L,N,T,cores,model,data,True,True)

    print("the simulation took {} seconds".format(time.clock()-start_time))

if __name__ == "__main__":

    run_simulation(L,N,T,cores,model)
