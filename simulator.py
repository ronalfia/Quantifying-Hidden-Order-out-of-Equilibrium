from clg import create_clg_lattice, parallel_update, clg_activity
from manna import create_manna_lattice, parallel_manna_update, manna_activity
from compression import cid

import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import time

Z=2
#==============================================================================
# create_multiple_realizations_data(L,N,T,R)
# this function employs parallel computation to distributelly compute multiple
# realizations at once. The output data is stored in csv files named
# realization_{i} where i is a sequential running index
#==============================================================================
def create_multiple_realizations(L,N,T,R,model = 'clg',observable = 'cid'):
    cores = mp.cpu_count()
    # if the numbers of cores is not a divisor of R return an error
    if (R%cores!=0):
        print("the number of realizations must be a multiplication of the"\
              "number of cores which is {}".format(cores))
        return

    for batch in range(0,R,cores):
        processes = []
        for batch_realization in range(cores):
            processes.append(mp.Process(target=create_realization,\
                       args=(L,N,T,str(batch_realization+batch),model,observable)))
            processes[batch_realization].start()

        while True in [processes[batch_realization].is_alive() \
                       for batch_realization in range(cores)]:
            time.sleep(15)
        print("finished computation of {} realizations".format((batch+cores)))

#==============================================================================
# create_realization(L,N,T)
# create one simulation of L sites with particles in N where N is an array of
# positive integers with maximum possible integer not ascending L. T is an array
# of time steps for which the system is sampled by measuring its cid.
# it creates a dataframe file with each column matching the different number of
# particle options and rows for the respective time steps named
# realization{signature}
#==============================================================================
def create_realization(L,N,T,signature,model = 'clg',observable = 'cid'):
    if (max(N) > L) and (model == 'clg'):
        print("Cannot initailize a density of particles larger than 1")
        return
    data = pd.DataFrame(index = create_index(T))

    for n in N:
        if (model == 'clg'):
            lattice = create_clg_lattice(n,L)
            if (observable == 'cid'):
                values = [cid(''.join(str(int(x)) for x in lattice))]
            elif (observable == 'activity'):
                values = [clg_activity(lattice)]
            for t in T:
                parallel_update(lattice,t,True)
                if (observable == 'cid'):
                    values.append(cid(''.join(str(int(x)) for x in lattice)))
                elif (observable == 'activity'):
                    values.append(clg_activity(lattice))

            data[str(n)] = values
            data.to_csv("realization{}.csv".format(signature))
        elif (model == 'manna'):
            lattice = create_manna_lattice(n,L)
            if (observable == 'cid'):
                values = [cid(lattice,model)]
            elif (observable == 'activity'):
                values = [manna_activity(lattice, Z)]
            for t in T:
                parallel_manna_update(lattice,t, Z)
                if (observable == 'cid'):
                    values.append(cid(lattice, model))
                elif (observable == 'activity'):
                    values.append(manna_activity(lattice, Z))
            data[str(n)] = values
            data.to_csv("realization{}.csv".format(signature))
        print("finished calculating density {} for signature {}".format(float(n)/L, signature))
    return

def compare_dynamical_rules(timesteps, timestep, length):
    activity_parallel = []
    activity_random = []

    lattice_for_parallel = create_clg_lattice(int(length*0.6),length)
    lattice_for_random = create_clg_lattice(int(length*0.6),length)

    for t in range(timesteps):
        if (t%timestep == 0):
            activity_parallel.append(clg_activity(lattice_for_parallel))
            activity_random.append(clg_activity(lattice_for_random))
        parallel_update(lattice_for_parallel,1)
        parallel_update(lattice_for_random,1,True)

    return activity_parallel,activity_random

#==============================================================================
# create_index(T)
# receive an array of time steps T and returns an index later to be used in the
# dataset storing values according to those timesteps
#==============================================================================
def create_index(T):
    index = [0]
    for t_loc,t_value in enumerate(T):
        index.append(index[t_loc]+t_value)
    return index
