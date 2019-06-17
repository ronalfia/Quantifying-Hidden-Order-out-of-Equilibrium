import numpy as np
from manna import create_manna_lattice,count_particles

def flatten_manna_configuration(configuration):
    ''' this function will create a lempel ziv compressible manna configuration

    states in a manna site can be multiple digits in a 10 base numerical system
    which makes their compression  be somewhat problematic. In order to avoid
    this difficulty this function will take any such number and transform it
    into the biggest possible one digit integer

    Args:
        configuration (numpy array): the microstate of a system

    Returns:
        numpy array of an easily compressed configuration
    '''
    flatten_configuration = configuration
    for site in range(len(configuration)):
        if configuration[site] > 9:
            flatten_configuration[site] = 9
    return flatten_configuration


def lz_78(configuration, model = 'clg'):
    ''' finds the different patterns of a configuration going from left to right

    returns the different patterns encountered in configuration while going from
    left to right. Assumes that the configuration is given as a string and was
    tests for configurations which are composed only of 0 or 1 values for each
    character. if there is an extra pattern while reaching the end of the string
    it will be counted as a different pattern

    Args:
        configuration (numpy array): the microstate of a system
        model (string): the name of the model implemented on that system

    Returns:
        list of the different patterns in that microstate
    '''
    if model == 'clg':
        if (type(configuration) == str):
            str_representation = configuration
        else:
            print("ERROR : Configuration Must be passed as a string to the compression algorithm")
            return
    elif model == 'manna':
        flatten_configuration = flatten_manna_configuration(configuration)
        str_representation = ''.join(str(int(x)) for x in flatten_configuration)
    seen_patterns = []

    current_pattern = ''

    for site in range(len(str_representation)):
        current_pattern += str_representation[site]
        if not (current_pattern in seen_patterns):
            seen_patterns.append(current_pattern)
            current_pattern = ''
    if not (current_pattern == ''):
        seen_patterns.append(current_pattern)
    return seen_patterns


    #=============================================================
# lz_78_number_of_patterns(configuration)
# returns the number of different patterns according to the lz_78 function
#==============================================================================
def lz_78_number_of_patterns(configuration,model = 'clg'):
    return len(lz_78(configuration,model))


def cid(configuration, model='clg', random_shuffle=False):
    ''' computes the cid of a configuration

    implemented for two models namely conserved lattice gas and manna model,
    this function computes the Computable Information Density of a configuration
    by using counting the number of different patterns and inserting them to an
    entropy operator and normalizing by a random reference

    Args:
        configuration (numpy array): the microstate of a system
        model (string): the name of the model implemented on that system

    Returns:
        float of the Computable Information Density for that microstate

    '''
    n_p, random_reference = 0.0, 0.0
    if model == 'clg':
        n_p = lz_78_number_of_patterns(configuration)
        random_reference = lz_78_number_of_patterns(''.join(str(x) for x in np.random.randint(2,size=len(configuration))))
        cid = (n_p*np.log2(n_p))/(random_reference*np.log2(random_reference))
        return cid
    elif model == 'manna':
        n_p = lz_78_number_of_patterns(configuration,'manna')
        if random_shuffle:
            random_reference = np.array(configuration)
            np.random.shuffle(random_reference)
            random_reference_str = ''.join(str(int(x)) for x in random_reference)
            random_reference = lz_78_number_of_patterns(random_reference_str)
            return (n_p * np.log2(n_p)) / (random_reference*np.log2(random_reference))
        else:
            return (n_p * np.log2(n_p)) / len(configuration)
    else:
        print("a compression scheme for {} has not been implemented".format(model))
        return
