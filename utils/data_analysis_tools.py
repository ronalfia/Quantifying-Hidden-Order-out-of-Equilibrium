import pandas as pd
import os
from main.simulator import create_index
import numpy as np
import matplotlib.pyplot as plt

def string_list(l):
    return '{'+ str(l[0]) +',' + "".join(str(x)+',' for x in l[1:len(l)-1])+ str(l[len(l)-1]) +'}'

def create_figure_caption(L,N,T,R, model = 'clg'):
    caption = string_list(N) + ' Particles \nhave been distributed on a discrete ring with'
    caption += ' {} sites. Those particles have been randomly distributed '.format(L)
    caption += 'for {} separated times.\nEach distribution has been individually'.format(R)
    if model == 'clg':
        caption += ' propagated according to the CLG rules to a total number '
    elif model == 'manna':
        caption += ' propagated according to the Manna rules to a total number '
    else:
        print("the model {} has not been implemented".format(model))
        return
    caption += 'of {} updates.\nThe CiD has been measured and plotted'.format(sum(T))
    caption += ' after the following time intervals - ' + string_list(T)
    caption += ' The different\nrealizations of each unique density enabled'
    caption += ' the computation of the standard deviation for each density. It was\n'
    caption += 'plotted as the error bar of that density.'
    return ""

def visualize_results(L, N, T, cores, model, data, y_caption_height=-0.4, analytical_result=False, activity=False):
    average_values = data.applymap(lambda x : sum(x)/len(x))
    standard_deviation = data.applymap(lambda x : np.std(x))

    rect = 0.25,0.25,0.5,0.5
    clg_figure = plt.figure(figsize=[12,10],facecolor='white')
    clg_figure.add_axes(rect,facecolor='grey')

    x_values = []

    for row_index,row in enumerate(average_values.iterrows()):
        x_values = [float(x) / L for x in list(pd.Series(row[1]).index)]
        y_values = pd.Series(row[1]).values
        yerr_values = pd.Series(standard_deviation.iloc[row_index]).values
        label_ = str(row[0])
        plt.errorbar(x = x_values,y = y_values, yerr = yerr_values, label = label_,fmt = 'o-')

    if analytical_result:
        plt.plot(x_values[len(x_values)/2:],[(2/x)*(2*x-1)*(1-x) for x in x_values[len(x_values)/2:]])

    plt.text(-0.25,y_caption_height,create_figure_caption(L,N,T,cores,model))
    plt.grid(axis='y',linewidth = 2)
    plt.grid(axis='x',linewidth = 1)
    if model == 'manna':
        plt.title('Manna Model on a Chain')
    else:
        plt.title('Conserved Lattice Gas on a Chain')
    plt.legend(title = 'Time Steps', shadow=True ,fontsize='large')
    plt.xlabel('Density')
    if activity:
        plt.ylabel('Activity')
    else:
        plt.ylabel('Computable Information Density')
    plt.savefig(str(model)+'_figure.png')

def analyze_data(N,T,R,remove=False):
    # concat realizations into one ensemble table
    data_table = insert_data(R,create_data_table(N,T))

    # save ensemble table
    data_table.to_csv('ensemble of {} realization.csv'.format(R))

    # remove realizations files
    if remove:
        remove_realization_files(R)

    return data_table

def create_data_table(N,T):
    data_table = pd.DataFrame(index=create_index(T),columns=N)
    for index in data_table.index:
        for column in data_table.columns:
            data_table[column][index] = []
    return data_table


def insert_data(R,data_table):
    for realization in ['realization'+str(r)+'.csv' for r in range(R)]:
        realization_table = pd.read_csv(realization)
        for index_realization,index_data in enumerate(data_table.index):
            for column in data_table.columns:
                data_table[column][index_data].append(realization_table[unicode(str(column),"utf-8")][index_realization])
    return data_table

def remove_realization_files(R):
    for realization in ['realization'+str(r)+'.csv' for r in range(R)]:
        os.remove(realization)
