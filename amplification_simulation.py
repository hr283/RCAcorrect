from __future__ import print_function
import numpy as np
import random
from pandas import DataFrame
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def simulation():
    orig_vec_length = 1000
    boundaries = [10, 40, 90, 180, 360, 600, 1000] # boundaries must be of length 7 (or edit the if clauses below!!
    record_vec = list(np.arange(0,orig_vec_length))
    amplify_n = 200 #this number of elements will be chosen each time to be amplified.
    rounds = 10
    percent_sequence = 20 #this number of reads will be chosen for sequencing

    # amplify
    r = 0
    vec_length = orig_vec_length
    while r < rounds:
        for elem in record_vec:
            rand = random.randint(1,int(vec_length/amplify_n))
            if rand == 1:
                record_vec.append(elem)

        vec_length = len(record_vec)
        r += 1

    # sequence
    reads = []
    for elem in record_vec:
        rand = random.randint(1,int(vec_length/percent_sequence))
        if rand == 1:
            reads.append(elem)

    # count
    orig_counts = np.zeros(len(boundaries))
    for elem in np.arange(0, orig_vec_length):
        if elem < boundaries[0]:
            orig_counts[0] += 1
        elif elem < boundaries[1]:
            orig_counts[1] += 1
        elif elem < boundaries[2]:
            orig_counts[2] += 1
        elif elem < boundaries[3]:
            orig_counts[3] += 1
        elif elem < boundaries[4]:
            orig_counts[4] += 1
        elif elem < boundaries[5]:
            orig_counts[5] += 1
        elif elem < boundaries[6]:
            orig_counts[6] += 1

    read_counts = np.zeros(len(boundaries))
    for elem in reads:
        if elem < boundaries[0]:
            read_counts[0] += 1
        elif elem < boundaries[1]:
            read_counts[1] += 1
        elif elem < boundaries[2]:
            read_counts[2] += 1
        elif elem < boundaries[3]:
            read_counts[3] += 1
        elif elem < boundaries[4]:
            read_counts[4] += 1
        elif elem < boundaries[5]:
            read_counts[5] += 1
        elif elem < boundaries[6]:
            read_counts[6] += 1

    data_to_plot = []
    for i, c in enumerate(orig_counts):
        orig_freq = float(c)/orig_vec_length
        read_freq = float(read_counts[i])/len(reads)
        data_to_plot.append({'True freq': orig_freq, 'Amplified freq': read_freq})

    return(DataFrame(data_to_plot))

# plotting
if __name__ == "__main__":
    dirpath = "/Users/hroberts/Desktop"
    random.seed(1)
    reps = 500
    left, width = 0.1, 0.8
    bottom, height = 0.1, 0.8
    rect_scatter = [left, bottom, width, height]

    fig = plt.figure(figsize=(6, 6))
    axScatter = plt.axes(rect_scatter)
    plt.title('Results of amplification simulations')

    amplified_freqs = []
    for i in range(reps):
        df = simulation()
        axScatter.plot( df['True freq'], df['Amplified freq'], '.', ms=4, mew=0, rasterized=False, color='purple' )
        amplified_freqs.append(df['Amplified freq'])

    amplified_freq_array = np.array(amplified_freqs)
    medians = [np.median(amplified_freq_array[:,i]) for i in range(np.shape(amplified_freq_array)[1])]
    axScatter.plot( df['True freq'], medians, '.', ms=6, mew=0, rasterized=False, color='orange')
    axScatter.plot([0, 0.5], [0, 0.5], 'k-', color = 'r')

    axScatter.set_xlabel('True freq')
    axScatter.set_ylabel('Amplified freq')
    axScatter.set_xlim(0, 0.5)
    axScatter.set_ylim(0, 0.5)

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['FreeSans']

    fig.savefig(dirpath + '/amplification_simulation.pdf')