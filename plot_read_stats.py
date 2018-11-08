from __future__ import print_function
import pysam
import argparse
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pandas import DataFrame
import matplotlib.patches as mpatches

"""
This script was used to make the read vs template length plot in Figure 3.
"""

def get_data_frame(bam_file_name, genome_length):
    data = []
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    for read in samfile.fetch():
        if not read.is_secondary and not read.is_supplementary and not read.is_unmapped:
            total_length=read.query_length
            mapped_length=read.query_alignment_length
            if mapped_length > genome_length:
                mapped_length = genome_length
            data.append({'Read length': total_length, 'Template length': mapped_length})

    return DataFrame(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('primary_bam', type=str, help="path to bamfile with primary mappings (to a concatenated ref)")
    parser.add_argument('dirpath', type=str, help="path to output directory")
    args = parser.parse_args()

    genome_length = 3300
    df = get_data_frame(args.primary_bam, genome_length)

    # classify reads as full length (containing a whole genome copy) or partial length, with some leniency.
    full_length = df['Template length'] >= genome_length - 100
    partial_length = df['Template length'] < genome_length - 100

    nullfmt = ticker.NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(figsize=(6, 6))
    axScatter = plt.axes(rect_scatter)
    plt.yscale('log')
    axHistx = plt.axes(rect_histx)
    plt.title('Read length vs template length')
    axHisty = plt.axes(rect_histy)
    plt.yscale('log')

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot
    axScatter.plot(df[full_length]['Template length'], df[full_length]['Read length'], '.', ms=4, mew=0, rasterized=True,
                   color='purple', label='Full')
    axScatter.plot(df[partial_length]['Template length'], df[partial_length]['Read length'], '.', ms=4, mew=0, rasterized=True,
                   color='rosybrown', label='Partial')

    axScatter.set_xlabel('Template length (bp)')
    axScatter.set_ylabel('Read length (bp)')
    axScatter.set_xlim(0, genome_length + 100)
    axScatter.set_ylim(100, 100000)

    full_data = mpatches.Patch(facecolor='purple', label='Full', edgecolor='black')
    partial_data = mpatches.Patch(facecolor='rosybrown', label='Partial', edgecolor='black')

    axScatter.legend(handles=[full_data, partial_data], loc='best')

    # the histograms
    xbins = np.linspace(0, genome_length + 100, 61)
    ybins = np.logspace(np.log10(100), np.log10(200000), 61)
    axHistx.hist([df[full_length]['Template length'], df[partial_length]['Template length']], bins=xbins, rwidth=100,
                 color=['purple', 'rosybrown'], stacked=True)
    axHisty.hist([df[full_length]['Read length'], df[partial_length]['Read length']], bins=ybins, rwidth=100,
                 orientation='horizontal', color=['purple', 'rosybrown'], stacked=True)

    axHistx.locator_params(nbins=4, axis='y')
    axHisty.locator_params(nbins=4, axis='x')

    axHistx.set_xlim(0, genome_length + 100)
    axHistx.set_ylim(0, 1000)
    axHisty.set_ylim(100, 100000)

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['FreeSans']

    fig.savefig(args.dirpath + '/template_vs_read_length.pdf')

