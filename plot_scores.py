#!/usr/bin/env python
"""
Plot diagram of TAD border scores distribution. 
On the horizontal axis scores from 1 to 10 are put. 
For each score value a number of TAD borders with this score is calculated.

Usage:
  plot_scores.py -d <directory_with_score_files> -o <output_filename> [--header <plot_header>]

Options:
  -h --help                        Show this screen.
  --version                        Show version.
  -d <directory_with_score_files>  Directory with score files.
  -o <output_filename>             A name for the output file with a bar plot.
  --header <plot_header>           A header for the plot. Default: 'Number of TAD borders'
"""


import sys

print

modules = ["docopt", "os", "numpy", "matplotlib"]
exit_flag = False
for module in modules:
    try:
        __import__(module)
    except ImportError:
        exit_flag = True
        sys.stderr.write("Error: Python module " + module + " is not installed.\n")

if exit_flag:
    sys.stderr.write("You can install these modules with a command: pip install <module>\n")
    sys.stderr.write("(Administrator privileges may be required.)\n")
    sys.exit(1)


from docopt import docopt
from os.path import join
from os.path import exists
from os.path import isdir
from os import listdir
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., height + 2, '%d'%int(height),
                ha='center', va='bottom')


if __name__ == '__main__':
    arguments = docopt(__doc__, version='plot_score_distrib 0.3')
    input_directory = arguments["-d"].rstrip('/')
    if not exists(input_directory):
        print "Error: Can't find input directory: no such directory '" + \
              input_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(input_directory):
        print "Error: Input directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    output_filename = arguments["-o"]
    plot_header = arguments["--header"]
    if plot_header == '':
        plot_header = 'Number of TAD borders'
    
    score_distrib = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
    score_files = listdir(input_directory)
    for file in score_files:
        file_full = join(input_directory, file)
        curr_score = 0
        with open(file_full, 'r') as src:
            for line in src:
                line_list = line.rstrip('\n').split()
                if line_list[0] == '#':
                    continue
                curr_score = int(float(line_list[3]))
                score_distrib[curr_score] += 1
            score_distrib[curr_score] -= 1 # Exclude the last border score
    print 'Score distribution:     ', score_distrib
    total_count = sum(score_distrib.values())
    print 'Total number of borders:', total_count
    print

    fig = plt.figure()
    ax = fig.add_subplot(111)
    n = 10
    values = []
    for key, value in score_distrib.iteritems():
        values.append(value)
    ind = np.arange(n)
    width = 0.7
    y_max = max(values) + 100 #800
    bars = ax.bar(ind, values, width, color='blue')
    autolabel(bars)
    ax.set_ylim(0, y_max)
    ax.set_xlabel('Scores')
    ax.set_ylabel('Number of borders')
    ax.set_title(plot_header)
    legend_text = 'Total number of borders: ' + str(total_count)
    ax.legend([bars[0]], [legend_text], loc='upper left')
    xTickMarks = [str(i) for i in range(1, 11)]
    ax.set_xticks(ind + width / 2)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames)
    plt.savefig(output_filename)

