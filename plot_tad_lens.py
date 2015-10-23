#!/usr/bin/env python

"""
Usage:
  plot_tad_lens.py -i <tad_bedfiles_list> [--header <plot_header> -l <labels_list> -o <output_filename> --no-fliers --show-means] 

Options:
  -h --help               Show this screen.
  --version               Show version.
  -i <tad_bedfiles_list>  One TADs BED filename or several TADs BED filenames separated by comma.
  --header <plot_header>  The header for the plot. Default: 'TADs length boxplot'.
  -l <lables_list>        The list of labels for boxes on the plot. Default: '1','2','3',... .
  -o <output_filename>    The output file for the boxplot. Default: 'TAD_lens_boxplot.png'.
  --no-fliers             Don't plot outliers.
  --show-means            Plot mean values.
"""

import sys

print

modules = ["docopt", "os", "matplotlib"]
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
from os.path import exists
from os.path import isfile
from sys import stdout
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_tad_lens(bed_filenames, plot_header, tad_labels, output_filename, no_fliers, show_means):
    print "Create TAD length boxplot for the following files:"
    stdout.flush()
    for bed_filename in bed_filenames:
        print bed_filename
    print
    stdout.flush()
    print "Plot header:", plot_header
    print "Plot outliers:", 'Yes' if not no_fliers else 'No'
    print "Plot means:", 'Yes' if show_means else 'No'
    print
    stdout.flush()
    print "TAD labels:"
    stdout.flush()
    for tad_label in tad_labels:
        print tad_label
    print
    stdout.flush()
    print "Output file:", output_filename
    print
    stdout.flush()

    # Plot TADs length boxplot
    print 'Plot TADs length boxplot...',
    stdout.flush()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    boxplot_data = []
    for bed_filename in bed_filenames:
        with open(bed_filename, 'r') as bed_file:
            tad_lens = []
            for line_number, line in enumerate(bed_file):
                if line_number == 0:
                    continue
                bed_fields = line.rstrip('\n').split('\t')
                start_coord = int(bed_fields[1])
                end_coord = int(bed_fields[2])
                tad_lens.append(end_coord - start_coord)
            boxplot_data.append(tad_lens[:])
            del(tad_lens[:])
    sym = '' if no_fliers else 'b.'
    ax.boxplot(boxplot_data, labels = tad_labels, showmeans = show_means, sym = sym, whis = [5, 95])
    ax.set_ylabel('TAD length, bp')
    ax.set_title(plot_header)
    plt.savefig(output_filename)
    print 'Finish.'
    stdout.flush()


if __name__ == '__main__':
    arguments = docopt(__doc__, version='plot_tad_lens 0.3')
    bed_filenames = arguments["-i"].split(",")
    for bed_filename in bed_filenames:
        if not exists(bed_filename):
            print "Error: Can't find TAD BED file: no such file '" + \
                  bed_filename + "'. Exit.\n"
            sys.exit(1)
        if not isfile(bed_filename):
            print "Error: TAD BED file must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
     
    if arguments["--header"] != None:
        plot_header = arguments["--header"]
    else:
        plot_header = 'TADs length boxplot'

    if arguments["-l"] != None:
        tad_labels = arguments["-l"].split(",")
        for tad_label in tad_labels:
            if tad_label == '':
                print "Error: labels must be non-empty strings. Exit."
                sys.exit(1)
    else:
        tad_labels = range(len(bed_filenames))

    if arguments["-o"] != None:
        output_filename = arguments["-o"]
    else:
        output_filename = 'TAD_lens_boxplot.png'

    if arguments["--no-fliers"]:
        no_fliers = True
    else:
        no_fliers = False

    if arguments["--show-means"]:
        show_means = True
    else:
        show_means = False

    plot_tad_lens(bed_filenames, plot_header, tad_labels, output_filename, no_fliers, show_means)
    print 'Processing is finished.'

