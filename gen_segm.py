#!/usr/bin/env python
"""
Generate segmentation of reference genome into combinations of markers (states of ChromHMM).
BED file with windows of fixed size is added into input directory.
Testing is necessary.

Usage:
  gen_segm.py -d <input_directory> -s <chromosome_sizes> -i <states_number_interval> -c <ChromHMM_directory> -p <common_prefix> -o <output_directory> [-g <reference_genome_tag> -w <window_size> -t <thread_number>]

Options:
  -h --help                    Show this screen.
  --version                    Show version.
  -d <input_directory>         Directory with BED files containing peaks for all marks to be processed.
  -s <chromosome_sizes>        Text file with chromosome sizes.
  -i <states_number_interval>  Interval of states number (a:b = [a,b]; :b = [1,b]; a, b = 1, 2, ...; a < b < 2n - 2, where n is a number of markers.
  -c <ChromHMM_directory>      ChromHMM directory.
  -p <common_prefix>           Common prefix of all peak BED files.
  -o <output_directory>        Output directory name.
  -g <reference_genome_tag>    Tag of reference genome. Default: hg19.
  -w <window_size>             Size of the windows to split reference genome. Default: 200.
  -t <thread_number>           Number of threads to be used by ChromHMM. Default: 4.
"""

import sys

print

modules = ["docopt", "os"]
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
from os.path import basename
from os.path import splitext
from os.path import join
from os.path import exists
from os.path import isfile
from os.path import isdir
from os import makedirs
from os import listdir
from os import system
from os import devnull
from os import getcwd
from os.path import dirname
from os.path import abspath
from sys import stdout


def gen_windows_bed(window_size, chromosome_sizes_file, input_directory, genome_tag):
    windows_bed_file = join(input_directory, genome_tag + '_w' + str(window_size) + '_windows.bed')
    command_line_list = ['bedtools', 'makewindows', '-g', chromosome_sizes_file, \
                         '-w', str(window_size)]
    command_line = 'bedtools makewindows -g ' + chromosome_size_file + ' -w ' + \
                   str(window_size) + ' > ' + windows_bed_file
    with open(windows_bed_file, 'w') as dst:
        code_makewindows = call(command_line_list, stdout=dst)
    if code_makewindows != 0:
        print "Error: Can't make windows BED file. Exit.\n"
        sys.exit(1)
    return windows_bed_file


def gen_emission_table(input_directory, windows_bed_file, chromosome_size_file, prefix):
    tadstates_directory = dirname(abspath(__file__))
    make_table_path = join(tadstate_directory, 'make_table.sh')
    code_make_table = call([make_table_path, prefix, windows_bed_file, chromosome_size_file])
    if code_make_table != 0:
        print "Error: Can't generate table. Exit.\n"
        sys.exit(1)

def gen_segmentation(input_directory, output_directory, chromhmm_directory, thread_number, \
                     start_num, end_num, genome_tag):
    for states_num in range(start_num, end_num + 1):
        out_directory = join(output_directory, 'Segm_' + str(states_num) + '_states')
        command_line_list = ['java', '-mx4000M', '-jar', \
                             join(chromhmm_directory, 'ChromHMM.jar'), 'LearnModel', '-p', \
                             str(thread_number), input_directory, out_directory, \
                             str(states_num), genome_tag]
        command_line = 'java -mx4000M -jar ' + join(chromhmm_directory, 'ChromHMM.jar') + \
                       ' LearnModel -p ' + str(thread_number) + ' ' + input_directory + \
                       ' ' + out_directory + ' ' + str(states_num) + ' ' + genome_tag
        print 'For ' + states_num + '...'
        print command_line
        stdout.flush()
        code_chromhmm = call(command_line_list)
        if code_chromhmm != 0:
            print "Error: Can't generate segmentation for " + str(states_num) + "states. "
            stdout.flush()


if __name__ == '__main__':
    arguments = docopt(__doc__, version='gen_segm 0.1')
    input_directory = arguments["-d"].rstrip('/')
    if not exists(input_directory):
        print "Error: Can't find input directory: no such directory '" + \
              input_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(input_directory):
        print "Error: Input directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    chromosome_size_file = arguments["-c"]
    if not exists(chromosome_size_file):
        print "Error: Can't find text file with chromosome sizes: no such file '" + \
              chromosome_size_file + "'. Exit.\n"
        sys.exit(1)
    if not isfile(chromosome_size_file):
        print "Error: File with chromosome sizes must be a regular file. " + \
              "Something else given. Exit.\n"
        sys.exit(1)

    state_interval = arguments["-i"]
    interval_list = state_interval.split(':')
    if len(interval_list) <= 1:
        print "Error: Interval of state numbers should contain ':'. Exit.\n"
        sys.exit(1)

    start_num_str = interval_list[0]
    end_num_str = interval_list[1]

    if start_num_str == '':
        start_num = 1
    else:
        try:
            start_num = int(start_num_str)
        except ValueError:
            print "Error: Start number of states must be an integer from 1 to 2n-2, " + \
                  "where n is a number of marks. Exit.\n"
            sys.exit(1)
        if start_num < 1 or start_num > 2 * mark_num - 2:
            print "Error: Start number of states must be an integer from 1 to 2n-2, " + \
                  "where n is a number of marks. Exit.\n"
            sys.exit(1)

    if end_num_str == '':
         print "Error: End number of states must be defined. Exit."
         sys.exit(1)
    else:
        try:
            end_num = int(end_num_str)
        except ValueError:
            print "Error: End number of states must be an integer from 1 to 2n-2, " + \
                  "where n is a number of marks. Exit.\n"
            sys.exit(1)
        if end_num < 1 or end_num > 2 * mark_num - 2:
            print "Error: End number of states must be an integer from 1 to 2n-2, " + \
                  "where n is a number of marks. Exit.\n"
            sys.exit(1)
        if end_num < start_num:
            print "Error: End number of states must be greater than start number. Exit.\n"
            sys.exit(1)

    chromhmm_directory = arguments["-c"].rstrip('/')
    if not exists(chromhmm_directory):
        print "Error: Can't find ChromHMM directory: no such directory '" + \
              chromhmm_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(chromhmm_directory):
        print "Error: ChromHMM directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    prefix = arguments["-p"]

    output_directory = arguments["-o"].rstrip('/')
    if exists(output_directory) and not isdir(output_directory):
        print "Error: Output directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)


    if arguments["-g"] != None:
        genome_tag = arguments["-g"]
    else:
        genome_tag = 'hg19'

    if arguments["-w"] != None:
        window_size_str = arguments["-w"]
        try:
            window_size = int(window_size_str)
        except ValueError:
            print "Error: Window size must be an integer greater than 0. Exit.\n"
            sys.exit(1)
        if window_size < 0:
            print "Error: Window size must be an integer greater than 0. Exit.\n"
            sys.exit(1)
    else:
        window_size = 200

    if arguments["-t"] != None:
        thread_number_str = arguments["-t"]
        try:
            thread_number = int(thread_number_str)
        except:
            print "Error: Thread number must be an integer greater than 0. Exit.\n"
            sys.exit(1)
    else:
        thread_number = 4

    print 'Check bedtools...'
    stdout.flush()
    with open(devnull, 'w') as FNULL:
        check_bedtools_code = call(['bedtools', '--help'], stdout=FNULL)
    if check_bedtools_code != 0:
        print 'Error: bedtools are not found. (It must be installed and available for ' + \
              "non-administrative user with command 'bedtools'). Exit."
        sys.exit(1)
    print 'OK'
    stdout.flush()

    if not exists(output_directory):
        makedirs(output_directory)

    print
    print 'Generate BED file with ' + str(window_size) + ' bp windows...'
    stdout.flush()
    windows_bed_file = gen_windows_bed(window_size, chromosome_size_file, input_directory, \
                                       genome_tag)
    print 'OK.'
    print 'Generate emission table...'
    stdout.flush()
    gen_emission_table(input_directory, windows_bed_file, chromosome_size_file, prefix)
    print 'OK.'
    print 'Generate segmenation...'
    stdout.flush()
    gen_segmentation(input_directory, output_directory, chromhmm_directory, thread_number, \
                     start_num, end_num, genome_tag)
    print 'OK.'
    stdout.flush()

    print
    print "Done!"
    print
    stdout.flush()
