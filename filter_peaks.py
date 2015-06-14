#!/usr/bin/env python
"""
Filter peaks by score value using config file.

Usage:
  filter_peaks.py -d <input_directory> -c <config_file> -o <output_directory>

Options:
  -h --help              Show this screen.
  --version              Show version.
  -d <input_directory>   Directory with BED files containing peaks for all marks to be processed.
  -c <config_file>       Text file with chosen calls and score thresholds.
  -o <output_directory>  Output directory name.
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
from sys import stdout


def filter_peaks(bed_file_full, output_directory, score_thr, caller_name):
    if caller_name == 'SICER':
        bed_file_full_filtered = (splitext(bed_file_full))[0] + (splitext(bed_file_full))[1] + \
                                '_greater' + str(int(score_thr)) + '.bed'
    else:
        bed_file_full_filtered = (splitext(bed_file_full))[0] + '_greater' + \
                                 str(int(score_thr)) + '.bed'
    with open(bed_file_full_filtered, 'w') as dst, open(bed_file_full, 'r') as src:
        for line in src:
            fields_list = line.rstrip('\n').split()
            if caller_name == 'SICER':
                score_str = fields_list[3]
            else:
                score_str = fields_list[4]
            score = float(score_str)
            if score > score_thr:
                dst.write(line)


if __name__ == '__main__':
    arguments = docopt(__doc__, version='filter_peaks 0.1')
    input_directory = arguments["-d"].rstrip('/')
    if not exists(input_directory):
        print "Error: Can't find input directory: no such directory '" + \
              input_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(input_directory):
        print "Error: Input directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    config_file = arguments["-c"]
    if not exists(config_file):
        print "Error: Can't find config file: no such file '" + \
              config_file + "'. Exit.\n"
        sys.exit(1)
    if not isfile(config_file):
        print "Error: Config file must be a regular file. Something else given. Exit.\n"
        sys.exit(1)

    output_directory = arguments["-o"].rstrip('/')
    if exists(output_directory) and not isdir(output_directory):
        print "Error: Output directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    if not exists(output_directory):
        makedirs(output_directory)

    with open(config_file, 'r') as src:
        for line in src:
            field_list = line.rstrip('\n').split()
            mark_name = field_list[0]
            caller_name = field_list[1]
            if caller_name == 'SICER':
                gap_size = int(field_list[2])
                score_thr_str = field_list[3]
            elif caller_name == 'MACS':
                score_thr_str = field_list[2]
            else:
                print "Error: Undefined caller name found in the config file: " + caller_name + \
                      ". It must be SICER or MACS. Exit.\n"
                sys.exit(1)
            if score_thr_str != '-':
                score_thr = float(score_thr_str)
            else:
                score_thr = 0
            for file in listdir(input_directory):
                if mark_name in file:
                    if caller_name == 'SICER':
                        if 'scoreisland' in file and 'G' + str(gap_size) in file:
                            bed_file_full = join(input_directory, file)
                            break
                    else: 
                        if 'MACS' in file:
                            bed_file_full = join(input_directory, file)
                            break
            print
            if score_thr == 0:
                print 'Pass ' + bed_file_full
                stdout.flush()
                continue
            print 'Filter ' + bed_file_full + ' by score >= ' + str(score_thr) + ' ...'
            stdout.flush()
            filter_peaks(bed_file_full, output_directory, score_thr, caller_name)
