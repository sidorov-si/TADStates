#!/usr/bin/env python
"""
Generate BED files with TAD borders of particular width.
TAD border of width 2r between two TADs is defined as a region consisting of 
r bp to the left and r bp to the right of the point separating the two TADs.

Usage:
  call_borders.py (-f <TADs_filename> | -d <input_directory>) -w <border_width> [-n <track_name_for_all_borders> -o <output_directory>]

Options:
  -h --help                        Show this screen.
  --version                        Show version.
  -f <TADs_filename>               Name of a BED file with TAD coordinates.
  -d <input_directory>             Name of a directory with BED files containing TAD coordinates.
  -w <border_width>                Border width (in bp).
  -n <track_name_for_all_borders>  A name for a track with all borders. Default: All_borders_<border_width>.
  -o <output_directory>            Output directory.
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
from os.path import exists
from os.path import join
from os.path import splitext
from os.path import basename
from os.path import isdir
from os.path import isfile
from os import listdir
from os import makedirs


def call_borders(filename, output_directory, half_width):
    name, ext = splitext(filename)
    file_basename = basename(name)
    if output_directory == '':
        file_name = name
    else:
        file_name = file_basename
    output_filename = join(output_directory, file_name + '_borders_' + str(border_width) + '.bed')
    
    prev_chrom_name = ''
    with open(filename, 'r') as src, open(output_filename, 'w') as dst:
        left_coord = ''
        i = -1
        for line in src:
            i += 1
            line_list = line.rstrip('\n').split('\t')
            if len(line_list) == 1:
                dst.write(line) # copy track line
                i -= 1
                continue
            if left_coord == '': # for the first line
                left_coord = int(line_list[2]) - half_width
                continue
            chrom_name = line_list[0]
            if chrom_name != prev_chrom_name and prev_chrom_name != '':
                left_coord = int(line_list[2]) - half_width
                i = 0
                prev_chrom_name = chrom_name
                continue
            right_coord = int(line_list[1]) + half_width
            border_name = chrom_name + '.border.' + str(i)
            score = 0 # Just to fill in the field
            strand = '.' # Just to fill in the field
            color = '0,255,0' # green
            border_line = chrom_name + '\t' + str(left_coord) + '\t' + str(right_coord) + '\t' + \
                          border_name + '\t' + str(score) + '\t' + strand + '\t' + \
                          str(left_coord) + '\t' + str(right_coord) + '\t' + color
            dst.write(border_line + '\n')
            left_coord = int(line_list[2]) - half_width
            prev_chrom_name = chrom_name


if __name__ == '__main__':
    arguments = docopt(__doc__, version='call_borders 0.3')
    if arguments["-f"] != None:
        filename = arguments["-f"]
        if not exists(filename):
            print "Error: Can't find BED file with TAD coordinates: no such file '" + \
                  filename + "'. Exit.\n"
            sys.exit(1)
        if not isfile(filename):
            print "Error: BED file with TAD coordinates must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
        input_directory = None
    else:
        filename = None
        input_directory = arguments["-d"].rstrip('/')
        if not exists(input_directory):
            print "Error: Can't find input directory: no such directory '" + \
                  input_directory + "'. Exit.\n"
            sys.exit(1)
        if not isdir(input_directory):
            print "Error: Input directory must be a directory:). Something else given. Exit.\n"
            sys.exit(1)

    try:
        border_width = int(arguments["-w"])
    except ValueError:
        print "Error: Border width must be an integer greater than 0. Exit.\n"
        sys.exit(1)
    if border_width <= 0:
        print "Error: Border width must be an integer greater than 0. Exit.\n"
    half_width = border_width / 2

    if arguments["-n"] != None:
        track_name = arguments["-n"]
    else:
        track_name = "All_borders_" + str(border_width)

    if arguments["-o"] != None:
        output_directory = arguments["-o"].rstrip('/')
    elif input_directory != None:
        output_directory = input_directory + '_borders_' + str(border_width)
    else:
        output_directory = ''

    if output_directory != '' and (not exists(output_directory)):
        makedirs(output_directory)

    if filename != None:
        call_borders(filename, output_directory, half_width)
    else:
        file_list = listdir(input_directory)
        for file in file_list:
            full_path = join(input_directory, file)
            if isdir(full_path):
                continue
            call_borders(full_path, output_directory, half_width)
        # merge BED files for individual chromosomes in one BED file
        filename_list = listdir(output_directory)
        genome_bed_filename = join(output_directory, 'All_borders_' + str(border_width) + '.bed')
        with open(genome_bed_filename, 'w') as dst:
            track_line = 'track name="' + track_name + '" visibility=1 itemRgb="On"'
            dst.write(track_line + '\n')
            for filename in sorted(filename_list):
                with open(filename, 'r') as src:
                    for i, line in enumerate(src):
                        if i == 0:
                            continue
                        dst.write(line)

