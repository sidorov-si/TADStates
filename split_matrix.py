#!/usr/bin/env python
"""
Split a whole-genome contact matrix into per-chromosome contact matrices.

Usage:
  split_matrix.py (-m <contact_matrix> | -d <input_directory>) -r <matrix_resolution> [-o <output_directory>]

Options:
  -h --help               Show this screen.
  --version               Show version.
  -m <contact_matrix>     Whole genome matrix of contact counts.
  -d <input_directory>    Directory with whole genome contact matrices of the same resolution.
  -r <matrix_resolution>  Matrix resolution, e. g., 100Kb or 20 Kb. (The same for all matrices in input directory).
  -o <output_directory>   Output directory.
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
from os.path import isdir
from os.path import join
from os.path import exists
from os.path import isfile
from os import listdir
from os import makedirs

def get_submatrices(matrix_filename, total_output_directory, matrix_resolution):
    with open(matrix_filename, 'r') as src:
        # parse header and create a file for each chromosome
        print "Parse header...",
        chrom_shifts = {}
        header = src.readline().strip(' \n\t')
        header_list = header.split('\t')
        chrom_list = []
        chrom_region_count = 0
        chrom_prev_region_count = 0
        chrom_prev_name = ''
        chrom_curr_name = (header_list[0].split('_'))[0]
        matrix_basename_list = (basename(matrix_filename).split('.'))[:-1]
        matrix_basename = '.'.join(matrix_basename_list)
        matrix_output_directory = join(total_output_directory, matrix_basename + '_per_chrom')
        if not exists(matrix_output_directory):
            makedirs(matrix_output_directory)
        for word in header_list:
            chrom_list.append(word)
            chrom_region_count += 1
            word_list = word.split('_')
            chrom_name = word_list[0]
            if chrom_name != chrom_curr_name:
                chrom_filename = join(matrix_output_directory, chrom_curr_name + \
                                 '_' + matrix_resolution + '.txt')
                with open(chrom_filename, 'w') as dst:
                    del chrom_list[-1]
                    chrom_region_count -= 1
                    chrom_header = '\t'.join(chrom_list)
                    dst.write(chrom_header + '\n')
                    if chrom_prev_name == '':
                        chrom_shifts[chrom_curr_name] = {'current' : chrom_prev_region_count,
                                                         'next' : chrom_region_count}
                    else:
                        curr_chrom_shift = chrom_shifts[chrom_prev_name]['next']
                        next_chrom_shift = curr_chrom_shift + chrom_region_count
                        chrom_shifts[chrom_curr_name] = {'current' : curr_chrom_shift,
                                                         'next' : next_chrom_shift}
                chrom_prev_region_count = chrom_region_count
                chrom_prev_name = chrom_curr_name
                chrom_curr_name = chrom_name
                del chrom_list[:]
                chrom_list.append(word)
                chrom_region_count = 1

        chrom_filename = join(matrix_output_directory, chrom_curr_name + '_' + \
                         matrix_resolution + '.txt')
        with open(chrom_filename, 'w') as dst:
            chrom_header = '\t'.join(chrom_list)
            dst.write(chrom_header + '\n')
            
        if chrom_prev_name == '':
            chrom_shifts[chrom_curr_name] = {'current' : chrom_prev_region_count,
                                             'next' : chrom_region_count}
        else:
            curr_chrom_shift = chrom_shifts[chrom_prev_name]['next']
            next_chrom_shift = curr_chrom_shift + chrom_region_count
            chrom_shifts[chrom_curr_name] = {'current' : curr_chrom_shift,
                                             'next' : next_chrom_shift}
        print "done."

        # cut out a submatrix for each chromosome and 
        # append it to the corresponding file line by line
        print "Cut out submatrices..."
        curr_chrom_name = ''
        line = src.readline().strip('\n')
        while True:
            if not line:
                break
            line_list = line.split('\t')
            curr_region_name = line_list[0]
            chrom_name = (curr_region_name.split('_'))[0]
            if curr_chrom_name != '':
                print "done."
            print "Process chromosome " + chrom_name + "...",
            # +1 because index 0 is for region name            
            first_index = chrom_shifts[chrom_name]['current'] + 1 
            # +1 to include the last element            
            last_index = chrom_shifts[chrom_name]['next'] + 1 
            curr_chrom_name = chrom_name
            chrom_filename = join(matrix_output_directory, chrom_name + '_' + \
                             matrix_resolution + '.txt')
            with open(chrom_filename, 'a') as dst:
                while line and (chrom_name == curr_chrom_name):
                    curr_region_list = line_list[first_index:last_index]
                    curr_region_line = '\t'.join(curr_region_list)
                    curr_region_line = curr_region_name + '\t' + curr_region_line
                    dst.write(curr_region_line + '\n')
                    line = src.readline().strip('\n')
                    if line:
                        line_list = line.split('\t')
                        curr_region_name = line_list[0]
                        chrom_name = (curr_region_name.split('_'))[0]
        print "done."
        print "Done."


if __name__ == '__main__':
    arguments = docopt(__doc__, version='split_matrix 0.3')
    if arguments["-m"] != None:
        matrix_filename = arguments["-m"]
        if not exists(matrix_filename):
            print "Error: Can't find file with contact matrix: no such file '" + \
                  matrix_filename + "'. Exit.\n"
            sys.exit(1)
        if not isfile(matrix_filename):
            print "Error: File with contact matrix must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
        input_directory = None
    else:
        input_directory = arguments["-d"].rstrip('/')
        if not exists(input_directory):
            print "Error: Can't find input directory: no such directory '" + \
                  input_directory + "'. Exit.\n"
            sys.exit(1)
        if not isdir(input_directory):
            print "Error: Input directory must be a directory:). Something else given. Exit.\n"
            sys.exit(1)
        matrix_filename = None

    try:
        matrix_resolution = int(arguments["-r"])
    except ValueError:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)
    if matrix_resolution <= 0:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)

    if arguments["-o"] != None:
        total_output_directory = arguments["-o"].rstrip('/')
    elif input_directory != None:
        total_output_directory = input_directory + '_per_chrom'
    else:
        total_output_directory = ''


    if total_output_directory != '' and (not exists(total_output_directory)):
        makedirs(total_output_directory)

    if matrix_filename != None: # there is only one contact matrix
        print 'Matrix:', matrix_filename
        get_submatrices(matrix_filename, total_output_directory, matrix_resolution)
    else: # there is a directory with contact matrices
        matrix_files = listdir(input_directory)
        for matrix_file in matrix_files:
            matrix_file_full = join(input_directory, matrix_file)
            if isdir(matrix_file_full):
                continue
            print 'Matrix:', matrix_file_full
            get_submatrices(matrix_file_full, total_output_directory, matrix_resolution)
