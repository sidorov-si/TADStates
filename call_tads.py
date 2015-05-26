#!/usr/bin/env python
"""
Generate BED files with TAD coordinates from contact matrices (or one contact matrix).
For TADs search TADbit library is used.

Usage:
  call_tads.py (-m <contact_matrix> [-c <chromosome_name>] | -d <input_directory>) -r <matrix_resolution> [-t <thread_number> -n <track_name_for_all_TADs> -o <output_directory>] 

Options:
  -h --help                     Show this screen.
  --version                     Show version.
  -m <contact_matrix>           Matrix of contact counts (.mat-file).
  -d <input_directory>          Directory with contact matrices (.mat-files).
  -c <chromosome_name>          Name of the chromosome. Determined from matrix file by default.
  -r <matrix_resolution>        Matrix resolution.
  -t <thread_number>            Number of threads for TADbit. Default: 4.
  -n <track_name_for_all_TADs>  A name for a track with all TADs. Default: All_TADs.
  -o <output_directory>         Output directory name. Default: directory that contains this script.
"""

import sys

print

modules = ["docopt", "pytadbit", "os"]
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
from pytadbit import Chromosome
from os.path import basename
from os.path import splitext
from os.path import join
from os.path import exists
from os.path import isfile
from os.path import isdir
from os import makedirs
from os import listdir
from re import search


def call_tads(matrix_filename, chrom_name):
    print "Processing chromosome " + chrom_name + '...'
    print "Contact matrix:" + matrix_filename
    chrom_number = search(r'\d+|X|Y', chrom_name).group(0)
    if len(chrom_number) == 1 and chrom_number != 'X' and chrom_number != 'Y':
        chrom_number = '0' + chrom_number
        chrom_id = 'chr' + chrom_number
    else:
        chrom_id = chrom_name
        
    matrix_basename = splitext(basename(matrix_filename))[0]
    output_txt_filename = join(txt_directory, chrom_id + '_TADs.txt')
    output_bed_filename = join(bed_directory, chrom_id + '_TADs.bed')
    filename_list.append(output_bed_filename)

    # Call TADs and write their borders in TADbit text format and in BED format
    chrom = Chromosome(name=chrom_name)
    chrom.add_experiment(chrom_name, hic_data=matrix_filename, resolution=matrix_resolution)
    chrom.find_tad(chrom_name, n_cpus=thread_number)
    chrom.experiments[chrom_name].write_tad_borders(savedata=output_txt_filename)
    with open(output_txt_filename, 'r') as src, open(output_bed_filename, 'w') as dst:
        track_line = 'track name="' + chrom_name + '_TADs" visibility=1 itemRgb="On"'
        dst.write(track_line + '\n')
        for i, line in enumerate(src):
            if line.split()[0] == '#':
                continue
            line_list = line.split()
            tad_name =  chrom_name + '.' + 'TAD' + '.' + str(i)
            # Coordinates in BED format are 0-based, 
            # and a region is presented by [x,y) interval.
            start_pos = (int(line_list[1]) - 1) * matrix_resolution
            end_pos = int(line_list[2]) * matrix_resolution
            score = 0 # Just to fill in the field
            strand = '.' # Just to fill in the field
            if i%2:
                color = '0,0,255' # blue
            else:
                color = '255,0,0' # red
            bed_line = chrom_name + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + \
                        tad_name + '\t' + str(score) + '\t' + strand + '\t' + \
                        str(start_pos) + '\t' + str(end_pos) + '\t' + color
            dst.write(bed_line + '\n')
    chrom_filename = join(tdb_directory, chrom_id + '.tdb')
    chrom.save_chromosome(chrom_filename, force=True)


def get_chrom_name(matrix_filename):
    with open(matrix_filename, 'r') as src:
        line_list = src.readline().strip().split('\t ')
        chrom_name = (line_list[0].split('_'))[0]
    return chrom_name


if __name__ == '__main__':
    arguments = docopt(__doc__, version='call_tads 0.5')
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
        if arguments["-c"] != None:
            chrom_name = arguments["-c"]
        else:
            chrom_name = None
        input_directory = None
    else:
        matrix_filename = None
        chrom_name = None
        input_directory = arguments["-d"].rstrip('/')
        if not exists(input_directory):
            print "Error: Can't find input directory: no such directory '" + \
                  input_directory + "'. Exit.\n"
            sys.exit(1)
        if not isdir(input_directory):
            print "Error: Input directory must be a directory:). Something else given. Exit.\n"
            sys.exit(1)

    try:
        matrix_resolution = int(arguments["-r"])
    except ValueError:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)
    if matrix_resolution <= 0:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)

    if arguments["-t"] != None:
         try:
             thread_number = int(arguments["-t"])
         except ValueError:
             print "Error: Thread number must be an integer greater than 0. Exit.\n"
             sys.exit(1)
         if thread_number <= 0:
             print "Error: Thread number must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
        thread_number = 4

    if arguments["-n"] != None:
        track_name = arguments["-n"]
    else:
        track_name = "All_TADs"

    if arguments["-o"] != None:
        output_directory = arguments["-o"].rstrip('/')
    else:
        if input_directory != None:
            output_directory = input_directory + '_TADs'
        else:
            output_directory = ''

    filename_list = []
    if output_directory != '':
        if not exists(output_directory):
            makedirs(output_directory)
    bed_directory = join(output_directory, 'BED')
    txt_directory = join(output_directory, 'TXT')
    tdb_directory = join(output_directory, 'TDB')
    if not exists(bed_directory):
        makedirs(bed_directory)
    if not exists(txt_directory):
        makedirs(txt_directory)
    if not exists(tdb_directory):
        makedirs(tdb_directory)

    if matrix_filename != None: # there is only one contact matrix
        if chrom_name == None:
            chrom_name = get_chrom_name(matrix_filename)
        call_tads(matrix_filename, chrom_name)
    else: # there is a directory with matrices
        files = listdir(input_directory)
        for file in files:
            if splitext(basename(file))[1] == '.mat':
                matrix_file_full = join(input_directory, file)
                chrom_name = get_chrom_name(matrix_file_full)
                call_tads(matrix_file_full, chrom_name)
        # merge BED files for individual chromosomes in one BED file
        genome_bed_filename = join(bed_directory, 'All_TADs.bed')
        with open(genome_bed_filename, 'w') as dst:
            track_line = 'track name="' + track_name + '" visibility=1 itemRgb="On"'
            dst.write(track_line + '\n')
            for filename in sorted(filename_list):
                with open(filename, 'r') as src:
                    for i, line in enumerate(src):
                        if i == 0:
                            continue
                        dst.write(line)
