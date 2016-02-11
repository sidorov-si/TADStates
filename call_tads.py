#!/usr/bin/env python
"""
Generate BED files with TAD coordinates from contact matrices 
(or one contact matrix). For TADs search TADbit library is used.
TADs can be called from one Hi-C experiment or several Hi-C 
experiments for the same chromosomes. In the latter case all
experiments must have the same resolution, names of all contact 
matrices for the same chromosome must begin with the same 
chromosome ID (e. g., chr1 or chr01). If several input directories
are set, then all of them must contain the same count of contact 
matrix files for the same set of chromosomes.

All files contained in the input directory or directories are 
considered to be contact matrix files.

Usage:
  call_tads.py -m <contact_matrix_list> -r <matrix_resolution> [-c <chromosome_name> -t <thread_number> -o <output_directory>] 
  call_tads.py -d <input_directories_list> -r <matrix_resolution> [-n <track_name_for_all_TADs> -t <thread_number> -o <output_directory>] 

Options:
  -h --help                     Show this screen.
  --version                     Show version.
  -m <contact_matrix>           One contact matrix (.mat-file) or a list of contact matrices for one chromosome. Several contact matrix filenames must be separated with comma.
  -d <input_directory>          Directory with contact matrices (.mat-files) of a list of directories with contact matrices. Several directory names must be separated with comma.
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


def call_tads(matrix_filenames, chrom_name):
    print
    print "Call TADs for chromosome " + chrom_name + '...'
    print "Contact matrices: "
    for matrix_filename in matrix_filenames:
        print matrix_filename
    chrom_number = search(r'\d+|X|Y', chrom_name).group(0)
    if len(chrom_number) == 1 and chrom_number != 'X' and chrom_number != 'Y':
        chrom_number = '0' + chrom_number
        chrom_id = 'chr' + chrom_number
    else:
        chrom_id = chrom_name
        
    output_txt_filename = join(txt_directory, chrom_id + '_TADs.txt')
    print 'Output TXT file:', output_txt_filename
    output_bed_filename = join(bed_directory, chrom_id + '_TADs.bed')
    print 'Output BED file:', output_bed_filename
    filename_list.append(output_bed_filename)

    """tads_2D_filename = join(png_directory, chrom_id + '_TADs_2D.png')
    print 'Output 2D TAD plot file:', tads_2D_filename
    tads_1D_filename = join(png_directory, chrom_id + '_TADs_1D.png')
    print 'Output 1D TAD plot file:', tads_1D_filename"""

    # Call TADs and write their borders in TADbit text format and in BED format
    chrom = Chromosome(name=chrom_name)
    if len(matrix_filenames) > 1: # several matrices for one chromosome
        combined_experiment_name = 'batch'
        experiment_names = []
        for matrix_index, matrix_filename in enumerate(matrix_filenames):
            experiment_name = splitext(basename(matrix_filename))[0] + '_' + str(matrix_index)
            experiment_names.append(experiment_name)
            combined_experiment_name += '_' + experiment_name 
            chrom.add_experiment(experiment_name, hic_data=matrix_filename, \
                                 resolution=matrix_resolution)
        chrom.find_tad(experiment_names, batch_mode = True, n_cpus=thread_number)
        chrom.experiments[combined_experiment_name].write_tad_borders(savedata=output_txt_filename)
        #chrom.visualize(combined_experiment_name, paint_tads=True, savefig=tads_2D_filename)
        #chrom.tad_density_plot(combined_experiment_name, savefig=tads_1D_filename)
    else: # only one matrix for one chromosome
        matrix_filename = matrix_filenames[0]
        experiment_name = splitext(basename(matrix_filename))[0]
        chrom.add_experiment(experiment_name, hic_data=matrix_filename, \
                             resolution=matrix_resolution)
        chrom.find_tad(experiment_name, n_cpus=thread_number)
        chrom.experiments[experiment_name].write_tad_borders(savedata=output_txt_filename)
        #chrom.visualize(experiment_name, paint_tads=True, savefig=tads_2D_filename)
        #chrom.tad_density_plot(experiment_name, savefig=tads_1D_filename)

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
    print 'Finish.'


def get_chrom_name(matrix_filename):
    with open(matrix_filename, 'r') as src:
        line_list = src.readline().strip().split('\t ')
        chrom_name = (line_list[0].split('_'))[0]
    return chrom_name


if __name__ == '__main__':
    arguments = docopt(__doc__, version='call_tads 0.7')
    if arguments["-m"] != None:
        matrix_filenames = arguments["-m"].split(",")
        for matrix_filename in matrix_filenames:
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
        input_directories = None
    elif arguments["-d"] != None:
        matrix_filenames = None
        chrom_name = None
        input_directories = arguments["-d"].split(",")
        for input_directory in input_directories:
            input_directory = input_directory.rstrip('/')
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
    #png_directory = join(output_directory, 'PNG')
    if not exists(bed_directory):
        makedirs(bed_directory)
    if not exists(txt_directory):
        makedirs(txt_directory)
    if not exists(tdb_directory):
        makedirs(tdb_directory)
    #if not exists(png_directory):
    #    makedirs(png_directory)

    if matrix_filenames != None: # there is only one chromosome to process
        if chrom_name == None:
            chrom_name = get_chrom_name(matrix_filenames[0])
        call_tads(matrix_filenames, chrom_name)
    else: # there is a directory with matrices or a list of directories with matrices
        # prepare a list of contact matrix filenames for each chromosome
        # all these lists are saved in a list
        prev_count = -1
        prev_directory = ''
        full_list_of_files = []
        for directory in input_directories:
            files = sorted(listdir(directory))
            if prev_count != -1 and len(files) != prev_count:
                print "Error: The number of contact matrices in directory " + directory + \
                      " is not equal to the number of contact matricws in directory" + \
                      prev_directory + ". All directories must contain the same number of " + \
                      "contact matrices (files)."
                sys.exit(1)
            else:
                prev_count = len(files)
                prev_directory = directory
                full_list_of_files.append(files)
        full_list_per_chr = [list(chr_tuple) for chr_tuple in zip(*full_list_of_files)]
        for matrix_filenames in full_list_per_chr: 
            # call TADs for each chromosome
            matrix_filenames_full = [join(directory, file) for directory, file in \
                                     zip(input_directories, matrix_filenames)]
            chrom_name = get_chrom_name(matrix_filenames_full[0])
            call_tads(matrix_filenames_full, chrom_name)
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

