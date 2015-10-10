#!/usr/bin/env python
"""
Calculate cross-window score (CWS) for each border between Hi-C windows 
(the width of a Hi-C window is equal to a contact matrix resolution, and a genome is split into back-to-back windows).
CWS of a border is a number of contacts that cross the border. If vicinity_size is set then during CWS calculation only those contacts are considered that connect regions (vicinity_size / 2) bp upstream and downstream the border. Othewise, CWS is calculated across the whole genome.
For each chromosome a BED file is created with all border CWS values.
All BED files are also concatinated to the whole genome BED file.
Finally, CWS is plotted for each chromosome, and the graph is stored into correspondent PNG file.
BED file with TADs can also be set as input, and TAD border coordinates will be brought from there and showed on plots.
With -R option you can set a specific region within the chromosome. Only this region will be plotted. Coordinates must be set in bp and be multiples of matrix resolution.

Usage:
  calc_cws.py (-m <contact_matrix> [-c <chromosome_name> -n <track_name> -R <chromosome_region> -T <BED_file_with_TADs>] | -d <input_directory> [-N <track_name_for_whole_genome_BED> -O <output_whole_genome_BED_file>]) -r <matrix_resolution> [-o <output_directory> -e <vicinity_size>]

Options:
  -h --help                             Show this screen.
  --version                             Show version.
  -m <contact_matrix>                   File with contact matrix.
  -d <input_directory>                  Directory with contact matrices (.mat-files).
  -c <chromosome_name>                  Name of the chromosome. Determined from matrix file by default.
  -r <matrix_resolution>                Matrix resolution.
  -e <vicinity_size>                    Vicinity size, bp. Default: CWS is calculated across the whole genome.
  -n <track_name>                       A name for the track for one chromosome. Default: 'chromosome_name'_CWS
  -R <chromosome_region>                A region on a chromosome (a:b = [a,b]; :b = [matrix_resolution,b]; a: = [a, chromosome_end]; : = the whole chromosome; a, b are set in bp, a < b, a and b are multiples of matrix_resolution).
  -T <BED_file_with_TADs>               BED file with TADs. TAD border coordinates are brought from there.
  -N <track_name_for_whole_genome_BED>  A name for the track with all TADs. Default: All_TADs.
  -o <output_directory>                 Output directory name. Default: directory that contains this script.
  -O <output_whole_genome_BED_file>     Output whole genome BED file. Is also put in the output directory. Default: generated out of the contact matrix filename.
"""

import sys

print

modules = ["docopt", "os", "numpy", "matplotlib", "re"]
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
from sys import stdout
from os import makedirs
from os import listdir
from re import search
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def calc_cws(contact_matrix_filename, chrom_name):
    print
    print 'Contact matrix file:', contact_matrix_filename
    print 'Matrix resolution:  ', matrix_resolution
    print 'BED track name:     ', track_name
    if vicinity_size != -1:
        print 'Vicinity size:      ', vicinity_size
    else:
        print 'Vicinity size:       Calculate CWS across the whole genome.'
    stdout.flush()
    print 'Region to plot:     ',
    if start_coord == None and end_coord == None:
        print 'The whole chromosome.'
    else:
        print 'From', str(start_coord), 'to', str(end_coord), 'bp.'
    chrom_number = search(r'\d+|X|Y', chrom_name).group(0)
    if len(chrom_number) == 1 and chrom_number != 'X' and chrom_number != 'Y':
        chrom_number = '0' + chrom_number
        chrom_id = 'chr' + chrom_number
    else:
        chrom_id = chrom_name
    output_bed_filename = join(bed_directory, chrom_id + '_CWS.bed')
    filename_list.append(output_bed_filename)
    output_png_filename = join(png_directory, chrom_id + '_CWS.png')
    print 'Output BED file:    ', output_bed_filename
    print 'Output PNG file:    ', output_png_filename
    stdout.flush()

    # Calculate CWS for all borders between windows
    with open(contact_matrix_filename, 'r') as infile:
        print
        print 'Calculate CWS for chromosome', chrom_name, '...',
        stdout.flush()
        header = infile.readline()
        first_line = (infile.readline().rstrip('\n').split())[1:]
        matrix_width = len(first_line)
        contact_matrix = numpy.empty(shape=(matrix_width, matrix_width))
        line_number = 0
        while True:
            raw_line = infile.readline()
            if not raw_line:
                break; # EOF is reached
            line = (raw_line.rstrip('\n').split())[1:]
            if first_line != None:
                contact_matrix[0] = first_line
                first_line = None
                line_number += 1
            contact_matrix[line_number] = line
            line_number += 1
        
        # calculate global CWS (over the whole genome)    
        cws = numpy.empty(shape=(matrix_width, matrix_width))
        n = matrix_width - 1 # the number of the last line in the contact matrix
        # for local CWS consider only regions k windows upstream and downstream
        k = vicinity_size / (2 * matrix_resolution) 
        # initial values
        cws[n, 0] = contact_matrix[n, 0]
        for j in range(1, n):
            cws[n, j] = cws[n, j - 1] + contact_matrix[n, j]
        for i in range(n - 1, 0, -1):
            cws[i, 0] = cws[i + 1, 0] + contact_matrix[i, 0]
        # dynamic programming over the cws matrix;
        # diagonals are calculated one by one
        for diagonal_index in range(1, n - 1):
            for i, j in zip(range(n - diagonal_index, n), range(1, diagonal_index + 1)):
                cws[i, j] = cws[i, j - 1] + cws[i + 1, j] - cws[i + 1, j - 1] + contact_matrix[i, j]
        
        if vicinity_size != -1:
            # calculate local CWS 
            # (number of contacts in the 2k-vicinity of a border that cross the border)
            for i, j in zip(range(1, n + 1), range(0, n)):
                if j - k + 1 < 0:
                    left_rect_sum = 0
                    left_bottom_rect_sum = 0
                else:
                    left_rect_sum = cws[i, j - k + 1]
                if i + k - 1 > n:
                    bottom_rect_sum = 0
                    left_bottom_rect_sum = 0
                else:
                    bottom_rect_sum = cws[i + k - 1, j]
                if not j - k + 1 < 0 and not i + k - 1 > n:
                    left_bottom_rect_sum = cws[i + k - 1, j - k + 1]
                cws[i, j] = cws[i, j] - left_rect_sum - bottom_rect_sum + left_bottom_rect_sum

        # CWS for all borders between windows
        result = [cws[i, j] for i, j in zip(range(1, n + 1), range(0, n))]
        print 'Finish.'
        stdout.flush()
        
        # Generate BED file with CWS for all borders between windows
        print 'Generate BED file with CWS for chromosome', chrom_name, '...',
        stdout.flush()
        with open(output_bed_filename, 'w') as dst:
            track_line = 'track name="' + track_name + '" visibility=1 itemRgb="On"'
            dst.write(track_line + '\n')
            for border_number, cws_value in enumerate(result):
                border_name = str(cws_value)
                # Coordinates in BED format are 0-based, 
                # and a region is presented by [x,y) interval.
                start_pos = (border_number + 1) * matrix_resolution - 1
                end_pos = (border_number + 1) * matrix_resolution
                score = cws_value
                strand = '.' # Just to fill in the field
                color = '0,0,0' # black
                bed_line = chrom_name + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + \
                           border_name + '\t' + str(score) + '\t' + strand + '\t' + \
                           str(start_pos) + '\t' + str(end_pos) + '\t' + color
                dst.write(bed_line + '\n')
        print 'Finish.'
        stdout.flush()

        # Generate PNG file with CWS plot
        print 'Generate PNG file with CWS plot for chromosome', chrom_name, '...',
        stdout.flush()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if tads_filename == None:
            plot_header = 'CWS for ' + chrom_name
        else:
            plot_header = 'CWS for ' + chrom_name + ' with marked TAD borders'
        if start_coord == None and end_coord == None:
            borders_count = len(result)
            ind = numpy.arange(matrix_resolution, (borders_count + 1) * matrix_resolution, matrix_resolution)
            ax.plot(ind, result, '.-', color = 'blue')
        else:
            start_number = start_coord / matrix_resolution - 1
            if end_coord != chr_len:
               end_number = end_coord / matrix_resolution - 1
            else:
               end_number = end_coord / matrix_resolution - 2
            region_result = result[start_number : end_number + 1]
            additional_value = matrix_resolution if end_coord != chr_len else 0
            ind = numpy.arange(start_coord, end_coord + additional_value, matrix_resolution)
            ax.plot(ind, region_result, '.-', color = 'blue')
            plot_header += '. Region: ' + str(start_coord) + '-' + str(end_coord) + ' bp'

        if tads_filename != None:
        # also plot TAD borders
            with open(tads_filename, 'r') as tads:
                tad_border_coords = []
                for i, line in enumerate(tads):
                    if i == 0:
                        continue # leave out the header
                    line_fields = line.rstrip('\n').split('\t')
                    border_coord = int(line_fields[2])
                    tad_border_coords.append(border_coord)
                tad_border_coords.pop()
                tad_border_numbers = [coord / matrix_resolution - 1 for coord in tad_border_coords]
                tad_border_cws = [result[border_number] for border_number in tad_border_numbers]
                ax.plot(tad_border_coords, tad_border_cws, '.', color = 'red')

        if start_coord == None:
            max_cws = max(result)
            last_border_number = len(result) - 1
            ax.set_xlim(0, (last_border_number + 2) * matrix_resolution)
            ax.set_ylim(0, max_cws * 1.05)
        else:
            ax.set_xlim(start_coord - matrix_resolution, end_coord + matrix_resolution)
            region_max_cws = max(region_result)
            ax.set_ylim(0, region_max_cws * 1.05)
        ax.set_xlabel('Chromosome coordinates, bp')
        ax.set_ylabel('CWS')
        ax.set_title(plot_header)
        plt.savefig(output_png_filename)
        print 'Finish.'
        stdout.flush()


def get_chrom_name(matrix_filename):
    with open(matrix_filename, 'r') as src:
        line_list = src.readline().strip().split('\t ')
        chrom_name = (line_list[0].split('_'))[0]
    return chrom_name


if __name__ == '__main__':
    arguments = docopt(__doc__, version='calc_cws 0.4')

    try:
        matrix_resolution = int(arguments["-r"])
    except ValueError:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)
    if matrix_resolution <= 0:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)

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
        if arguments["-n"] != None:
            track_name = arguments["-n"]
        else:
            track_name = None
       
        if arguments["-R"] != None:
            with open(matrix_filename, 'r') as matrix_file:
                header = matrix_file.readline().rstrip('\n').split('\t')
                chr_len = len(header) * matrix_resolution # approximate chromosome length
            chr_region = arguments["-R"]
            region_list = chr_region.split(':')
            if len(region_list) <= 1:
                print "Error: Chromosome region value should contain ':'. Exit.\n"
                sys.exit(1)
            start_coord_str = region_list[0]
            end_coord_str = region_list[1]
            if start_coord_str == '':
                start_coord = matrix_resolution
            else:
                try:
                    start_coord = int(start_coord_str)
                except ValueError:
                    print "Error: Start coordinate of the chromosome region must be an integer. Exit.\n"
                    sys.exit(1)
                if start_coord < matrix_resolution or start_coord >= chr_len:
                    print "Error: Start coordinate must be >= matrix_resolution and < " + str(chr_len) + ". Exit.\n"
                    sys.exit(1)
                if start_coord % matrix_resolution != 0:
                    print "Error: Start coordinate must be a multiple of matrix resolution ( " + \
                          str(matrix_resolution) + "). Exit.\n"
                    sys.exit(1)
            
            if end_coord_str == '':
                end_coord = chr_len
            else:
                try:
                    end_coord = int(end_coord_str)
                except ValueError:
                    print "Error: End coordinate of the chromosome region must be an integer. Exit.\n"
                    sys.exit(1)
                if end_coord <= matrix_resolution or end_coord > chr_len:
                    print "Error: End coordinate must be > matrix_resolution and <= " + str(chr_len) + ". Exit.\n"
                    sys.exit(1)
                if end_coord % matrix_resolution != 0:
                    print "Error: End coordinate must be a multiple of matrix resolution ( " + \
                          str(matrix_resolution) + "). Exit.\n"
                    sys.exit(1)
                if end_coord <= start_coord:
                    print "Error: End coordinate must be greater than start coordinate. Exit.\n"
                    sys.exit(1)
        else:
            start_coord = None
            end_coord = None

        if arguments["-T"] != None:
            tads_filename = arguments["-T"]
            if not exists(tads_filename):
                print "Error: Can't find BED file with TADs: no such file '" + \
                      tads_filename + "'. Exit.\n"
                sys.exit(1)
            if not isfile(tads_filename):
                print "Error: BED file with TADs must be a regular file. " + \
                      "Something else given. Exit.\n"
                sys.exit(1)        
        input_directory = None
        output_wg_bed_filename = None
        all_track_name = None
    else:
        matrix_filename = None
        chrom_name = None
        tads_filename = None
        start_coord = None
        end_coord = None
        input_directory = arguments["-d"].rstrip('/')
        if not exists(input_directory):
            print "Error: Can't find input directory: no such directory '" + \
                  input_directory + "'. Exit.\n"
            sys.exit(1)
        if not isdir(input_directory):
            print "Error: Input directory must be a directory:). Something else given. Exit.\n"
            sys.exit(1)
        all_track_name = None
        if arguments["-N"] != None:
            all_track_name = arguments["-N"]
        else:
            all_track_name = "All_border_CWS"
        if arguments["-O"] != None:
            output_wg_bed_filename = arguments["-O"]
        else:
            output_wg_bed_filename = 'All_borders_CWS.bed'

    if arguments["-o"] != None:
        output_directory = arguments["-o"].rstrip('/')
    else:
        if input_directory != None:
            output_directory = input_directory + '_CWS'
        else:
            output_directory = ''
    if arguments["-e"] != None:
        try:
            vicinity_size = int(arguments["-e"])
        except ValueError:
            print "Error: Vicinity size must be an integer greater than 0. Exit.\n"
            sys.exit(1)
        if vicinity_size <= 0:
            print "Error: Vicinity size must be an integer greater than 0. Exit.\n"
            sys.exit(1)
    else:
        vicinity_size = -1

    filename_list = []
    if output_directory != '':
        if not exists(output_directory):
            makedirs(output_directory)
    bed_directory = join(output_directory, 'BED_CWS')
    png_directory = join(output_directory, 'PNG_CWS')
    if not exists(bed_directory):
        makedirs(bed_directory)
    if not exists(png_directory):
        makedirs(png_directory)

    if input_directory != None:
        print 'Input directory:             ', input_directory
        stdout.flush()
    if output_directory != None:
        print 'Output directory:            ', output_directory
        stdout.flush()
    if output_wg_bed_filename != None:
        print 'Output whole genome BED file:', output_wg_bed_filename
        stdout.flush()
    if all_track_name != None:
        print 'Whole genome BED track name: ', all_track_name
        stdout.flush()

    if matrix_filename != None: # there is only one contact matrix
        if chrom_name == None:
            chrom_name = get_chrom_name(matrix_filename)
        if track_name == None:
            track_name = chrom_name + '_CWS'
        calc_cws(matrix_filename, chrom_name)
    else: # there is a directory with matrices
        print
        print 'Calculate CWS for all chromosomes in the input directory...'
        stdout.flush()
        files = listdir(input_directory)
        for file in files:
            if splitext(basename(file))[1] == '.mat':
                matrix_file_full = join(input_directory, file)
                chrom_name = get_chrom_name(matrix_file_full)
                track_name = chrom_name + '_CWS'
                calc_cws(matrix_file_full, chrom_name)
        print 
        print 'All chromosomes are processed.'
        stdout.flush()
        # merge BED files for individual chromosomes in one BED file
        print 'Generate whole genome BED file with CWS...',
        stdout.flush()
        genome_bed_filename = join(bed_directory, output_wg_bed_filename)
        with open(output_wg_bed_filename, 'w') as dst:
            track_line = 'track name="' + all_track_name + '" visibility=1 itemRgb="On"'
            dst.write(track_line + '\n')
            for filename in sorted(filename_list):
                with open(filename, 'r') as src:
                    for i, line in enumerate(src):
                        if i == 0:
                            continue
                        dst.write(line)
        print 'Finish.'
        stdout.flush()
    print 'Processing is finished.'
    stdout.flush()

