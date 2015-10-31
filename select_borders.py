#!/usr/bin/env python
"""
Generate BED files with TAD borders that have scores from the particular interval [a,b].
The script is convinient for selecting 'weak' or 'strong' borders.
Borders with high score are called 'strong'. They tend to be conservative and enable 
small amount of interactions to occur over them.

Usage:
  select_borders.py (-t <TADs_filename> -s <scores_filename> [-n <track_name>] | -T <TADs_directory> -S <scores_directory> [-N <track_name_for_all_borders> --all]) -w <border_width> -i <score_interval> [-x <name_suffix> -o <output_directory> -O <output_directory_for_total_track> --scores-as-names --color]

Options:
  -h --help                              Show this screen.
  --version                              Show version.
  -t <TADs_filename>                     Name of a BED file with TAD coordinates in one chromosome.
  -s <scores_filename>                   Name of a TXT file with border scores in one chromosome.
  -n <track_name>                        Name of a track of selected borders. Default: Borders_<border_width>_<scores>.
  -T <TADs_directory>                    Name of a directory containing only BED files with TAD coordinates.
  -S <scores_directory>                  Name of a directory containing only TXT files with border scores.
  -N <track_name_for_all_borders>        Name of a track of all selected borders. Default: All_borders_<border_width>_<scores>.
  --all                                  Consider all files from TADs_directory. Default: consider only chr* files.
  -w <border_width>                      Border width (in bp).
  -i <score_interval>                    Score interval (a:b = [a,b]; :b = [1,b]; a: = [a, 10]; : = [1,10]; a,b in {1,2,...,10}, a <= b).
  -o <output_directory>                  Output directory.
  -O <output_directory_for_total_track>  Output directory for the track with all borders. (It is copied here from output directory).
  -x <name_suffix>                       Suffix for all filenames and tracknames. Default: ''.
  --scores-as-names                      Set scores as border names.
  --color                                Color TAD borders according to their score.
"""


import sys

print

modules = ["docopt", "os", "shutil"]
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
from shutil import copy


score_rgb_dict = {1:'0,0,255', 2:'0,255,255', 3:'0,255,170', 4:'0,255,0', \
                  5:'85,255,0', 6:'170,255,0', 7:'255,255,0', 8:'255,170,0', \
                  9:'255,85,0', 10:'255,0,0'}


def select_borders(tads_filename, scores_filename, start_score, end_score, \
                   chr_output_directory, half_width, scores_as_names, name_suffix, \
                   track_name):
    name, ext = splitext(tads_filename)
    file_basename = basename(name)
    if chr_output_directory == '':
        file_name = name
    else:
        file_name = file_basename
    output_filename = join(chr_output_directory, file_name + '_borders_' + \
                           str(border_width) + '_' + str(start_score) + \
                           '-' + str(end_score) + name_suffix + '.bed')
    with open(tads_filename, 'r') as tads, open(scores_filename, 'r') as scores, \
         open(output_filename, 'w') as dst:
        tads_track_line = tads.readline()
        scores_comment = scores.readline()
        zipped = zip(tads, scores)
        left_coord = ''
        i = 0
        for line_tads, line_scores in zipped:
            i += 1
            line_tads_list = line_tads.rstrip('\n').split()
            line_scores_list = line_scores.rstrip('\n').split()
            chrom_name = line_tads_list[0]
            if i == 1:
                track_line = 'track name="' + track_name + '" visibility=1 itemRgb="On"'
                dst.write(track_line + '\n')
            if left_coord == '':
                left_coord = int(line_tads_list[2]) - half_width
            else:
                right_coord = int(line_tads_list[1]) + half_width
                if scores_as_names:
                    border_name = str(score)
                else:
                    border_name = chrom_name + '.border.' + str(i - 1)
                strand = '.' # Just to fill in the field
                if color_borders:
                    color = score_rgb_dict[score]
                else:
                    color = '0,255,0' # yellow
                border_line = chrom_name + '\t' + str(left_coord) + '\t' + \
                              str(right_coord) + '\t' + border_name + '\t' + \
                              str(score) + '\t' + strand + '\t' + str(left_coord) + '\t' + \
                              str(right_coord) + '\t' + color
                dst.write(border_line + '\n')
                left_coord = int(line_tads_list[2]) - half_width

            score = int(float(line_scores_list[3]))
            if score < start_score or score > end_score:
                left_coord = ''
            

if __name__ == '__main__':
    arguments = docopt(__doc__, version='select_borders 0.6')
    try:
        border_width = int(arguments["-w"])
    except ValueError:
        print "Error: Border width must be an integer greater than 0. Exit.\n"
        sys.exit(1)
    if border_width <= 0:
        print "Error: Border width must be an integer greater than 0. Exit.\n"
        sys.exit(1)
    half_width = border_width / 2

    score_interval = arguments["-i"]
    interval_list = score_interval.split(':')
    if len(interval_list) <= 1:
        print "Error: Score interval should contain ':'. Exit.\n"
        sys.exit(1)

    start_score_str = interval_list[0]
    end_score_str = interval_list[1]

    if start_score_str == '':
        start_score = 1
    else:
        try:
            start_score = int(start_score_str)
        except ValueError:
            print "Error: Start score must be an integer from {1,2,...,10}. Exit.\n"
            sys.exit(1)
        if start_score < 1 or start_score > 10:
            print "Error: Start score must be an integer from {1,2,...,10}. Exit.\n"
            sys.exit(1)

    if end_score_str == '':
        end_score = 10
    else:
        try:
            end_score = int(end_score_str)
        except ValueError:
            print "Error: End score must be an integer from {1,2,...,10}, >= start score. Exit.\n"
            sys.exit(1)
        if end_score < 1 or end_score > 10:
            print "Error: End score must be an integer from {1,2,...,10}, >= start score. Exit.\n"
            sys.exit(1)
        if end_score < start_score:
            print "Error: End score must be an integer from {1,2,...,10}, >= start score. Exit.\n"
            sys.exit(1)

    if arguments["-x"] != None:
        name_suffix = arguments["-x"]
    else:
        name_suffix = ''

    if arguments["-t"] != None:
        tads_filename = arguments["-t"]
        if not exists(tads_filename):
            print "Error: Can't find BED file with TADs: no such file '" + \
                  tads_filename + "'. Exit.\n"
            sys.exit(1)
        if not isfile(tads_filename):
            print "Error: BED file with TADs must be a regular file. Something else given. Exit.\n"
            sys.exit(1)
        scores_filename = arguments["-s"]
        if not exists(scores_filename):
            print "Error: Can't find file with TAD border scores: no such file '" + \
                  scores_filename + "'. Exit.\n"
            sys.exit(1)
        if not isfile(scores_filename):
            print "Error: File with TAD border scores must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
        tads_directory = None
        scores_directory = None
        if arguments["-n"] != None:
            track_name = arguments["-n"]
        else:
            track_name = 'Borders_' + str(border_width) + '_' + str(start_score) + \
                         '-' + str(end_score) + name_suffix
    else:
        tads_filename = None
        scores_filename = None
        tads_directory = arguments["-T"].rstrip('/')
        if not exists(tads_directory):
            print "Error: Can't find directory with TAD BED files: no such directory '" + \
                  tads_directory + "'. Exit.\n"
            sys.exit(1)
        if not isdir(tads_directory):
            print "Error: Directory with TAD BED files must be a directory:). " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
        scores_directory = arguments["-S"].rstrip('/')
        if not exists(scores_directory):
            print "Error: Can't find directory with files containing TAD border scores: " + \
                  "no such directory '" + scores_directory + "'. Exit.\n"
            sys.exit(1)
        if not isdir(scores_directory):
            print "Error: Directory with files containing TAD border scores must be " + \
                  "a directory:). Something else given. Exit.\n"
            sys.exit(1)

        if arguments["-N"]:
            wg_track_name = arguments["-N"]
        else:
            wg_track_name = 'All_borders_' + str(border_width) + '_' + str(start_score) + \
                         '-' + str(end_score) + name_suffix
        track_name = 'Borders_' + str(border_width) + '_' + str(start_score) + \
                         '-' + str(end_score) + name_suffix

        if arguments["--all"]:
            all = True
        else:
            all = False

    if arguments["--scores-as-names"]:
        scores_as_names = True
    else:
        scores_as_names = False

    if arguments["--color"]:
        color_borders = True
    else:
        color_borders = False

    if arguments["-o"] != None:
        output_directory = arguments["-o"].rstrip('/')
    elif tads_directory != None:
        suffix = str(border_width) + '_' + str(start_score) + '-' + str(end_score)
        output_directory = tads_directory + '_borders_' + suffix
    else:
        output_directory = ''

    if output_directory != '' and (not exists(output_directory)):
        makedirs(output_directory)

    chr_output_directory = join(output_directory, 'per_chr')
    if not exists(chr_output_directory):
        makedirs(chr_output_directory)

    output_directory_total = None
    if arguments["-O"] != None:
        output_directory_total = arguments["-O"]
        if not exists(output_directory_total):
            makedirs(output_directory_total)

    if tads_filename != None:
        select_borders(tads_filename, scores_filename, start_score, end_score, \
                       chr_output_directory, half_width, scores_as_names, name_suffix, \
                       track_name)
    else:
        files_list = sorted(listdir(tads_directory))
        if all:
            tads_list = files_list[:]
        else:
            tads_list = [file for file in files_list if file.startswith('chr')]
        scores_list = sorted(listdir(scores_directory))
        zipped = zip(tads_list, scores_list)
        for tads, scores in zipped:
            tads_full_path = join(tads_directory, tads)
            scores_full_path = join(scores_directory, scores)
            select_borders(tads_full_path, scores_full_path, start_score, end_score, \
                           chr_output_directory, half_width, scores_as_names, name_suffix, \
                           track_name)
        # merge BED files for individual chromosomes in one BED file
        filename_list = listdir(chr_output_directory)
        files_list = [join(chr_output_directory, f) for f in filename_list]
        genome_bed_filename = join(output_directory, 'All_borders_' + str(border_width) + \
                                   '_' + str(start_score) + '-' + str(end_score) + \
                                   name_suffix + '.bed')
        with open(genome_bed_filename, 'w') as dst:
            track_line = 'track name="' + wg_track_name + '" visibility=1 itemRgb="On"'
            dst.write(track_line + '\n')
            for filename in sorted(files_list):
                with open(filename, 'r') as src:
                    for i, line in enumerate(src):
                        if i == 0:
                            continue
                        dst.write(line)

        if output_directory_total != None:
            copy(genome_bed_filename, output_directory_total)

