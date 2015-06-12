#!/usr/bin/env python
"""
Calculate enrichment of regions with states using ChromHMM.

Usage:
  calc_enr.py (-r <regions_file> | -R <directory_with_regions_files>) (-s <segmentation_directory> | -S <directory_with_segmentation_directories>) -c <ChromHMM_directory> -o <output_directory>

Options:
  -h --help                                     Show this screen.
  --version                                     Show version.
  -r <regions_file>                             BED file with regions to calc enrichments for. 
  -R <directory_with_regions_files>             Directory with BED files containing regions to calc enrichments for.          
  -s <segmentation_directory>                   Directory with segmentation produced by ChromHMM.
  -S <directory_with_segmentation_directories>  Directory with directories containing segmentations.
  -c <ChromHMM_directory>                       ChromHMM directory.
  -o <output_directory>                         Output directory name.
"""

import sys

print

modules = ["docopt", "os", "subprocess"]
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
from os.path import isdir
from os.path import isfile
from os import makedirs
from os import listdir
from subprocess import call
from sys import stdout


def calc_enr(regions, segm_files, png_directory, svg_directory, txt_directory, chromhmm_directory):
    regions_part = splitext(basename(regions))[0]
    for segm_file in segm_files:
        print
        print 'Calc enrichment for', basename(regions), 'and', basename(segm_file), '...'
        stdout.flush()
        segm_part = splitext(basename(segm_file))[0]
        prefix = regions_part + '_' + segm_part
        command_line_list = ['java', '-mx1600M', '-jar', join(chromhmm_directory, 'ChromHMM.jar'), \
                             'OverlapEnrichment', segm_file, regions, prefix]
        code = call(command_line_list)
        if code != 0:
            print 'Something went wrong!'
        else:
            print 'Done.'
            call(['mv', prefix + '.png', png_directory])
            call(['mv', prefix + '.svg', svg_directory])
            call(['mv', prefix + '.txt', txt_directory])


if __name__ == '__main__':
    arguments = docopt(__doc__, version='calc_enr 0.2')
    if arguments["-r"] != None:
        regions = arguments["-r"]
        if not exists(regions):
            print "Error: Can't find BED file with regions: no such file '" + \
                  regions + "'. Exit.\n"
            sys.exit(1)
        if not isfile(regions):
            print "Error: BED file with regions must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
    else:
        regions = arguments["-R"].rstrip('/')
        if not exists(regions):
            print "Error: Can't find directory with region BED files: no such directory '" + \
                  regions + "'. Exit.\n"
            sys.exit(1)
        if not isdir(regions):
            print "Error: Directory with region BED files must be a directory:). " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
            
    if arguments["-s"] != None:
        segm_dir = arguments["-s"].rstrip('/')
        if not exists(segm_dir):
            print "Error: Can't find directory with segmentation: no such directory '" + \
                  segm_dir + "'. Exit.\n"
            sys.exit(1)
        if not isdir(segm_dir):
            print "Error: Directory with segmentation must be a directory:). " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
        segm_directory = None
    else:
        segm_dir = None
        segm_directory = arguments["-S"].rstrip('/')
        if not exists(segm_directory):
            print "Error: Can't find directory with directories containing segmentations: " + \
                  "no such directory '" + segm_directory + "'. Exit.\n"
            sys.exit(1)
        if not isdir(segm_directory):
            print "Error: Directory with directories containing segmentations must " + \
                  "be a directory:). Something else given. Exit.\n"
            sys.exit(1)

    chromhmm_directory = arguments["-c"].rstrip('/')
    if not exists(chromhmm_directory):
        print "Error: Can't find ChromHMM directory: no such directory '" + \
              chromhmm_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(chromhmm_directory):
        print "Error: ChromHMM directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    output_directory = arguments["-o"].rstrip('/')

    if not exists(output_directory):
        makedirs(output_directory)
    png_directory = join(output_directory, 'PNG')
    svg_directory = join(output_directory, 'SVG')
    txt_directory = join(output_directory, 'TXT')
    if not exists(png_directory):
        makedirs(png_directory)
    if not exists(svg_directory):
        makedirs(svg_directory)
    if not exists(txt_directory):
        makedirs(txt_directory)

    if segm_dir != None: # there is only one segmentation
        segm_dirs = [segm_dir]
    else:
        segm_dirnames = sorted(listdir(segm_directory))
        segm_dirs = [join(segm_directory, d) for d in segm_dirnames]

    segm_files = []
    for dir in segm_dirs:
        filenames_list = listdir(dir)
        segm_filenames = [f for f in filenames_list if 'segments' in f]
        segm_files.extend([join(dir, f) for f in segm_filenames])

    calc_enr(regions, segm_files, png_directory, svg_directory, txt_directory, chromhmm_directory)
