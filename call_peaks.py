#!/usr/bin/env python
"""
Call peaks for each mark with SICER and MACS v1.4.
Input BAM files are REPLACED with sorted BAM files.
BED files with aligned reads are added to BAM files.

Usage:
  call_peaks.py -d <input_directory> -s <SICER_directory> -g <IGVTools_directory> -o <output_directory> [-i <reference_genome_tag> -r <redundancy_threshold_for_SICER> -w <window_size_for_SICER> -f <fragment_size_for_SICER> -l <effective_genome_fraction_for_SICER> -e <e-value_for_SICER>]

Options:
  -h --help                                 Show this screen.
  --version                                 Show version.
  -d <input_directory>                      Directory with BAM files for all marks to be processed.
  -s <SICER_directory>                      Directory with SICER scripts.
  -g <IGVTools_directory>                   Directory with IGVTools.
  -o <output_directory>                     Output directory name.
  -i <reference_genome_tag>                 Reference genome tag. Only 'hg19' is currently supported. Default: hg19.
  -r <redundancy_threshold_for_SICER>       Value of redundancy threshold option for SICER. Default: 1.
  -w <window_size_for_SICER>                Value of window size option for SICER. Default: 200.
  -f <fragment_size_for_SICER>              Value of fragment size option for SICER. Default: 250.
  -l <effective_genome_fraction_for_SICER>  Value of effective genome fraction option for SICER. Default: 0.77.
  -e <e-value_for_SICER>                    E-value for SICER. Default: 100.
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
from os.path import isfile
from os.path import isdir
from os import makedirs
from os import listdir
from os import system
from os import devnull
from os import getcwd
from sys import stdout
from subprocess import call


def sort_bam(bam_file_full):
    print
    print 'Sort BAM file', bam_file_full, '...'
    stdout.flush()
    bam_full_filename = (splitext(bam_file_full))[0]
    output_filename = bam_full_filename + '_sorted.bam'
    output_base = bam_full_filename + '_sorted'
    command_line_list = ['time', 'samtools', 'sort', bam_file_full, output_base]
    command_line = 'time samtools sort ' + bam_file_full + ' ' + output_base
    print command_line
    stdout.flush()
    code = call(command_line_list)
    if code != 0:
        print 'Error: Something went wrong!'
        print 'Omit file', filename
        return ''
    else:
        print 'Done.'
        call(['rm', bam_file_full])
        print bam_file_full + ' was replaced with ' + output_filename
        stdout.flush()
        return output_filename


def create_bed_aligned(bam_file_full_sorted):
    print
    print 'Create BED file from sorted BAM file', bam_file_full_sorted, '...'
    stdout.flush()
    bam_sorted_full_filename = (splitext(bam_file_full_sorted))[0]
    bed_file_full = bam_sorted_full_filename + '.bed'
    command_line_list = ['time', 'bedtools', 'bamtobed', '-i', bam_file_full_sorted]
    command_line = 'time bedtools bamtobed ' + bam_file_full_sorted + ' > ' + bed_file_full
    print command_line
    stdout.flush()
    with open(bed_file_full, 'w') as out:
        code = call(command_line_list, stdout=out)
    if code != 0:
        print 'Error: Something went wrong!'
        print 'Omit file', bam_file_full_sorted
        return ''
    else:
        print 'Done.'
        print bed_file_full + ' was created.'
        stdout.flush()
        return bed_file_full


def create_tdf(bam_file_full_sorted, genome_tag, igvtools_directory, peaks_directory, \
               output_directory):
    print
    print 'Create TDF file from sorted BAM file', bam_file_full_sorted, '...'
    stdout.flush()
    tdf_filename = (splitext(basename(bam_file_full_sorted)))[0] + '.tdf'
    tdf_file_full = join(peaks_directory, tdf_filename)
    command_line_list = ['time', join(igvtools_directory, 'igvtools'), 'count', '-z', '5', \
                         '-w', '50', bam_file_full_sorted, tdf_file_full, genome_tag]
    command_line = 'time ' + join(igvtools_directory, 'igvtools') + ' count -z 5 -w 50 ' + \
                   bam_file_full_sorted + ' ' + tdf_file_full + ' ' + genome_tag
    print command_line
    stdout.flush()
    code = call(command_line_list)
    if code != 0:
        print 'Error: Something went wrong!'
        print 'Omit file', bam_file_full_sorted
        return False
    else:
        print 'Done.'
        print tdf_file_full + ' was created.'
    stdout.flush()
    return True
    
    
def call_peaks(bed_file_full, bam_file_full_sorted, sicer_directory, \
               input_directory, output_directory, genome_tag, \
               red_thr, window_size, fragment_size, egf, e_value):
    bed_file = basename(bed_file_full)
    mark_name = (splitext(bed_file))[0]
    # call peaks with SICER
    for gap_size in [200, 600]:
        print 
        print 'Call peaks with SICER, gap_size = ', gap_size
        stdout.flush()
        directory_name = mark_name
        out_dir = join(output_directory, directory_name + '_SICER_w' + str(window_size) + \
                                         '_g' + str(gap_size))
        makedirs(out_dir)
        command_line_list = ['time', join(sicer_directory, 'SICER-rb.sh'), input_directory, \
                             bed_file, out_dir, genome_tag, str(red_thr), str(window_size), \
                             str(fragment_size), str(egf), str(gap_size), str(e_value)]
        command_line = 'time ' + join(sicer_directory, 'SICER-rb.sh') + ' ' + input_directory + \
                       ' ' + bed_file + ' ' + out_dir + ' ' + genome_tag + ' ' + str(red_thr) + \
                       ' ' + str(window_size) + ' ' + str(fragment_size) + ' ' + str(egf) + \
                       ' ' + str(gap_size) + ' ' + str(e_value)
        print command_line
        stdout.flush()
        code = call(command_line_list)
        if code != 0:
            print 'Error: Something went wrong!'
            print 'Omit file', bed_file_full
        else:
            print 'Done.'
        stdout.flush()
    # call peaks with MACS
    print
    print 'Call peaks with MACS...'
    stdout.flush()
    prefix = mark_name + '_MACS'
    if genome_tag == 'hg19':
        macs_genome_tag = 'hs' # only hg19 tag is currently supported
    command_line_list = ['macs14', '-t', bam_file_full_sorted, '-f', 'BAM', \
                         '-g', macs_genome_tag, '-n', prefix]
    command_line = 'macs14 -t ' + bam_file_full_sorted + ' -f BAM -g ' + macs_genome_tag + \
                   ' -n ' + prefix
    print command_line
    stdout.flush()
    code = call(command_line_list)
    if code != 0:
        print 'Error: Something went wrong!'
        print 'Omit file', bam_file_full_sorted
        stdout.flush()
    else:
        print 'Done.'
        stdout.flush()
        macs_out_directory = join(output_directory, prefix)
        makedirs(macs_out_directory)
        working_directory = getcwd()
        call('mv ' + working_directory + '/' + prefix + '* ' + macs_out_directory, shell=True)
 

def copy_bed_files(output_directory, peaks_directory):
    print
    print 'Copy BED files with peaks to Peaks_BED_TDF...'
    stdout.flush()
    dir_list = listdir(output_directory)
    for dir in dir_list:
        if not isdir(join(output_directory, dir)):
            continue
        if 'SICER' in dir:
            command_line = 'cp ' + output_directory + '/' + dir + '/*scoreisland ' + \
                           peaks_directory
            print command_line
            stdout.flush()
            call(command_line, shell=True)
        elif 'MACS' in dir:
            command_line = 'cp ' + output_directory + '/' + dir + '/*peaks.bed ' + peaks_directory
            print command_line
            stdout.flush()
            call(command_line, shell=True)
    print 'Done.'
    stdout.flush()


if __name__ == '__main__':
    arguments = docopt(__doc__, version='call_peaks 0.1')
    input_directory = arguments["-d"].rstrip('/')
    if not exists(input_directory):
        print "Error: Can't find input directory: no such directory '" + \
              input_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(input_directory):
        print "Error: Input directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    sicer_directory = arguments["-s"].rstrip('/')
    if not exists(sicer_directory):
        print "Error: Can't find SICER directory: no such directory '" + \
              sicer_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(sicer_directory):
        print "Error: SICER directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    igvtools_directory = arguments["-g"].rstrip('/')
    if not exists(igvtools_directory):
        print "Error: Can't find IGVTools directory: no such directory '" + \
              igvtools_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(igvtools_directory):
        print "Error: IGVTools directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    output_directory = arguments["-o"].rstrip('/')
    if exists(output_directory) and not isdir(output_directory):
        print "Error: Output directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    if arguments["-i"] != None:
        genome_tag = arguments["-i"]
    else:
        genome_tag = 'hg19'

    # Only hg19 tag is currently supported
    if genome_tag != 'hg19':
        print 'Error: Only hg19 tag is currently supported. Exit.\n'
        sys.exit(1)

    if arguments["-r"] != None:
         try:
             red_thr = float(arguments["-r"])
         except ValueError:
             print "Error: Redundancy threshold must be an integer greater than 0. Exit.\n"
             sys.exit(1)
         if red_thr <= 0:
             print "Error: Redundancy threshold must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
        red_thr = 1

    if arguments["-r"] != None:
         try:
             red_thr = int(arguments["-r"])
         except ValueError:
             print "Error: Redundancy threshold must be an integer greater than 0. Exit.\n"
             sys.exit(1)
         if red_thr <= 0:
             print "Error: Redundancy threshold must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
        red_thr = 1

    if arguments["-w"] != None:
         try:
             window_size = int(arguments["-w"])
         except ValueError:
             print "Error: Window size must be an integer greater than 0. Exit.\n"
             sys.exit(1)
         if window_size <= 0:
             print "Error: Window size must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
        window_size = 200

    if arguments["-f"] != None:
         try:
             fragment_size = int(arguments["-f"])
         except ValueError:
             print "Error: Fragment size must be an integer greater than 0. Exit.\n"
             sys.exit(1)
         if fragment_size <= 0:
             print "Error: Fragment size must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
        fragment_size = 250

    if arguments["-l"] != None:
         try:
             egf = float(arguments["-l"])
         except ValueError:
             print "Error: Effective genome fraction must be a number between 0 and 1. Exit.\n"
             sys.exit(1)
         if egf <= 0 or egf > 1:
             print "Error: Effective genome fraction must be a number between 0 and 1. Exit.\n"
             sys.exit(1)
    else:
        egf = 0.77

    if arguments["-e"] != None:
         try:
             e_value = int(arguments["-e"])
         except ValueError:
             print "Error: E-value must be an integer greater than 0. Exit.\n"
             sys.exit(1)
         if e_value <= 0:
             print "Error: E-value must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
        e_value = 100


    print 'Check samtools...'
    stdout.flush()
    with open(devnull, 'w') as FNULL:
        check_samtools_code = call(['samtools', 'flags'], stderr=FNULL)
    if check_samtools_code != 0:
        print 'Error: samtools are not found. (It must be installed and available for ' + \
              "non-administrative user with command 'samtools'). Exit."
        sys.exit(1)
    print 'OK'
    stdout.flush()

    print 'Check bedtools...'
    stdout.flush()
    with open(devnull, 'w') as FNULL:
        check_bedtools_code = call(['bedtools', '--help'], stdout=FNULL)
    if check_bedtools_code != 0:
        print 'Error: bedtools are not found. (It must be installed and available for ' + \
              "non-administrative user with command 'samtools'). Exit."
        sys.exit(1)
    print 'OK'
    stdout.flush()

    print 'Check MACS v1.4...'
    stdout.flush()
    with open(devnull, 'w') as FNULL:
        check_macs14_code = call(['macs14', '--help'], stdout=FNULL)
    if check_macs14_code != 0:
        print 'Error: MACS v1.4 is not found. (It must be installed and available for ' + \
              "non-administrative user with command 'macs14'). Exit."
        sys.exit(1)
    print 'OK'
    stdout.flush()

    if not exists(output_directory):
        makedirs(output_directory)
    peaks_directory = join(output_directory, 'Peaks_BED_TDF')
    makedirs(peaks_directory)

    bam_files_list = [file for file in listdir(input_directory) if splitext(file)[1] == '.bam']
    print
    print 'BAM files to process:'
    for file in bam_files_list:
        print file
    print
    stdout.flush()
    for bam_file in bam_files_list:
        bam_file_full = join(input_directory, bam_file)
        bam_file_full_sorted = sort_bam(bam_file_full)
        if bam_file_full_sorted == '':
            continue
        bed_file_full = create_bed_aligned(bam_file_full_sorted)
        if bed_file_full == '':
            continue
        is_created = create_tdf(bam_file_full_sorted, genome_tag, igvtools_directory, \
                                peaks_directory, output_directory)
        if not is_created:
            continue
        call_peaks(bed_file_full, bam_file_full_sorted, sicer_directory, \
                   input_directory, output_directory, genome_tag, \
                   red_thr, window_size, fragment_size, egf, e_value)
        copy_bed_files(output_directory, peaks_directory)
    working_directory = getcwd()
    print 'Move chr.list and igv.log to output directory...'
    command_line_list = ['mv', join(working_directory, 'chr.list'), \
                         join(working_directory, 'igv.log'), output_directory]
    command_line = 'mv chr.list igv.log ' + output_directory
    print command_line
    stdout.flush()
    call(command_line_list)
    print
    print 'Processing is finished!'
    print
    stdout.flush()
