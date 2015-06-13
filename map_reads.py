#!/usr/bin/env python
"""
Align ChIP-Seq reads on reference genome.

Usage:
  map_reads.py -d <input_directory> -i <genome_index> -o <output_directory> [-t <thread_number>] 

Options:
  -h --help              Show this screen.
  --version              Show version.
  -d <input_directory>   Directory with FASTQ files containing single reads for all marks to be processed.
  -i <genome_index>      Bowtie index prefix (with path) for reference genome.
  -o <output_directory>  Output directory name.
  -t <thread_number>     Number of threads for bowtie. Default: 4.

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
from sys import stdout
from subprocess import call


def map_reads(filename_full, filename, ref_index, output_directory):
    print
    print 'Processing file', filename, '...'
    stdout.flush()
    output_filename = join(output_directory, filename + '.sam')
    sam_files_list.append(output_filename)
    command_line_list = ['bowtie', '-m', '1', '-v', '3', '--best', '--strata', '-t', '-S', '-p', \
                         '32', ref_index, filename_full]
    command_line = 'bowtie -m 1 -v 3 --best --strata -t -S -p 32 ' + ref_index + ' ' + \
                   filename_full + ' > ' + output_filename
    print 'Start alignment:'
    print command_line
    stdout.flush()
    with open(output_filename, 'w') as out:
        code = call(command_line_list, stdout=out)
    if code != 0:
        print 'Error: Something went wrong!'
        print 'Omit file', filename
    else:
        print 'Done.'
        print 'Alignment was written to', output_filename
    stdout.flush()


def sam_to_bam(sam_files_list):
    print 
    stdout.flush()
    print 'Convert SAM files to BAM...'
    for sam_file_full in sam_files_list:
        bam_file_full = (splitext(sam_file_full))[0] + '.bam'
        command_line_list = ['time', 'samtools', 'view', '-bS', '-o',  bam_file_full, \
                             sam_file_full]
        command_line = 'time samtools view -bS -o ' + bam_file_full + ' ' + sam_file_full
        print 'Convert ', sam_file_full
        stdout.flush()
        code = call(command_line_list)
        if code != 0:
            print 'Error: Something went wrong!'
            print 'Omit file', sam_file_full
            stdout.flush()
        else:
            print 'Done.'
            print 'File', sam_file_full, ' was converted to', bam_file_full
            bam_files_list.append(bam_file_full)
            command_line_list_rm = ['rm', sam_file_full]
            command_line = 'rm ' + sam_file_full
            print 'Remove', sam_file_full, ':'
            stdout.flush()
            call(command_line_list_rm)

def calc_map_stat(bam_files_list, output_directory):
    print
    stdout.flush()
    stat_filename = join(output_directory, 'align_stats.txt')
    print 'Calculate alignment statistics...'
    stdout.flush()
    for bam_file_full in bam_files_list:
        command_line_list = ['samtools', 'flagstat', bam_file_full]
        command_line = 'samtools flagstat ' + bam_file_full + ' > ' + stat_filename
        print 'Calculate statistics for ' + bam_file_full + ':'
        print command_line
        stdout.flush()
        with open(stat_filename, 'a') as out:
            out.write(bam_file_full + ':\n')
            out.write('-' * (len(bam_file_full) + 1) + '\n')
            out.flush()
            code = call(command_line_list, stdout=out)
            out.write('\n')
        if code != 0:
            print 'Error: Something went wrong!'
            print 'Omit file', file
        else:
            print 'OK.'
    print 'Done.'
    print 'Alignment statistics for all marks were written to ' + stat_filename
    stdout.flush()


if __name__ == '__main__':
    arguments = docopt(__doc__, version='map_reads 0.1')
    input_directory = arguments["-d"].rstrip('/')
    if not exists(input_directory):
        print "Error: Can't find input directory: no such directory '" + \
              input_directory + "'. Exit.\n"
        sys.exit(1)
    if not isdir(input_directory):
        print "Error: Input directory must be a directory:). Something else given. Exit.\n"
        sys.exit(1)

    ref_index = arguments["-i"]

    output_directory = arguments["-o"].rstrip('/')
    if exists(output_directory) and not isdir(output_directory):
        print "Error: Output directory must be a directory:). Something else given. Exit.\n"
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

    print 'Check bowtie...'
    stdout.flush()
    with open(devnull, 'w') as FNULL:
        check_bowtie_code = call(['bowtie', '-h'], stdout=FNULL)
    if check_bowtie_code != 0:
        print 'Error: bowtie is not found. (It must be installed and available for ' + \
              "non-administrative user with command 'bowtie'). Exit."
        sys.exit(1)
    print 'OK'
    stdout.flush()

    print 'Check samtools...'
    no_stat = False
    stdout.flush()
    with open(devnull, 'w') as FNULL:
        check_samtools_code = call(['samtools', 'flags'], stderr=FNULL)
    if check_samtools_code != 0:
        print 'Error: samtools are not found. (It must be installed and available for ' + \
              "non-administrative user with command 'samtools'). Exit."
        sys.exit(1)
    print 'OK'
    stdout.flush()

    if not exists(output_directory):
        makedirs(output_directory)

    sam_files_list = []
    bam_files_list = []
    files = listdir(input_directory)
    for file in files:
        filename = splitext(basename(file))[0]        
        filename_full = join(input_directory, file)
        map_reads(filename_full, filename, ref_index, output_directory)
    sam_to_bam(sam_files_list)
    calc_map_stat(bam_files_list, output_directory)
