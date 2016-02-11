#!/usr/bin/env python
"""
Generate contact matrices from Hi-C reads in FASTQ format.
If quality assessment and trimming are necessary for Hi-C reads, they must be 
made before launching the script on these reads.

Usage:
  make_matrices.py (-1 <left_fastq_reads> -2 <right_fastq_reads> -p <reads_prefix> | -r <fastq_reads>) -g <fasta_ref_genome> -i <ref_genome_index> -o <output_directory> [-e <restriction_enzyme> -c <list_of_chromosome_names> -R <matrix_resolution> -t <threads_number> --clean --tmp-dir <temp_directory>]

Options:
  -h --help                      Show this screen.
  --version                      Show version.
  -1 <left_fastq_reads>          FASTQ file with left (forward) Hi-C reads.
  -2 <right_fastq_reads>         FASTQ file with right (reverse) Hi-C reads.
  -p <output_prefix>             Prefix for output files with information about both left and right Hi-C reads.
  -r <fastq_reads>               FASTQ file with both left (forward) and right (reverse) Hi-C reads.
  -g <fasta_ref_genome>          FASTA file with reference genome.
  -i <ref_genome_index>          Reference genome index for GEMtools.
  -e <restriction_enzyme>        Name of the restriction enzyme. Default: HindIII. The full list of possible restriction enzyme names see https://github.com/3DGenomes/tadbit/blob/master/_pytadbit/mapping/restriction_enzymes.py
  -c <list_of_chromosome_names>  Text file with chromosome names (one name per line, without empty spaces). Default: chr1, ..., chr22, chrX, chrY (Homo sapiens kariotype).
  -R <matrix_resolution>         Resolution of the output Hi-C matrices (bp). Default: 100000.
  -o <output_directory>          Output directory.
  -t <threads_number>            Number of threads for read mapping. Default: 8.
  --clean                        Remove all SAM and TSV files from the output directory (they're necessary for only the counstruction of contact matrices).
  --tmp-dir <temp-directory>     Temp directory. Default: /tmp.
"""

# Some code from the TADbit tutorial is used here

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
from pytadbit.mapping.mapper import iterative_mapping
from pytadbit.parsers.sam_parser import parse_sam
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.mapping.mapper import get_intersection
from pytadbit.mapping.filter import filter_reads
from pytadbit.mapping.filter import apply_filter
from pytadbit.mapping.analyze import hic_map
from os.path import basename
from os.path import splitext
from os.path import join
from os.path import exists
from os.path import isfile
from os import rename
from os import remove
from sys import stdout


def calc_left_right_ranges(reads_fastq): 
    # called when both left and right reads are stored in one file
    with open(reads_fastq, 'r') as src:
        read_id = src.readline()
        read_seq = src.readline().rstrip('\n')
        read_len = len(read_seq) / 2
        print 'Read lenght:', read_len
        stdout.flush()
        min_len = 20
        step = 5
        range_len = (read_len - min_len + 1) / step + 1
        # calc range_start_left
        range_start_left = [1] * range_len
        # calc range_stop_left
        first_left_stop_pos = min_len
        last_left_stop_pos = first_left_stop_pos + (range_len - 1) * step
        range_stop_left = [pos for pos in range(first_left_stop_pos, last_left_stop_pos + 1, step)]
        # calc range_start_right
        right_start_pos = read_len + 1
        range_start_right = [right_start_pos] * range_len
        # calc range_stop_right
        first_right_stop_pos = right_start_pos + min_len - 1
        last_right_stop_pos = first_right_stop_pos + (range_len - 1) * step
        range_stop_right = [pos for pos in range(first_right_stop_pos, \
                            last_right_stop_pos + 1, step)]
        return range_start_left, range_stop_left, range_start_right, range_stop_right


def calc_range(left_reads_fastq): 
    # called when reads are separated in two files
    # to calc ranges we need only one of the two files here
    with open(left_reads_fastq, 'r') as src:
        read_id = src.readline()
        read_seq = src.readline().rstrip('\n')
        read_len = len(read_seq)
        print 'Read lenght:', read_len
        stdout.flush()
        min_len = 20
        step = 5
        range_len = (read_len - min_len + 1) / step + 1
        # calc range_start
        range_start = [1] * range_len
        # calc range_stop
        first_stop_pos = min_len
        last_stop_pos = first_stop_pos + (range_len - 1) * step
        range_stop = [pos for pos in range(first_stop_pos, last_stop_pos + 1, step)]
        range_start_left = range_start
        range_start_right = range_start[:]
        range_stop_left = range_stop
        range_stop_right = range_stop[:]
        return range_start_left, range_stop_left, range_start_right, range_stop_right


def count_lines(filename):
    with open(filename, 'r') as src:
        for i, line in enumerate(src):
            pass
    return i + 1


def add_resolution(chromosomes, resolution, output_directory):
    stdout.flush()
    for chr in chromosomes:
        chr_filename = join(output_directory, chr + '.mat')
        chr_new_filename = join(output_directory, chr + '_' + str(resolution) + '.mat')
        rename(chr_filename, chr_new_filename)
        stdout.flush()


def add_headers(chromosomes, resolution, output_directory):
    stdout.flush()
    for chr in chromosomes:
        chr_filename = join(output_directory, chr + '_' + str(resolution) + '.mat')
        matrix_size = count_lines(chr_filename)
        id_list = [chr + '_' + str(i) for i in range(1, matrix_size + 1)]
        header = '\t'.join(id_list)
        chr_filename_tmp = chr_filename + '.tmp'
        with open(chr_filename, 'r') as src, open(chr_filename_tmp, 'w') as dst:
            stdout.flush()
            dst.write(header + '\n')
            for id, line in zip(id_list, src):
                dst.write(id + '\t' + line)
        rename(chr_filename_tmp, chr_filename)


def make_matrices(left_reads_fastq, right_reads_fastq, reads_fastq, genome_fasta, genome_index, \
                  output_directory, output_prefix, enzyme, res, chromosomes, threads_number, \
                  clean_tmp, tmp_dir):

    print 'Begin to process reads.'

    left_reads = ''
    right_reads = ''
    if reads_fastq != '': # left and right reads are stored in one file
        range_start_left, range_stop_left, \
        range_start_right, range_stop_right = calc_left_right_ranges(reads_fastq)
        print 'Reads:                     ', reads_fastq
        left_reads = reads_fastq
        right_reads = reads_fastq
    else: # left and right reads are stored separately
        range_start_left, range_stop_left, \
        range_start_right, range_stop_right = calc_range(left_reads_fastq)
        print 'Left reads:                ', left_reads_fastq
        print 'Right reads:               ', right_reads_fastq
        print 'Output prefix:             ', output_prefix
        left_reads = left_reads_fastq
        right_reads = right_reads_fastq

    print 'Reference genome FASTA:    ', genome_fasta
    print 'Reference genome GEM index:', genome_index
    print 'Output directory:          ', output_directory
    print 'Temp directory:            ', tmp_dir
    print 'Enzyme:                    ', enzyme
    print 'Resolution:                ', res, 'bp'
    print 'Number of threads:         ', threads_number
    print 'Start pos for left reads:  ', range_start_left
    print 'Stop pos for left reads:   ', range_stop_left
    print 'Start pos for right reads: ', range_start_right
    print 'Stop pos for right reads:  ', range_stop_right
    stdout.flush()

    # map left reads to reference genome
    out_sam_left_name = splitext(basename(left_reads))[0] + '_left.sam'
    out_sam_left_path = join(output_directory, out_sam_left_name)
    print 'Iterative mapping of left reads (using ' + str(threads_number) + ' threads)...'
    stdout.flush()
    sams_left = iterative_mapping(genome_index, left_reads, out_sam_left_path, \
                                  range_start_left, range_stop_left, nthreads=threads_number,
                                  temp_dir=tmp_dir)
    print 'Done.'
    stdout.flush()

    # map right reads to reference genome
    out_sam_right_name = splitext(basename(right_reads))[0] + '_right.sam'
    out_sam_right_path = join(output_directory, out_sam_right_name)
    print 'Iterative mapping of right reads (using ' + str(threads_number) + ' threads)...'
    stdout.flush()
    sams_right = iterative_mapping(genome_index, right_reads, out_sam_right_path, \
                                   range_start_right, range_stop_right, nthreads=threads_number,
                                   temp_dir=tmp_dir)
    print 'Done.'
    stdout.flush()

    # load reference genome sequence
    print 'Load reference genome sequence...'
    stdout.flush()
    chroms = chromosomes[:]
    genome_seq = parse_fasta(genome_fasta, chr_names=chroms)
    print 'Done.'
    stdout.flush()

    # create files with information about every left and right read 
    # and about their placement with respect to restriction sites
    tsv_left_name = splitext(basename(left_reads))[0] + '_left.tsv'
    tsv_left = join(output_directory, tsv_left_name)
    tsv_right_name = splitext(basename(right_reads))[0] + '_right.tsv'
    tsv_right = join(output_directory, tsv_right_name)
    print 'Get information about restriction sites and reads placement...'
    stdout.flush()
    parse_sam(sams_left, sams_right, tsv_left, tsv_right, genome_seq, enzyme, \
              verbose=True, ncpus=8)
    print 'Done.'
    stdout.flush()

    # create file with both left and right reads that uniquelly mapped to reference genome
    if reads_fastq != '': # left and right reads are stored in one file
        common_reads_prefix = splitext(basename(reads_fastq))[0]
    else: # left and right reads are stored separately
        common_reads_prefix = output_prefix
    uniq_reads_name = common_reads_prefix + '_both_map_uniq.tsv'
    uniq_reads = join(output_directory, uniq_reads_name)
    print 'Merge info about left and right reads in one file...'
    stdout.flush()
    get_intersection(tsv_left, tsv_right, uniq_reads, verbose=True)
    print 'Done.'
    stdout.flush()

    # find read IDs that are filtered by default TADbit filters
    print 'Mask reads...'
    stdout.flush()
    masked = filter_reads(uniq_reads)
    print 'Done.'
    stdout.flush()

    # apply all filters (exclude reads that were filtered)
    print 'Filter masked reads...'
    stdout.flush()
    filtered_reads_name = common_reads_prefix + '_filtered.tsv'
    filtered_reads = join(output_directory, filtered_reads_name)
    apply_filter(uniq_reads, filtered_reads, masked)
    print 'Done.'
    stdout.flush()

    # create matrices (one matrix per chromosome)
    print 'Create Hi-C maps (one per chromosome)...'
    stdout.flush()
    hic_map(filtered_reads, resolution=res, by_chrom='intra', savedata=output_directory)
    print 'Done.'
    stdout.flush()
    print 'Add resolution (' + str(resolution) + ') to matrix filenames...'
    stdout.flush()
    add_resolution(chromosomes, resolution, output_directory)
    print 'Done.'
    stdout.flush()
    print 'Add headers to matrix files...'
    stdout.flush()
    add_headers(chromosomes, resolution, output_directory)
    print 'Done.'
    stdout.flush()
    if clean_tmp: # Remove all SAM and TSV files from the output directory
        print 'Remove SAM and TSV files from the output directory.'
        stdout.flush()
        remove(out_sam_left_path + '*')
        remove(out_sam_right_path + '*')
        remove(join(output_directory, '*.tsv'))
        print 'Done.'
        stdout.flush()


def get_chromosome_names(chromosome_names_file):
    chrom_names = []
    if not exists(chromosome_names_file):
        print "Error: Can't find text file with chromosome names: no such file '" + \
              chromosome_names_file + "'. Exit.\n"
        sys.exit(1)
    if not isfile(chromosome_names_file):
        print "Error: File with chromosome names must be a regular file. " + \
              "Something else given. Exit.\n"
        sys.exit(1)
    with open(chromosome_names_file, 'r') as src:
         for line in src:
             name = line.rstrip('\n')
             chrom_names.append(name)
    return chrom_names


if __name__ == '__main__':
    arguments = docopt(__doc__, version='make_matrices 0.8')
    if arguments["-r"] != None: # left and right reads are stored in one file
        reads_fastq = arguments["-r"]
        if not exists(reads_fastq):
            print "Error: Can't find FASTQ file with reads: no such file '" + \
                  reads_fastq + "'. Exit.\n"
            sys.exit(1)
        if not isfile(reads_fastq):
            print "Error: FASTQ file with reads must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
        left_reads_fastq = ''
        right_reads_fastq = ''
        output_prefix = ''
    else: # left and right reads are stored separately
        left_reads_fastq = arguments["-1"]
        if not exists(left_reads_fastq):
            print "Error: Can't find FASTQ file with left (forward) reads: no such file '" + \
                  left_reads_fastq + "'. Exit.\n"
            sys.exit(1)
        if not isfile(left_reads_fastq):
            print "Error: FASTQ file with left (forward) reads must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)

        right_reads_fastq = arguments["-2"]
        if not exists(right_reads_fastq):
            print "Error: Can't find FASTQ file with right (reverse) reads: no such file '" + \
                  right_reads_fastq + "'. Exit.\n"
            sys.exit(1)
        if not isfile(right_reads_fastq):
            print "Error: FASTQ file with reads must be a regular file. " + \
                  "Something else given. Exit.\n"
            sys.exit(1)
        output_prefix = arguments["-p"]
        reads_fastq = ''

    genome_fasta = arguments["-g"]
    if not exists(genome_fasta):
        print "Error: Can't find FASTA file with reference genome sequence: " + \
              "no such file '" + genome_fasta + "'. Exit.\n"
        sys.exit(1)
    if not isfile(genome_fasta):
        print "Error: FASTA file with reference genome sequence must be a regular file. " + \
              "Something else given. Exit.\n"
        sys.exit(1)

    genome_index = arguments["-i"]
    if not exists(genome_index):
        print "Error: Can't find GEM index file for reference genome: no such file '" + \
              genome_index + "'. Exit.\n"
        sys.exit(1)
    if not isfile(genome_index):
        print "Error: GEM index file for reference genome must be a regular file. " + \
        "Something else given. Exit.\n"
        sys.exit(1)

    output_directory = arguments["-o"].rstrip('/')

    if arguments["-e"] != None:
         enzyme = arguments["-e"]
    else:
         enzyme = 'HindIII'

    if arguments["-R"] != None:
         try:
             resolution = int(arguments["-R"])
         except ValueError:
             print "Error: Resolution must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
         resolution = 100000

    if arguments["-c"] != None:
         chromosomes = get_chromosome_names(arguments["-c"]) # Exit if file not found
    else:
         chromosomes = ['chr' + str(c) for c in range(1, 23) + ['X', 'Y']]

    if arguments["-t"] != None:
         try:
             threads_number = int(arguments["-t"])
         except ValueError:
             print "Error: Thread number must be an integer greater than 0. Exit.\n"
             sys.exit(1)
    else:
         threads_number = 8

    if arguments["--clean"]:
        clean_tmp = True
    else:
        clean_tmp = False

    if arguments["--tmp-dir"] != None:
        tmp_dir = arguments["--tmp-dir"]
    else:
        tmp_dir = "/tmp"

    make_matrices(left_reads_fastq, right_reads_fastq, reads_fastq, genome_fasta, genome_index, \
                  output_directory, output_prefix, enzyme, resolution, chromosomes, threads_number, \
                  clean_tmp, tmp_dir)

