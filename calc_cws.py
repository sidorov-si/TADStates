#!/usr/bin/env python
"""
Calculate cross-window score (CWS) for each border between Hi-C windows 
(the width of a Hi-C window is equal to a contact matrix resolution, 
and a genome is split into back-to-back windows).

CWS of a border is a number of contacts that cross the border. If vicinity_size 
is set then during CWS calculation only those contacts are considered 
that connect regions (vicinity_size / 2) bp upstream and downstream 
the border. Otherwise, CWS is calculated across the whole chromosome.
Vicinity size must be a multiple of (2 * matrix_resolution).

For each chromosome a BED file is created with all border CWS values.
All BED files are also concatinated to the whole genome BED file.
Finally, CWS is plotted for each chromosome, and the graph is stored 
into correspondent PNG file.

Either a BED file with TADs or a BED file with TAD borders can be set
as input. In the latter case, border scores will be showed in the plot,
unless --no-labels option is set. Zero Avg CWS on Scores vs Avg CWS plot 
most probably means that there are no TAD borders with such score.

With -R option you can set a specific region within the chromosome. 
Only this region will be plotted. Coordinates must be set in bp and 
be multiples of matrix resolution.

Usage:
  calc_cws.py -m <contact_matrix> -r <matrix_resolution> [-c <chromosome_name> -n <track_name> -R <chromosome_region> (-b <BED_file_with_TAD_borders> [--labels]) -o <output_directory> -e <vicinity_size> -s <name_suffix>]
  calc_cws.py -d <input_directory> -r <matrix_resolution> [-B <directory_with_TAD_borders_BED_files> -N <track_name_for_whole_genome_BedGraph> -O <whole_genome_BedGraph_filename> -o <output_directory> -e <vicinity_size> -s <name_suffix>]

Options:
  -h --help                                  Show this screen.
  --version                                  Show version.
  -m <contact_matrix>                        File with contact matrix.
  -d <input_directory>                       Directory with contact matrices (.mat-files).
  -c <chromosome_name>                       Name of the chromosome. Determined from matrix file by default.
  -r <matrix_resolution>                     Matrix resolution.
  -e <vicinity_size>                         Vicinity size, bp. Default: CWS is calculated across the whole genome.
  -n <track_name>                            A name for the track for one chromosome. Default: 'chromosome_name'_CWS
  -R <chromosome_region>                     A region on a chromosome (a:b = [a,b]; :b = [matrix_resolution,b]; a: = [a, chromosome_end]; : = the whole chromosome; a, b are set in bp, a < b, a and b are multiples of matrix_resolution).
  -b <BED_file_with_TAD_borders>             BED file with scored TAD borders. 'Score' field must present a border score.
  -B <directory_with_TAD_borders_BED_files>  Directory with BED files of TAD borders (one file per chromosome).
  -N <track_name_for_whole_genome_BedGraph>  A track name for the whole genome BedGraph. Default: Whole_genome_CWS.
  -o <output_directory>                      Output directory name. Default: directory that contains this script.
  -O <whole_genome_BedGraph_filename>        Output whole genome BedGraph filename. Is also put in the output directory.
  -s <name_suffix>                           Suffix that will be added to all filenames and tracknames. It should contain only symbols allowed for filenames. Default: ''.
  --labels                                   Label each TAD border with its score.
"""

import sys

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


def bp_to_KMbp(bp_count):
    if bp_count % 1000 == 0:
        if bp_count % 1000000 == 0:
	    return str(bp_count / 1000000) + 'Mbp'
        else:
            return str(bp_count / 1000) + 'Kbp'
    else:
        return str(bp_count) + 'bp'


def plot_borders(start_score, end_score, color, ax, tad_border_cws, tad_border_scores, \
                 tad_border_coords):
    tad_border_cws = [cws_value for border_score, cws_value in \
                                zip(tad_border_scores, tad_border_cws) \
                                if border_score >= start_score and border_score <= end_score]
    tad_border_coords = [border_coord for border_score, border_coord in \
                                      zip(tad_border_scores, tad_border_coords) \
                                      if border_score >= start_score and border_score <= end_score]
    ax.plot(tad_border_coords, tad_border_cws, '.', color = color, ms = 10)
    if labels: # annotate each border with a score label
        label_values = [str(score) for score in tad_border_scores \
                        if score >= start_score and score <= end_score]
        for label, x, y in zip(label_values, tad_border_coords, tad_border_cws):
            ax.annotate(
                label, 
                xy = (x, y), xytext = (0, -15),
                textcoords = 'offset points', va = 'bottom'
            )


def autolabel(rects, ax):
    # attach some text labels to bars
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., height + 0.2, '%d' % int(height),
                ha = 'center', va = 'bottom')


def calc_cws(contact_matrix_filename, chrom_name, borders_filename, \
             whole_genome_analysis, last_chr):

    global wg_boxplot
    global wg_score_cws
    global wg_borders_in_mins
    global wg_borders_out_mins
    global wg_borders_in_vic_mins
    global wg_borders_out_vic_mins
    print
    print 'Contact matrix file:'
    print '   ', contact_matrix_filename
    if borders_filename != None:
        print 'TAD borders file:'
        print '   ', borders_filename
    print 'Matrix resolution:'
    print '   ', bp_to_KMbp(matrix_resolution)
    print 'BedGraph track name:'
    print '   ', track_name
    if vicinity_size != -1:
        print 'Vicinity size:'
        print '   ', bp_to_KMbp(vicinity_size)
    else:
        print 'Vicinity size: '
        print '    The whole chromosome.'
    stdout.flush()
    print 'Region to plot:'
    if start_coord == None and end_coord == None:
        print '    The whole chromosome.'
    else:
        print '    From', bp_to_KMbp(start_coord), 'to', bp_to_KMbp(end_coord)
    chrom_number = search(r'\d+|X|Y', chrom_name).group(0)
    if len(chrom_number) == 1 and chrom_number != 'X' and chrom_number != 'Y':
        chrom_number = '0' + chrom_number
        chrom_id = 'chr' + chrom_number
    else:
        chrom_id = chrom_name
    if vicinity_size != -1:
        vicinity_infix = bp_to_KMbp(vicinity_size)
    else:
        vicinity_infix = 'Whole'
    output_bedgraph_filename = join(bedgraph_directory, chrom_id + '_CWS' + '_vic' + \
                                    vicinity_infix + name_suffix + '.bedGraph')
    filename_list.append(output_bedgraph_filename)
    output_png_filename = join(png_directory, chrom_id + '_CWS' + '_vic' + \
                               vicinity_infix + name_suffix + '.png')
    print 'Output BedGraph file:'
    print '   ', output_bedgraph_filename
    print 'Output PNG file (CWS):'
    print '   ', output_png_filename
    if borders_filename != None:
        output_png_boxplot = join(png_directory, chrom_id + '_Scores-CWS' + '_vic' + \
                                  vicinity_infix + name_suffix + '.png') 
        output_png_avgplot = join(png_directory, chrom_id + '_Scores-CWS_avg' + '_vic' + \
                                  vicinity_infix + name_suffix + '.png')
        output_png_barplot = join(png_directory, chrom_id + '_Borders_in_mins' + '_vic' + \
                                  vicinity_infix + name_suffix + '.png')
        output_png_barplot_vic = join(png_directory, chrom_id + '_Borders_in_prox_mins' + \
                                      '_vic' + vicinity_infix + name_suffix + '.png')
        print 'Output PNG file (Scores vs CWS):'
        print '   ', output_png_boxplot
        print 'Output PNG file (Scores vs Avg CWS):'
        print '   ', output_png_avgplot
        print 'Output PNG file (Borders in CWS mins):'
        print '   ', output_png_barplot
        print 'Output PNG file (Borders in CWS mins proximities):'
        print '   ', output_png_barplot_vic
        stdout.flush()

    # Calculate CWS for all borders between windows
    with open(contact_matrix_filename, 'r') as infile:
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
                if j - k < 0:
                    left_rect_sum = 0
                    left_bottom_rect_sum = 0
                else:
                    left_rect_sum = cws[i, j - k]
                if i + k > n:
                    bottom_rect_sum = 0
                    left_bottom_rect_sum = 0
                else:
                    bottom_rect_sum = cws[i + k, j]
                if not j - k < 0 and not i + k > n:
                    left_bottom_rect_sum = cws[i + k, j - k]
                cws[i, j] = cws[i, j] - left_rect_sum - bottom_rect_sum + left_bottom_rect_sum

        # CWS for all borders between windows
        result = [cws[i, j] for i, j in zip(range(1, n + 1), range(0, n))]
        print 'Finish.'
        stdout.flush()
        
        # Generate BedGraph file with CWS for all borders between windows
        print 'Generate BedGraph file with CWS for chromosome', chrom_name, '...',
        stdout.flush()
        with open(output_bedgraph_filename, 'w') as dst:
            bedgraph_track_line = 'track type=bedGraph name=' + bedgraph_track_name
            dst.write(bedgraph_track_line + '\n')
            for border_number, cws_value in enumerate(result):
                border_name = str(cws_value)
                # Coordinates in BedGraph format are 0-based, 
                # and a region is presented by [x,y) interval.
                start_pos = (border_number + 1) * matrix_resolution - 1
                end_pos = (border_number + 1) * matrix_resolution + 1
                score = cws_value
                bedgraph_line = chrom_name + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + \
                           str(score)
                dst.write(bedgraph_line + '\n')
        print 'Finish.'
        stdout.flush()

        # Generate PNG file with CWS plot
        print 'Generate PNG file with CWS plot for chromosome', chrom_name, '...',
        stdout.flush()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if borders_filename == None:
            plot_header = 'CWS for ' + chrom_name
        else:
            plot_header = 'CWS for ' + chrom_name + ' with TAD borders'
        plot_header += '. Vicinity:'
        if vicinity_size != -1:
            plot_header += ' ' + bp_to_KMbp(vicinity_size)
        else:
            plot_header += ' whole ' + chrom_name
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
            plot_header += '. Region: ' + bp_to_KMbp(start_coord) + ' - ' + bp_to_KMbp(end_coord)
        print 'Finish.'
        stdout.flush()

        if borders_filename != None:
            # Plot TAD borders colored according to their scores
            message = 'Color '
            if labels:
                message += 'and label '
            message += 'TAD borders for chromosome ' + chrom_name + ' ...'
            print message,
            stdout.flush()
            with open(borders_filename, 'r') as borders:
                tad_border_coords = []
                tad_border_scores = []
                for i, line in enumerate(borders):
                    if i == 0:
                        continue # leave out the header
                    line_fields = line.rstrip('\n').split('\t')
                    border_coord = (int(line_fields[1]) + int(line_fields[2])) / 2
                    border_score = int(line_fields[4])
                    tad_border_coords.append(border_coord)
                    tad_border_scores.append(border_score)
                tad_border_numbers = [coord / matrix_resolution - 1 for coord in tad_border_coords]
                tad_border_cws = [result[border_number] for border_number in tad_border_numbers]
                # Select the weakest borders and paint them green
                plot_borders(1, 3, 'green', ax, tad_border_cws, tad_border_scores, tad_border_coords)
                # Select the borders with the middle strength and paint them yellow
                plot_borders(4, 6, 'yellow', ax, tad_border_cws, tad_border_scores, tad_border_coords)
                # Select the strong borders and paint them orange
                plot_borders(7, 9, 'orange', ax, tad_border_cws, tad_border_scores, tad_border_coords)
                # Select the strongest borders and paint them red
                plot_borders(10, 10, 'red', ax, tad_border_cws, tad_border_scores, tad_border_coords)
            print 'Finish.'
            stdout.flush()

        # Save CWS plot to file
        print 'Save CWS plot for', chrom_name, 'to file ...',
        stdout.flush()
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
        ax.cla()
        print 'Finish.'
        stdout.flush()
 
        if borders_filename != None:
            # Plot TAD border scores vs CWS (boxplot)
            print "Generate PNG file with 'TAD border scores vs CWS' plot for chromosome", \
                  chrom_name, "...",
            stdout.flush()
            boxplot_data = []
            for score in range(1, 11):
                boxplot_data.append([cws for cws, border_score in \
                                     zip(tad_border_cws, tad_border_scores) \
                                     if border_score == score])
            if whole_genome_analysis:
                if not wg_boxplot:
                    wg_boxplot = boxplot_data
                else:
                    wg_boxplot = [s1 + s2 for s1, s2 in zip(wg_boxplot, boxplot_data)]
            ax.boxplot(boxplot_data, 0, 'b.', whis = [5, 95])
            boxplot_header = 'TAD border scores vs CWS for ' + chrom_name + '. Vicinity:'
            if vicinity_size != -1:
                boxplot_header += ' ' + bp_to_KMbp(vicinity_size)
            else:
                boxplot_header += ' whole ' + chrom_name
            ax.set_xlabel('TAD border scores, bp')
            ax.set_ylabel('CWS')
            ax.set_title(boxplot_header)
            plt.savefig(output_png_boxplot)
            ax.cla()
            print 'Finish.'
            stdout.flush()
            if whole_genome_analysis and last_chr:
                print "Generate whole genome PNG file with 'TAD border scores vs CWS' plot...",
                ax.boxplot(wg_boxplot, 0, 'b.', whis = [5, 95])
                wg_boxplot_header = 'Whole genome TAD border scores vs CWS. Vicinity:'
                if vicinity_size != -1:
                    wg_boxplot_header += ' ' + bp_to_KMbp(vicinity_size)
                else:
                    wg_boxplot_header += ' whole chr'
                ax.set_xlabel('TAD border scores, bp')
                ax.set_ylabel('CWS')
                ax.set_title(wg_boxplot_header)
                plt.savefig(wg_output_png_boxplot)
                ax.cla()
                print 'Finish.'
                stdout.flush()               

        if borders_filename != None:
            # Plot TAD border scores vs avg CWS 
            print "Generate PNG file with 'TAD border scores vs avg CWS' plot for chromosome", \
                  chrom_name, "...",
            stdout.flush()
            avgplot_data = []
            chr_score_cws = []
            for score in range(1, 11):
                current_cws = [cws for cws, border_score in zip(tad_border_cws, tad_border_scores) \
                                   if border_score == score]
                if whole_genome_analysis:
                    chr_score_cws.append(current_cws)
                if len(current_cws) != 0:
                    avg_cws = sum(current_cws) / len(current_cws)
                else:
                    avg_cws = 0
                avgplot_data.append(avg_cws)
            if whole_genome_analysis:
                if not wg_score_cws:
                    wg_score_cws = chr_score_cws
                else:
                    wg_score_cws = [s1 + s2 for s1, s2 in zip(wg_score_cws, chr_score_cws)]
            ax.plot(range(1, 11), avgplot_data, '.-')
            avgplot_header = 'TAD border scores vs avg CWS for ' + chrom_name + '. Vicinity:'
            if vicinity_size != -1:
                avgplot_header += ' ' + bp_to_KMbp(vicinity_size)
            else:
                avgplot_header += ' whole ' + chrom_name
            ax.set_xlabel('TAD border scores, bp')
            ax.set_ylabel('Avg CWS')
            ax.set_title(avgplot_header)
            plt.savefig(output_png_avgplot)
            ax.cla()
            print 'Finish.'
            stdout.flush()
            if whole_genome_analysis and last_chr:
                print "Generate whole genome PNG file with 'TAD border score vs avg CWS' plot...",
                wg_avgplot_data = []
                for sublist in wg_score_cws:
                    score_avg_cws = sum(sublist) / len(sublist) if len(sublist) != 0 else 0
                    wg_avgplot_data.append(score_avg_cws)
                ax.plot(range(1, 11), wg_avgplot_data, '.-')
                wg_avgplot_header = 'Whole genome TAD border scores vs avg CWS. Vicinity:'
                if vicinity_size != -1:
                    wg_avgplot_header += ' ' + bp_to_KMbp(vicinity_size)
                else:
                    wg_avgplot_header += ' whole chr'
                ax.set_xlabel('TAD border scores, bp')
                ax.set_ylabel('Avg CWS')
                ax.set_title(wg_avgplot_header)
                plt.savefig(wg_output_png_avgplot)
                ax.cla()
                print 'Finish.'
                stdout.flush()

        if borders_filename != None:
            # Plot TAD border counts in CWS local minimums and out of them
            print 'Plot TAD border counts in CWS local minimums and out of them for chromosome',
            print chrom_name, '...',
            stdout.flush()
            with open(borders_filename, 'r') as borders:
                tad_border_coords = []
                tad_border_scores = []
                for i, line in enumerate(borders):
                    if i == 0:
                        continue # leave out the header
                    line_fields = line.rstrip('\n').split('\t')
                    border_coord = (int(line_fields[1]) + int(line_fields[2])) / 2
                    border_score = int(line_fields[4])
                    tad_border_coords.append(border_coord)
                    tad_border_scores.append(border_score)
                tad_border_numbers = [coord / matrix_resolution - 1 for coord in tad_border_coords]
                tad_border_cws = [result[border_number] for border_number in tad_border_numbers]
                # Count borders in CWS local minimums and out of them 
                borders_in_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
                borders_out_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
                for cws, score, number in zip(tad_border_cws, tad_border_scores, tad_border_numbers):
                    if number > 0 and number < len(result) - 1:
                        if cws < result[number - 1] and cws < result[number + 1]:
                            borders_in_mins[score] += 1
                        else:
                            borders_out_mins[score] += 1
                if whole_genome_analysis:
                    if not wg_borders_in_mins:
                        wg_borders_in_mins = borders_in_mins
                    else:
                        wg_borders_in_mins = {key: wg_borders_in_mins.get(key) + \
                                                   borders_in_mins.get(key) \
                                                   for key in set(wg_borders_in_mins)}
                    if not wg_borders_out_mins:
                        wg_borders_out_mins = borders_out_mins
                    else:
                        wg_borders_out_mins = {key: wg_borders_out_mins.get(key) + \
                                                   borders_out_mins.get(key) \
                                                   for key in set(wg_borders_out_mins)}
                # Plot borders_in_mins and borders_out_mins
                in_mins = tuple(value for key, value in borders_in_mins.items())
                max_value_1 = max(in_mins)
                ind = numpy.arange(10) # the x locations for the groups
                width = 0.35           # the width of the bars
                rects1 = ax.bar(ind, in_mins, width, color = 'r')
                out_mins = tuple(value for key, value in borders_out_mins.items())
                max_value_2 = max(out_mins)
                rects2 = ax.bar(ind + width, out_mins, width, color = 'y')
                max_value = max(max_value_1, max_value_2)
                ax.set_xlabel('Scores')
                ax.set_ylabel('Number of TAD borders')
                plt.ylim(ymax = 1.2 * max_value)
                barplot_header = 'TAD borders and CWS local mins for ' \
                                 + chrom_name + '. Vicinity:'
                if vicinity_size != -1:
                    barplot_header += ' ' + bp_to_KMbp(vicinity_size)
                else:
                    barplot_header += ' whole ' + chrom_name
                ax.set_title(barplot_header)
                ax.set_xticks(ind + width)
                ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
                ax.legend((rects1[0], rects2[0]), ('In local mins', 'Out of local mins'), \
                          loc='upper left')
                ax.legend()
                autolabel(rects1, ax)
                autolabel(rects2, ax)
                plt.savefig(output_png_barplot)
                ax.cla()
                plt.autoscale()
                print 'Finish.'
                stdout.flush()
                if whole_genome_analysis and last_chr:
                    print "Plot TAD border counts in CWS local minimums and out of them " + \
                          "for the whole genome ...",
                    wg_in_mins = tuple(value for key, value in wg_borders_in_mins.items())
                    wg_in_mins_max_value = max(wg_in_mins)
                    ind = numpy.arange(10) # the x locations for the groups
                    width = 0.35           # the width of the bars
                    wg_rects1 = ax.bar(ind, wg_in_mins, width, color = 'r')
                    wg_out_mins = tuple(value for key, value in wg_borders_out_mins.items())
                    wg_out_mins_max_value = max(wg_out_mins)
                    wg_rects2 = ax.bar(ind + width, wg_out_mins, width, color = 'y')
                    wg_mins_max_value = max(wg_in_mins_max_value, wg_out_mins_max_value)
                    plt.ylim(ymax = 1.2 * wg_mins_max_value)
                    ax.set_xlabel('Scores')
                    ax.set_ylabel('Number of TAD borders')
                    wg_barplot_header = 'TAD borders and CWS local mins ' + \
                                        'for the whole genome. Vicinity:'
                    if vicinity_size != -1:
                        wg_barplot_header += ' ' + bp_to_KMbp(vicinity_size)
                    else:
                        wg_barplot_header += ' whole chr'
                    ax.set_title(wg_barplot_header)
                    ax.set_xticks(ind + width)
                    ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
                    ax.legend((wg_rects1[0], wg_rects2[0]), ('In local mins', 'Out of local mins'), \
                              loc='upper left')
                    ax.legend()
                    autolabel(wg_rects1, ax)
                    autolabel(wg_rects2, ax)
                    plt.savefig(wg_output_png_barplot)
                    ax.cla()
                    plt.autoscale()
                    print 'Finish.'
                    stdout.flush()

            # Plot TAD border counts in some proximity of CWS local minimums and out of them
            print 'Plot TAD border counts in a matrix-resolution proximity of CWS local minimums ' + \
                  'and out of them for chromosome', chrom_name, '...',
            stdout.flush()
            with open(borders_filename, 'r') as borders:
                tad_border_coords = []
                tad_border_scores = []
                for i, line in enumerate(borders):
                    if i == 0:
                        continue # leave out the header
                    line_fields = line.rstrip('\n').split('\t')
                    border_coord = (int(line_fields[1]) + int(line_fields[2])) / 2
                    border_score = int(line_fields[4])
                    tad_border_coords.append(border_coord)
                    tad_border_scores.append(border_score)
                tad_border_numbers = [coord / matrix_resolution - 1 for coord in tad_border_coords]
                tad_border_cws = [result[border_number] for border_number in tad_border_numbers]
                # Count borders in matrix-resolution of CWS local minimums and out of them 
                borders_in_vic_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
                borders_out_vic_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
                for cws, score, number in zip(tad_border_cws, tad_border_scores, tad_border_numbers):
                    if number > 1 and number < len(result) - 2:
                        if cws < result[number - 1] and cws < result[number + 1] or \
                           result[number - 1] < result[number - 2] and result[number - 1] < cws or \
                           result[number + 1] < cws and result[number + 1] < result[number + 2]:
                            borders_in_vic_mins[score] += 1
                        else:
                            borders_out_vic_mins[score] += 1
                if whole_genome_analysis:
                    if not wg_borders_in_vic_mins:
                        wg_borders_in_vic_mins = borders_in_vic_mins
                    else:
                        wg_borders_in_vic_mins = {key: wg_borders_in_vic_mins.get(key) + \
                                                   borders_in_vic_mins.get(key) \
                                                   for key in set(wg_borders_in_vic_mins)}
                    if not wg_borders_out_vic_mins:
                        wg_borders_out_vic_mins = borders_out_vic_mins
                    else:
                        wg_borders_out_vic_mins = {key: wg_borders_out_vic_mins.get(key) + \
                                                   borders_out_vic_mins.get(key) \
                                                   for key in set(wg_borders_out_vic_mins)}

                # Plot borders_in_vic_mins and borders_out_vic_mins
                in_vic_mins = tuple(value for key, value in borders_in_vic_mins.items())
                max_value_vic_1 = max(in_vic_mins)
                ind = numpy.arange(10) # the x locations for the groups
                width = 0.35        # the width of the bars
                rects1 = ax.bar(ind, in_vic_mins, width, color = 'r')
                out_vic_mins = tuple(value for key, value in borders_out_vic_mins.items())
                max_value_vic_2 = max(out_vic_mins)
                rects2 = ax.bar(ind + width, out_vic_mins, width, color = 'y')
                max_vic_value = max(max_value_vic_1, max_value_vic_2)
                plt.ylim(ymax = 1.2 * max_vic_value)
                ax.set_xlabel('Scores')
                ax.set_ylabel('Number of TAD borders')
                barplot_header = 'TAD borders and CWS local mins proximities for ' \
                                 + chrom_name + '. Vicinity:'
                if vicinity_size != -1:
                    barplot_header += ' ' + bp_to_KMbp(vicinity_size)
                else:
                    barplot_header += ' whole ' + chrom_name
                ax.set_title(barplot_header)
                ax.set_xticks(ind + width)
                ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
                ax.legend((rects1[0], rects2[0]), \
                          ('In local mins proximities', 'Out of local mins proximities'), \
                          loc='upper left')
                ax.legend()
                autolabel(rects1, ax)
                autolabel(rects2, ax)
                plt.savefig(output_png_barplot_vic)
                ax.cla()
                plt.autoscale()
                print 'Finish.'
                stdout.flush()
                if whole_genome_analysis and last_chr:
                    print "Plot TAD border counts in a matrix-resolution proximity of CWS " + \
                          "local minimums and out of them for the whole genome ...",
                    wg_in_vic_mins = tuple(value for key, value in wg_borders_in_vic_mins.items())
                    max_vic_value_1 = max(wg_in_vic_mins)
                    ind = numpy.arange(10) # the x locations for the groups
                    width = 0.35           # the width of the bars
                    wg_vic_rects1 = ax.bar(ind, wg_in_vic_mins, width, color = 'r')
                    wg_out_vic_mins = tuple(value for key, value in wg_borders_out_vic_mins.items())
                    max_vic_value_2 = max(wg_out_vic_mins)
                    wg_vic_rects2 = ax.bar(ind + width, wg_out_vic_mins, width, color = 'y')
                    max_value_vic = max(max_vic_value_1, max_vic_value_2)
                    plt.ylim(ymax = 1.2 * max_value_vic)
                    ax.set_xlabel('Scores')
                    ax.set_ylabel('Number of TAD borders')
                    wg_vic_barplot_header = 'TAD borders and CWS local mins proximities' + \
                                            'for the whole genome. Vicinity:'
                    if vicinity_size != -1:
                        wg_vic_barplot_header += ' ' + bp_to_KMbp(vicinity_size)
                    else:
                        wg_vic_barplot_header += ' whole chr'
                    ax.set_title(wg_barplot_header)
                    ax.set_xticks(ind + width)
                    ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
                    ax.legend((wg_vic_rects1[0], wg_vic_rects2[0]), \
                              ('In local mins proximities', 'Out of local mins proximities'), \
                              loc='upper left')
                    ax.legend()
                    autolabel(wg_vic_rects1, ax)
                    autolabel(wg_vic_rects2, ax)
                    plt.savefig(wg_output_png_barplot_vic)
                    ax.cla()
                    plt.autoscale()
                    print 'Finish.'
                    stdout.flush()


def get_chrom_name(matrix_filename):
    with open(matrix_filename, 'r') as src:
        line_list = src.readline().strip().split('\t ')
        chrom_name = (line_list[0].split('_'))[0]
    return chrom_name


if __name__ == '__main__':
    arguments = docopt(__doc__, version='calc_cws 1.1')

    try:
        matrix_resolution = int(arguments["-r"])
    except ValueError:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)
    if matrix_resolution <= 0:
        print "Error: Matrix resolution must be an integer greater than 0. Exit.\n"
        sys.exit(1)

    if arguments["-e"] != None:
        try:
            vicinity_size = int(arguments["-e"])
        except ValueError:
            print "Error: Vicinity size must be an integer greater than 0. Exit.\n"
            sys.exit(1)
        if vicinity_size <= 0:
            print "Error: Vicinity size must be an integer greater than 0. Exit.\n"
            sys.exit(1)
        if not vicinity_size % (2 * matrix_resolution) == 0:
            sys.exit(1)
    else:
        vicinity_size = -1

    if arguments["-s"] != None:
        name_suffix = arguments["-s"]
    else:
        name_suffix = ''

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

        borders_filename = None
        labels = None
        if arguments["-b"] != None:
            borders_filename = arguments["-b"]
            if not exists(borders_filename):
                print "Error: Can't find BED file with TAD borders: no such file '" + \
                      borders_filename + "'. Exit.\n"
                sys.exit(1)
            if not isfile(borders_filename):
                print "Error: BED file with TAD borders must be a regular file. " + \
                      "Something else given. Exit.\n"
                sys.exit(1)
            labels = arguments["--labels"]

        input_directory = None
        borders_directory = None
        output_wg_bedgraph_filename = None
        all_track_name = None
    else:
        matrix_filename = None
        chrom_name = None
        borders_filename = None
        labels = None
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
            all_track_name = "All_border_CWS" + '_vic' +  bp_to_KMbp(vicinity_size) + name_suffix 
        if arguments["-O"] != None:
            output_wg_bedgraph_filename = arguments["-O"]
        else:
            output_wg_bedgraph_filename = 'All_borders_CWS' + '_vic' + \
                                          bp_to_KMbp(vicinity_size) + name_suffix + '.bedGraph'
        borders_directory = None
        if arguments["-B"] != None:
            borders_directory = arguments["-B"]
            if not exists(borders_directory):
                print "Error: Can't find directory with TAD border BED files: " + \
                      "no such directory '" + borders_directory + "'. Exit.\n"
                sys.exit(1)
            if not isdir(borders_directory):
                print "Error: Directory with TAD border BED files must be a directory:). " + \
                      "Something else given. Exit.\n"
                sys.exit(1)

    if arguments["-o"] != None:
        output_directory = arguments["-o"].rstrip('/')
    else:
        if input_directory != None:
            output_directory = input_directory + '_CWS'
        else:
            output_directory = ''

    if arguments["-O"] == None and output_wg_bedgraph_filename != None:
        output_wg_bedgraph_filename = join(output_directory, output_wg_bedgraph_filename)

    filename_list = []
    if output_directory != '':
        if not exists(output_directory):
            makedirs(output_directory)
    bedgraph_directory = join(output_directory, 'BedGraph_CWS')
    png_directory = join(output_directory, 'PNG_CWS')
    if not exists(bedgraph_directory):
        makedirs(bedgraph_directory)
    if not exists(png_directory):
        makedirs(png_directory)
    if borders_directory != None:
        wg_output_png_boxplot = join(png_directory, 'Scores-CWS' + '_vic' + \
                                     bp_to_KMbp(vicinity_size) + name_suffix + '.png')
        wg_output_png_avgplot = join(png_directory, 'Scores-CWS_avg' + '_vic' + \
                                     bp_to_KMbp(vicinity_size) + name_suffix + '.png')
        wg_output_png_barplot = join(png_directory, 'Borders_in_mins' + '_vic' + \
                                     bp_to_KMbp(vicinity_size) + name_suffix + '.png')
        wg_output_png_barplot_vic = join(png_directory, 'Borders_in_prox_mins' + \
                                         '_vic' + bp_to_KMbp(vicinity_size) + name_suffix + '.png')

    print
    if input_directory != None:
        print 'Input directory:'
        print '   ', input_directory
        stdout.flush()
    if borders_filename != None:
        print 'BED file with TAD borders:'
        print '   ', borders_filename
        stdout.flush()
    if borders_directory != None:
        print 'Directory with TAD border BED files:'
        print '   ', borders_directory
        stdout.flush()
    if output_directory != None and output_directory != '':
        print 'Output directory:'
        print '   ', output_directory
        stdout.flush() 
    if name_suffix != '':
        print 'Suffix for all filenames and tracknames:'
        print '   ', name_suffix
        stdout.flush() 
    if output_wg_bedgraph_filename != None:
        print 'Output whole genome BedGraph file:'
        print '   ', output_wg_bedgraph_filename
        stdout.flush()
    if all_track_name != None:
        print 'Whole genome BedGraph track name:'
        print '   ', all_track_name
        stdout.flush()
    if borders_directory != None:
        print 'Whole genome output PNG file (Scores vs CWS):'
        print '   ', wg_output_png_boxplot
        print 'Whole genome output PNG file (Scores vs Avg CWS):'
        print '   ', wg_output_png_avgplot
        print 'Whole genome output PNG file (Borders in CWS mins):'
        print '   ', wg_output_png_barplot
        print 'Whole genome output PNG file (Borders in CWS mins proximities):'
        print '   ', wg_output_png_barplot_vic
        stdout.flush() 

    if matrix_filename != None: # there is only one contact matrix
        if chrom_name == None:
            chrom_name = get_chrom_name(matrix_filename)
        if track_name == None:
            track_name = chrom_name + '_CWS' + '_vic' + bp_to_KMbp(vicinity_size) + name_suffix
        bedgraph_track_name = track_name + '_bedGraph'
        whole_genome_analysis = False
        last_chr = False
        calc_cws(matrix_filename, chrom_name, borders_filename, whole_genome_analysis, last_chr)
    else: # there is a directory with matrices
        print
        print 'Calculate CWS for all chromosomes in the input directory...'
        stdout.flush()
        wg_boxplot = []
        wg_score_cws = []
        wg_borders_in_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        wg_borders_out_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        wg_borders_in_vic_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        wg_borders_out_vic_mins = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
        whole_genome_analysis = True
        matrix_filenames = sorted(listdir(input_directory))
        borders_filenames = sorted(listdir(borders_directory))
        filenames_zip = zip(matrix_filenames, borders_filenames)
        for number, (matrix_filename, borders_filename) in enumerate(filenames_zip):
            matrix_filename_full = join(input_directory, matrix_filename)
            borders_filename_full = join(borders_directory, borders_filename)
            chrom_name = get_chrom_name(matrix_filename_full)
            track_name = chrom_name + '_CWS' + '_vic' + bp_to_KMbp(vicinity_size) + name_suffix
            bedgraph_track_name = track_name + '_bedGraph'
            last_chr = True if number == len(filenames_zip) - 1 else False
            calc_cws(matrix_filename_full, chrom_name, borders_filename_full, \
                     whole_genome_analysis, last_chr)
        print 
        print 'All chromosomes are processed.'
        stdout.flush()
        # merge BedGraph files for individual chromosomes in one BedGraph file
        print 'Generate whole genome BedGraph file with CWS...',
        stdout.flush()
        with open(output_wg_bedgraph_filename, 'w') as dst:
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

