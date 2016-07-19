#!/usr/bin/env python

#
# Copyright 2016, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType


# 1) Align 20 million DNA-seq 101bp reads with HISAT2
#    time ../../../hisat2 -p 3 -x ../../indexes/HISAT2/genome -U 20M_1.fq -S alignment.sam
#    real 9m9.279s
#
# 2) Sort alignments using samtools
#
#    time samtools sort -@ 3 -o alignment.bam alignment.sam
#    real 3m31.291s
#
#    time ../../../hisat2 -p 3 -x ../../indexes/HISAT2/genome -U 20M_1.fq | samtools sort -@ 3 -o alignment.bam -
#    real 12m6.543s
#
# 3) Sort alignments with my methods hisat2_sort_alignments.py
#
#    time ../../../hisat2 -p 3 -x ../../indexes/HISAT2/genome -U 20M_1.fq | ./hisat2_sort_alignments.py -p 3
#    samtools view -bS -


"""
"""
def sort_alignments(input_fname,
                    temp_fname,
                    threads,
                    interval_size,
                    verbose):
    seq_list = []
    seq2file = {}
    first_alignment = True
    if input_fname == '-':
        input_file = sys.stdin
    else:
        input_file = open(input_fname)
    for line in input_file:
        if line.startswith("@"):
            print line,
            if line.startswith("@SQ"):
                _, seq_name, seq_length = line.strip().split()
                assert len(seq_name) > 3 and len(seq_length) > 3
                seq_name, seq_length = seq_name[3:], int(seq_length[3:])
                seq_list.append([seq_name, seq_length])
            continue
        if first_alignment:
            first_alignment = False
            seq_list.append(["unmapped", 0])
            for seq_name, seq_length in seq_list:
                tmp_fname = "%s.%s.tmp.sam" % (temp_fname, seq_name)
                tmp_file = open(temp_fname, 'w')
                seq2file[seq_name] = [temp_fname, tmp_file]
            
        read_name, flag, seq_name, pos = line.split()[:4]
        flag, pos = int(flag), int(pos)
        if flag & 0x4 != 0:
            print >> seq2file["unmapped"][1], line,
        else:
            print >> seq2file[seq_name][1], line,

    for seq_name, tmp_file in seq2file.items():

    for line in header:
        print line,


    for seq_name, tmp_file in seq2
    # Sort alignments
    """
    sambam_cmd = ["sort",
                  "view",
                  "-bS",
                  "-"]
    sambam_proc = subprocess.Popen(sambam_cmd,
                                   stdin=align_proc.stdout,
                                   stdout=open("hla_input_unsorted.bam", 'w'),
                                   stderr=open("/dev/null", 'w'))
    sambam_proc.communicate()
    bamsort_cmd = ["samtools",
                   "sort",
                   "hla_input_unsorted.bam",
                   "-o", "hla_input.bam"]
    bamsort_proc = subprocess.Popen(bamsort_cmd,
                                    stderr=open("/dev/null", 'w'))
    bamsort_proc.communicate()

    bamindex_cmd = ["samtools",
                    "index",
                    "hla_input.bam"]
    bamindex_proc = subprocess.Popen(bamindex_cmd,
                                     stderr=open("/dev/null", 'w'))
    bamindex_proc.communicate()
    """

    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Sort SAM alignments from HISAT2")
    parser.add_argument("-i", "--input",
                        dest="input_fname",
                        type=str,
                        default="-",
                        help="Input SAM file name (default: -, which is stdin)")
    parser.add_argument("--temp-fname",
                        dest="temp_fname",
                        type=str,
                        default="test_hisat2",
                        help="Temporary file name header")
    parser.add_argument("-p/--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads (default: 1")
    parser.add_argument("--interval-size",
                        dest="interval_size",
                        type=int,
                        default=30000000,
                        help="Interval size (default: 30000000, 30Mbps")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    sort_alignments(args.input_fname,
                    args.temp_fname,
                    args.threads,
                    args.interval_size,
                    args.verbose)
