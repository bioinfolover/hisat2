#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
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


import sys, os, subprocess, re
import inspect, random
import math
from argparse import ArgumentParser, FileType


def detectMisalignment(HLA_allele, seq, read, pos, seq_stats):
    seq_length = len(sequence)
    position = pos
    
    condition = True
    count = 0
    while(position < seq_length and count < len(read)):
        if HLA_allele[position] != read[count]:
            condition = False
            
            
            seq_stats[position][read[position]] += 1
        else:
            
            seq_stats[position]['I'] += 1
        
        count += 1
        position += 1
        
    return False


def checkRead(read, position, sequence, forward ):
    mis_positions = {}
    sub_sequence = sequence[position:position + len(read)]
    
    


"""
"""
def simulate_reads(HLAs,
                   test_HLA_list,
                   simulate_interval,
                   error_rate = -1):
                   
    HLA_reads_1, HLA_reads_2 = [], []
    for test_HLA_names in test_HLA_list:
        gene = test_HLA_names[0].split('*')[0]
        # ref_allele = refHLAs[gene]
        # ref_seq = HLAs[gene][ref_allele]

        # Simulate reads from two HLA alleles
        def simulate_reads_impl(seq, simulate_interval = 1, frag_len = 250, read_len = 100):
            comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            reads_1, reads_2 = [], []
            for i in range(0, len(seq) - frag_len + 1, simulate_interval):
                reads_1.append(seq[i:i+read_len])
                tmp_read_2 = reversed(seq[i+frag_len-read_len:i+frag_len])
                read_2 = ""
                for s in tmp_read_2:
                    if s in comp_table:
                        read_2 += comp_table[s]
                    else:
                        read_2 += s
                reads_2.append(read_2)
            return reads_1, reads_2
        # Simulate error in reads
        def simulate_error(reads, error_rate):
            newreads = []
            n_list = ['A','T','C','G']
            for read in reads:
                newread = ""
                for c in read:
                    if (error_rate >= random.uniform(0,1)):
                        i = n_list.index(c)
                        r = range(0,i) + range(i+1,4)
                        r = random.choice(r)
                        c = n_list[r]
                    newread += c
                newreads.append(newread)
            return newreads

        for test_HLA_name in test_HLA_names:
            HLA_seq = HLAs[gene][test_HLA_name]
            tmp_reads_1, tmp_reads_2 = simulate_reads_impl(HLA_seq, simulate_interval)
            if error_rate <= 0:
                HLA_reads_1 += tmp_reads_1
                HLA_reads_2 += tmp_reads_2
            else:
                HLA_reads_1 += simulate_error(tmp_reads_1, error_rate)
                HLA_reads_2 += simulate_error(tmp_reads_2, error_rate)

    # Write reads into a fasta read file
    def write_reads(reads, idx, label = ""):
        read_file = open('hla_input_%d.fa' % idx, 'w')
        for read_i in range(len(reads)):
            print >> read_file, ">%d" % (read_i + 1) + label
            print >> read_file, reads[read_i]
        read_file.close()
    write_reads(HLA_reads_1, 1, 'a')
    write_reads(HLA_reads_2, 2, 'b')


"""
Align reads, and sort the alignments into a BAM file
"""
def align_reads(ex_path,
                aligner,
                index_type,
                read_fname,
                fastq,
                threads,
                verbose):
    if aligner == "hisat2":
        hisat2 = os.path.join(ex_path, "hisat2")
        aligner_cmd = [hisat2,
                       "--no-unal",
                       "--mm"]
        if index_type == "linear":
            aligner_cmd += ["-k", "10"]
        aligner_cmd += ["-x", "hla.%s" % index_type]
    elif aligner == "bowtie2":
        aligner_cmd = [aligner,
                       "--no-unal",
                       "-k", "10",
                       "-x", "hla"]
    else:
        assert False
    assert len(read_fname) in [1,2]
    aligner_cmd += ["-p", str(threads)]
    if not fastq:
        aligner_cmd += ["-f"]
    if len(read_fname) == 1:
        aligner_cmd += ["-U", read_fname[0]]
    else:
        aligner_cmd += ["-1", "%s" % read_fname[0],
                        "-2", "%s" % read_fname[1]]

    if verbose:
        print >> sys.stderr, ' '.join(aligner_cmd)
    align_proc = subprocess.Popen(aligner_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=open("/dev/null", 'w'))

    sambam_cmd = ["samtools",
                  "view",
                  "-bS",
                  "-"]
    sambam_proc = subprocess.Popen(sambam_cmd,
                                   stdin=align_proc.stdout,
                                   stdout=open("hla_input_unsorted.bam", 'w'),
                                   stderr=open("/dev/null", 'w'))
    sambam_proc.communicate()
    if index_type == "graph":
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

        os.system("rm hla_input_unsorted.bam")            
    else:
        os.system("mv hla_input_unsorted.bam hla_input.bam")


"""
""" 
def normalize(prob):
    total = sum(prob.values())
    for allele, mass in prob.items():
        prob[allele] = mass / total

        
"""
"""
def prob_diff(prob1, prob2):
    diff = 0.0
    for allele in prob1.keys():
        if allele in prob2:
            diff += abs(prob1[allele] - prob2[allele])
        else:
            diff += prob1[allele]
    return diff


"""
"""
def HLA_prob_cmp(a, b):
    if a[1] != b[1]:
        if a[1] < b[1]:
            return 1
        else:
            return -1
    assert a[0] != b[0]
    if a[0] < b[0]:
        return -1
    else:
        return 1


"""
"""
def single_abundance(HLA_cmpt,
                     HLA_length):
    def normalize2(prob, length):
        total = 0
        for allele, mass in prob.items():
            assert allele in length
            total += (mass / length[allele])
        for allele, mass in prob.items():
            assert allele in length
            prob[allele] = mass / length[allele] / total

    HLA_prob, HLA_prob_next = {}, {}
    for cmpt, count in HLA_cmpt.items():
        alleles = cmpt.split('-')
        for allele in alleles:
            if allele not in HLA_prob:
                HLA_prob[allele] = 0.0
            HLA_prob[allele] += (float(count) / len(alleles))

    # normalize2(HLA_prob, HLA_length)
    normalize(HLA_prob)
    def next_prob(HLA_cmpt, HLA_prob, HLA_length):
        HLA_prob_next = {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            alleles_prob = 0.0
            for allele in alleles:
                assert allele in HLA_prob
                alleles_prob += HLA_prob[allele]
            for allele in alleles:
                if allele not in HLA_prob_next:
                    HLA_prob_next[allele] = 0.0
                HLA_prob_next[allele] += (float(count) * HLA_prob[allele] / alleles_prob)
        # normalize2(HLA_prob_next, HLA_length)
        normalize(HLA_prob_next)
        return HLA_prob_next

    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob, HLA_length)
        diff = prob_diff(HLA_prob, HLA_prob_next)
        HLA_prob = HLA_prob_next
        iter += 1
    for allele, prob in HLA_prob.items():
        allele_len = HLA_length[allele]
        HLA_prob[allele] /= float(allele_len)
    normalize(HLA_prob)
    HLA_prob = [[allele, prob] for allele, prob in HLA_prob.items()]
    HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
    return HLA_prob

    
"""
"""
def joint_abundance(HLA_cmpt,
                    HLA_length):
    allele_names = set()
    for cmpt in HLA_cmpt.keys():
        allele_names |= set(cmpt.split('-'))
    
    HLA_prob, HLA_prob_next = {}, {}
    for cmpt, count in HLA_cmpt.items():
        alleles = cmpt.split('-')
        for allele1 in alleles:
            for allele2 in allele_names:
                if allele1 < allele2:
                    allele_pair = "%s-%s" % (allele1, allele2)
                else:
                    allele_pair = "%s-%s" % (allele2, allele1)
                if not allele_pair in HLA_prob:
                    HLA_prob[allele_pair] = 0.0
                HLA_prob[allele_pair] += (float(count) / len(alleles))

    if len(HLA_prob) <= 0:
        return HLA_prob

    # Choose top allele pairs
    def choose_top_alleles(HLA_prob):
        HLA_prob_list = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
        HLA_prob_list = sorted(HLA_prob_list, cmp=HLA_prob_cmp)
        HLA_prob = {}
        best_prob = HLA_prob_list[0][1]
        for i in range(len(HLA_prob_list)):
            allele_pair, prob = HLA_prob_list[i]
            if prob * 2 <= best_prob:
                break                        
            HLA_prob[allele_pair] = prob
        normalize(HLA_prob)
        return HLA_prob
    HLA_prob = choose_top_alleles(HLA_prob)

    def next_prob(HLA_cmpt, HLA_prob):
        HLA_prob_next = {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            prob = 0.0
            for allele in alleles:
                for allele_pair in HLA_prob.keys():
                    if allele in allele_pair:
                        prob += HLA_prob[allele_pair]
            for allele in alleles:
                for allele_pair in HLA_prob.keys():
                    if not allele in allele_pair:
                        continue
                    if allele_pair not in HLA_prob_next:
                        HLA_prob_next[allele_pair] = 0.0
                    HLA_prob_next[allele_pair] += (float(count) * HLA_prob[allele_pair] / prob)
        normalize(HLA_prob_next)
        return HLA_prob_next

    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob)
        diff = prob_diff(HLA_prob, HLA_prob_next)
        HLA_prob = HLA_prob_next
        HLA_prob = choose_top_alleles(HLA_prob)
        iter += 1

    HLA_prob = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
    HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
    return HLA_prob


"""
"""
def HLA_typing(ex_path,
               simulation,
               reference_type,
               hla_list,
               partial,
               refHLAs,
               HLAs,
               HLA_names,
               HLA_lengths,
               refHLA_loci,
               Vars,
               Var_list,
               Links,
               exclude_allele_list,
               aligners,
               num_mismatch,
               fastq,
               read_fname,
               alignment_fname,
               threads,
               enable_coverage,
               best_alleles,
               detect_allele,
               HLA_fnames,
               verbose):

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    def lower_bound(Var_list, pos):
        low, high = 0, len(Var_list)
        while low < high:
            m = (low + high) / 2
            m_pos = Var_list[m][0]
            if m_pos < pos:
                low = m + 1
            elif m_pos > pos:
                high = m
            else:
                assert m_pos == pos
                while m > 0:
                    if Var_list[m-1][0] < pos:
                        break
                    m -= 1
                return m
        return low        
            
    if simulation:
        test_passed = {}
    for aligner, index_type in aligners:
        if index_type == "graph":
            print >> sys.stderr, "\n\t\t%s %s on %s" % (aligner, index_type, reference_type)
        else:
            print >> sys.stderr, "\n\t\t%s %s" % (aligner, index_type)

        if alignment_fname == "":
            # Align reads, and sort the alignments into a BAM file
            align_reads(ex_path,
                        aligner,
                        index_type,
                        read_fname,
                        fastq,
                        threads,
                        verbose)
            
        for test_HLA_names in hla_list:
            if simulation:
                gene = test_HLA_names[0].split('*')[0]
            else:
                gene = test_HLA_names
            ref_allele = refHLAs[gene]
            ref_seq = HLAs[gene][ref_allele]
            ref_exons = refHLA_loci[gene][-1]

            # Read alignments
            alignview_cmd = ["samtools",
                             "view"]
            if alignment_fname == "":
                alignview_cmd += ["hla_input.bam"]
            else:
                if not os.path.exists(alignment_fname + ".bai"):
                    os.system("samtools index %s" % alignment_fname)
                alignview_cmd += [alignment_fname]
            base_locus = 0
            if index_type == "graph":
                if reference_type == "gene":
                    alignview_cmd += ["%s" % ref_allele]
                else:
                    assert reference_type in ["chromosome", "genome"]
                    _, chr, left, right, _ = refHLA_loci[gene]
                    base_locus = left
                    alignview_cmd += ["%s:%d-%d" % (chr, left + 1, right + 1)]

                bamview_proc = subprocess.Popen(alignview_cmd,
                                                stdout=subprocess.PIPE,
                                                stderr=open("/dev/null", 'w'))

                sort_read_cmd = ["sort", "-k", "1", "-n"]
                alignview_proc = subprocess.Popen(sort_read_cmd,
                                                  stdin=bamview_proc.stdout,
                                                  stdout=subprocess.PIPE,
                                                  stderr=open("/dev/null", 'w'))
            else:
                alignview_proc = subprocess.Popen(alignview_cmd,
                                             stdout=subprocess.PIPE,
                                             stderr=open("/dev/null", 'w'))
            
            
            #print Vars['A']
            
            # Count alleles
            HLA_counts, HLA_cmpt = {}, {}
            coverage = [0 for i in range(len(ref_seq) + 1)]
            num_reads, total_read_len = 0, 0
            prev_read_id = None
            prev_exon = False
            if index_type == "graph":
                # Cigar regular expression
                cigar_re = re.compile('\d+\w')
                for line in alignview_proc.stdout:
                    cols = line.strip().split()
                    read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
                    read_seq, qual = cols[9], cols[10]
                    num_reads += 1
                    total_read_len += len(read_seq)
                    flag, pos = int(flag), int(pos)
                    pos -= (base_locus + 1)
                    if pos < 0:
                        continue

                    if flag & 0x4 != 0:
                        continue

                    NM, Zs, MD = "", "", ""
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("Zs"):
                            Zs = col[5:]
                        elif col.startswith("MD"):
                            MD = col[5:]
                        elif col.startswith("NM"):
                            NM = int(col[5:])

                    if NM > num_mismatch:
                        continue

                    # daehwan - for debugging purposes
                    debug = False
                    if read_id in ["2339"] and False:
                        debug = True
                        print "read_id: %s)" % read_id, pos, cigar_str, "NM:", NM, MD, Zs
                        print "            ", read_seq

                    vars = []
                    if Zs:
                        vars = Zs.split(',')

                    assert MD != ""
                    MD_str_pos, MD_len = 0, 0
                    read_pos, left_pos = 0, pos
                    right_pos = left_pos
                    cigars = cigar_re.findall(cigar_str)
                    cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
                    cmp_list = []
                    for i in range(len(cigars)):
                        cigar_op, length = cigars[i]
                        if cigar_op == 'M':
                            # Update coverage
                            if enable_coverage:
                                if right_pos + length < len(coverage):
                                    coverage[right_pos] += 1
                                    coverage[right_pos + length] -= 1
                                elif right_pos < len(coverage):
                                    coverage[right_pos] += 1
                                    coverage[-1] -= 1

                            first = True
                            MD_len_used = 0
                            while True:
                                if not first or MD_len == 0:
                                    if MD[MD_str_pos].isdigit():
                                        num = int(MD[MD_str_pos])
                                        MD_str_pos += 1
                                        while MD_str_pos < len(MD):
                                            if MD[MD_str_pos].isdigit():
                                                num = num * 10 + int(MD[MD_str_pos])
                                                MD_str_pos += 1
                                            else:
                                                break
                                        MD_len += num
                                # Insertion or full match followed
                                if MD_len >= length:
                                    MD_len -= length
                                    cmp_list.append(["match", right_pos + MD_len_used, length - MD_len_used])
                                    break
                                first = False
                                read_base = read_seq[read_pos + MD_len]
                                MD_ref_base = MD[MD_str_pos]
                                MD_str_pos += 1
                                assert MD_ref_base in "ACGT"
                                cmp_list.append(["match", right_pos + MD_len_used, MD_len - MD_len_used])
                                cmp_list.append(["mismatch", right_pos + MD_len, 1])
                                MD_len_used = MD_len + 1
                                MD_len += 1
                                # Full match
                                if MD_len == length:
                                    MD_len = 0
                                    break
                        elif cigar_op == 'I':
                            cmp_list.append(["insertion", right_pos, length])
                        elif cigar_op == 'D':
                            if MD[MD_str_pos] == '0':
                                MD_str_pos += 1
                            assert MD[MD_str_pos] == '^'
                            MD_str_pos += 1
                            while MD_str_pos < len(MD):
                                if not MD[MD_str_pos] in "ACGT":
                                    break
                                MD_str_pos += 1
                            cmp_list.append(["deletion", right_pos, length])
                        elif cigar_op == 'S':
                            cmp_list.append(["soft", right_pos, length])
                        else:                    
                            assert cigar_op == 'N'
                            cmp_list.append(["intron", right_pos, length])

                        if cigar_op in "MND":
                            right_pos += length

                        if cigar_op in "MIS":
                            read_pos += length

                    exon = False
                    for exon in ref_exons:
                        exon_left, exon_right = exon
                        if right_pos <= exon_left or pos > exon_right:
                            continue
                        else:
                            exon = True
                            break

                    if right_pos > len(ref_seq):
                        continue

                    def add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, exon = True):
                        max_count = max(HLA_count_per_read.values())
                        cur_cmpt = set()
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            if allele in exclude_allele_list:
                                continue                                
                            cur_cmpt.add(allele)                    
                            if not allele in HLA_counts:
                                HLA_counts[allele] = 1
                            else:
                                HLA_counts[allele] += 1

                        if len(cur_cmpt) == 0:
                            return

                        # daehwan - for debugging purposes                            
                        alleles = ["", ""]
                        # alleles = ["B*40:304", "B*40:02:01"]
                        allele1_found, allele2_found = False, False
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            if allele == alleles[0]:
                                allele1_found = True
                            elif allele == alleles[1]:
                                allele2_found = True
                        if allele1_found != allele2_found:
                            print alleles[0], HLA_count_per_read[alleles[0]]
                            print alleles[1], HLA_count_per_read[alleles[1]]
                            if allele1_found:
                                print ("%s\tread_id %s - %d vs. %d]" % (alleles[0], prev_read_id, max_count, HLA_count_per_read[alleles[1]]))
                            else:
                                print ("%s\tread_id %s - %d vs. %d]" % (alleles[1], prev_read_id, max_count, HLA_count_per_read[alleles[0]]))
                            print read_seq

                        cur_cmpt = sorted(list(cur_cmpt))
                        cur_cmpt = '-'.join(cur_cmpt)
                        add = 1
                        if partial and not exon:
                            add *= 0.2
                        if not cur_cmpt in HLA_cmpt:
                            HLA_cmpt[cur_cmpt] = add
                        else:
                            HLA_cmpt[cur_cmpt] += add

                    if read_id != prev_read_id:
                        if prev_read_id != None:
                            add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, prev_exon)

                        HLA_count_per_read = {}
                        for HLA_name in HLA_names[gene]:
                            if HLA_name.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[HLA_name] = 0

                    def add_count(var_id, add):
                        assert var_id in Links
                        alleles = Links[var_id]
                        for allele in alleles:
                            if allele.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[allele] += add
                            # daehwan - for debugging purposes
                            if debug:
                                if allele in ["DQA1*05:05:01:01", "DQA1*05:05:01:02"]:
                                    print allele, add, var_id

                    # Decide which allele(s) a read most likely came from
                    # also sanity check - read length, cigar string, and MD string
                    for var_id, data in Vars[gene].items():
                        var_type, var_pos, var_data = data
                        if var_type != "deletion":
                            continue
                        if left_pos >= var_pos and right_pos <= var_pos + int(var_data):
                            add_count(var_id, -1)                            
                    ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
                    cigar_match_len, MD_match_len = 0, 0            
                    for cmp in cmp_list:
                        type = cmp[0]
                        length = cmp[2]
                        if type == "match":
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos + length <= var_pos:
                                    break
                                if ref_pos <= var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":
                                        if ref_pos < var_pos and ref_pos + length > var_pos + len(var_data):
                                            add_count(var_id, -1)
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, Links[var_id]
                                    elif var_type == "deletion":
                                        del_len = int(var_data)
                                        if ref_pos < var_pos and ref_pos + length > var_pos + del_len:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, Links[var_id], -1, Vars[gene][var_id]
                                            # Check if this might be one of the two tandem repeats (the same left coordinate)
                                            cmp_left, cmp_right = cmp[1], cmp[1] + cmp[2]
                                            test1_seq1 = ref_seq[cmp_left:cmp_right]
                                            test1_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos + del_len:cmp_right + del_len]
                                            # Check if this happens due to small repeats (the same right coordinate - e.g. 19 times of TTTC in DQA1*05:05:01:02)
                                            cmp_left -= read_pos
                                            cmp_right += (len(read_seq) - read_pos - cmp[2])
                                            test2_seq1 = ref_seq[cmp_left+int(var_data):cmp_right]
                                            test2_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos+int(var_data):cmp_right]
                                            if test1_seq1 != test1_seq2 and test2_seq1 != test2_seq2:
                                                add_count(var_id, -1)
                                    else:
                                        if debug:
                                            print cmp, var_id, Links[var_id], -1
                                        add_count(var_id, -1)
                                var_idx += 1

                            read_pos += length
                            ref_pos += length
                            cigar_match_len += length
                            MD_match_len += length
                        elif type == "mismatch":
                            read_base = read_seq[read_pos]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "single":
                                        if var_data == read_base:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, 1, var_data, read_base, Links[var_id]

                                            # daehwan - for debugging purposes
                                            if False:
                                                read_qual = ord(qual[read_pos])
                                                add_count(var_id, (read_qual - 60) / 60.0)
                                            else:
                                                add_count(var_id, 1)
                                        # daehwan - check out if this routine is appropriate
                                        # else:
                                        #    add_count(var_id, -1)
                                var_idx += 1

                            cmp_MD += ("%d%s" % (MD_match_len, ref_seq[ref_pos]))
                            MD_match_len = 0
                            cigar_match_len += 1
                            read_pos += 1
                            ref_pos += 1
                        elif type == "insertion":
                            ins_seq = read_seq[read_pos:read_pos+length]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            # daehwan - for debugging purposes
                            if debug:
                                print left_pos, cigar_str, MD, vars
                                print ref_pos, ins_seq, Var_list[gene][var_idx], Vars[gene][Var_list[gene][var_idx][1]]
                                # sys.exit(1)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":                                
                                        if var_data == ins_seq:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, 1, Links[var_id]
                                            add_count(var_id, 1)
                                var_idx += 1

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dI" % length)
                        elif type == "deletion":
                            del_len = length
                            # Deletions can be shifted bidirectionally
                            temp_ref_pos = ref_pos
                            while temp_ref_pos > 0:
                                last_bp = ref_seq[temp_ref_pos + del_len - 1]
                                prev_bp = ref_seq[temp_ref_pos - 1]
                                if last_bp != prev_bp:
                                    break
                                temp_ref_pos -= 1
                            var_idx = lower_bound(Var_list[gene], temp_ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if temp_ref_pos < var_pos:
                                    first_bp = ref_seq[temp_ref_pos]
                                    next_bp = ref_seq[temp_ref_pos + del_len]
                                    if first_bp == next_bp:
                                        temp_ref_pos += 1
                                        continue
                                    else:
                                        break
                                if temp_ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "deletion":
                                        var_len = int(var_data)
                                        if var_len == length:
                                            if debug:
                                                print cmp, var_id, 1, Links[var_id]
                                                print ref_seq[var_pos - 10:var_pos], ref_seq[var_pos:var_pos+int(var_data)], ref_seq[var_pos+int(var_data):var_pos+int(var_data)+10]
                                            add_count(var_id, 1)
                                var_idx += 1

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_MD += ("%d" % MD_match_len)
                            MD_match_len = 0
                            cmp_cigar_str += ("%dD" % length)
                            cmp_MD += ("^%s" % ref_seq[ref_pos:ref_pos+length])
                            ref_pos += length
                        elif type == "soft":
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dS" % length)
                        else:
                            assert type == "intron"
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_cigar_str += ("%dN" % length)
                            ref_pos += length                    
                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                    cmp_MD += ("%d" % MD_match_len)
                    if read_pos != len(read_seq) or \
                            cmp_cigar_str != cigar_str or \
                            cmp_MD != MD:
                        print >> sys.stderr, "Error:", cigar_str, MD
                        print >> sys.stderr, "\tcomputed:", cmp_cigar_str, cmp_MD
                        print >> sys.stderr, "\tcmp list:", cmp_list
                        assert False            

                    prev_read_id = read_id
                    prev_exon = exon

                if num_reads <= 0:
                    continue

                if prev_read_id != None:
                    add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read)

                # Coverage
                # it is not used by the default
                if enable_coverage:
                    assert num_reads > 0
                    read_len = int(total_read_len / float(num_reads))
                    coverage_sum = 0
                    for i in range(len(coverage)):
                        if i > 0:
                            coverage[i] += coverage[i-1]
                        coverage_sum += coverage[i]
                    coverage_avg = coverage_sum / float(len(coverage))
                    assert len(ref_seq) < len(coverage)
                    for i in range(len(ref_seq)):
                        coverage_threshold = 1.0 * coverage_avg
                        if i < read_len:
                            coverage_threshold *= ((i+1) / float(read_len))
                        elif i + read_len > len(ref_seq):
                            coverage_threshold *= ((len(ref_seq) - i) / float(read_len))
                        if coverage[i] >= coverage_threshold:
                            continue
                        pseudo_num_reads = (coverage_threshold - coverage[i]) / read_len
                        var_idx = lower_bound(Var_list[gene], i + 1)
                        if var_idx >= len(Var_list[gene]):
                            var_idx = len(Var_list[gene]) - 1
                        cur_cmpt = set()
                        while var_idx >= 0:
                            var_pos, var_id = Var_list[gene][var_idx]
                            var_type, _, var_data = Vars[gene][var_id]
                            if var_type == "deletion":
                                del_len = int(var_data)
                                if i < var_pos:
                                    break
                                if i + read_len < var_pos + int(var_data):
                                    assert var_id in Links
                                    cur_cmpt = cur_cmpt.union(set(Links[var_id]))
                            var_idx -= 1
                        if cur_cmpt:
                            cur_cmpt = '-'.join(list(cur_cmpt))
                            if not cur_cmpt in HLA_cmpt:
                                HLA_cmpt[cur_cmpt] = 0
                            HLA_cmpt[cur_cmpt] += pseudo_num_reads
            else:
                assert index_type == "linear"
                def add_alleles(alleles):
                    if not allele in HLA_counts:
                        HLA_counts[allele] = 1
                    else:
                        HLA_counts[allele] += 1

                    cur_cmpt = sorted(list(alleles))
                    cur_cmpt = '-'.join(cur_cmpt)
                    if not cur_cmpt in HLA_cmpt:
                        HLA_cmpt[cur_cmpt] = 1
                    else:
                        HLA_cmpt[cur_cmpt] += 1

                prev_read_id, prev_AS = None, None
                alleles = set()
                for line in alignview_proc.stdout:
                    cols = line[:-1].split()
                    read_id, flag, allele = cols[:3]
                    flag = int(flag)
                    if flag & 0x4 != 0:
                        continue
                    if not allele.startswith(gene):
                        continue
                    if allele.find("BACKBONE") != -1:
                        continue

                    AS = None
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("AS"):
                            AS = int(col[5:])
                    assert AS != None
                    if read_id != prev_read_id:
                        if alleles:
                            if aligner == "hisat2" or \
                                    (aligner == "bowtie2" and len(alleles) < 10):
                                add_alleles(alleles)
                            alleles = set()
                        prev_AS = None
                    if prev_AS != None and AS < prev_AS:
                        continue
                    prev_read_id = read_id
                    prev_AS = AS
                    alleles.add(allele)

                if alleles:
                    add_alleles(alleles)

            HLA_counts = [[allele, count] for allele, count in HLA_counts.items()]
            def HLA_count_cmp(a, b):
                if a[1] != b[1]:
                    return b[1] - a[1]
                assert a[0] != b[0]
                if a[0] < b[0]:
                    return -1
                else:
                    return 1
            HLA_counts = sorted(HLA_counts, cmp=HLA_count_cmp)
            for count_i in range(len(HLA_counts)):
                count = HLA_counts[count_i]
                if simulation:
                    found = False
                    for test_HLA_name in test_HLA_names:
                        if count[0] == test_HLA_name:
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (count: %d)" % (count_i + 1, test_HLA_name, count[1])
                            found = True
                            """
                            if count_i > 0 and HLA_counts[0][1] > count[1]:
                                print >> sys.stderr, "Warning: %s ranked first (count: %d)" % (HLA_counts[0][0], HLA_counts[0][1])
                                assert False
                            else:
                                test_passed += 1
                            """
                    if count_i < 5 and not found:
                        print >> sys.stderr, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])
                else:
                    print >> sys.stderr, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])
                    if count_i >= 9:
                        break
            print >> sys.stderr

            HLA_prob = single_abundance(HLA_cmpt, HLA_lengths[gene])
            
            # Begin read alignment check for novel alleles    
            if index_type == "graph" and detect_allele:
                #num_alleles = 1
                alleles_found = []
                
                try:
                    os.mkdir("./Genotyping".format(i))
                except:
                    pass  
                
                #Identify closest matching alleles from results.    
                if len(test_HLA_names) == 2:
                    HLA_prob_temp = joint_abundance(HLA_cmpt, HLA_lengths[gene])
                    alleles_found = HLA_prob_temp[0][0].split('-')
                else:
                    alleles_found.append(HLA_prob[0][0])
                allele_file = {}
                allele_alignment_data = {}
                
                #Setup alignment results for both alleles
                for allele_det in alleles_found:
                    #print HLAs
                    gene_name = allele_det.split('*')[0]
                    allele_sequence = HLAs[gene_name][allele_det]
                    allele_alignment_data[allele_det] = []
                    for i in range(0,HLA_lengths[gene_name][allele_det]):
                        allele_alignment_data[allele_det].append({})
                    
                    first_time = False
                    try:
                        os.mkdir("./Genotyping/{0}".format(allele_det))
                        first_time = True
                    except:
                        pass  
                    
                    #print "./Genotyping/Allele{0}/{1}.fa".format(i,allele_det)
                    
                    seq_dest = "./Genotyping/{0}/{0}.fa".format(allele_det)
                    
                    #Perform read alignment to closest alleles. Store results in Genotyping folder
                    HLA_hisat2_graph_index_fnames = ["./Genotyping/{0}/hla.graph.{1}.ht2".format(allele_det,(i+1)) for i in range(8)]
                    HLA_hisat2_graph_index_fnames.append(seq_dest)
                    if not check_files(HLA_hisat2_graph_index_fnames):
                        seq_file = open(seq_dest,"w")
                        seq_file.write(">{0}\n".format(allele_det))
                        
                        for s in range(0, len(allele_sequence), 60):
                            seq_file.write(allele_sequence[s:s+60] + "\n")
                        seq_file.close()
                        
                        seq_file = open("./Genotyping/{0}/Sequence".format(allele_det),"w")
                        seq_file.write(allele_sequence)
                        seq_file.close()
                        
                        hisat2_build = os.path.join(ex_path, "hisat2-build")
                        build_cmd = [hisat2_build,
                                     "-p", str(threads),
                                     "--snp", HLA_fnames[3],
                                     "--haplotype", HLA_fnames[4] ,
                                     seq_dest,
                                     "./Genotyping/{0}/hla.graph".format(allele_det)]
                        if verbose:
                            print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
                        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                        proc.communicate()        
                                        
                    hisat2_path = os.path.join(ex_path, "hisat2")
                    aligner_cmd = [hisat2_path,
                                    "--mm",
                                    "-x", "./Genotyping/{0}/hla.graph".format(allele_det),
                                    "-p","{0}".format(1),
                                    "-f",
                                    "-1", "hla_input_1.fa",
                                    "-2", "hla_input_2.fa",
                                    "-S", "./Genotyping/{0}/Output".format(allele_det)]
                    if verbose:
                        print >> sys.stderr, ' '.join(aligner_cmd)
                    align_proc = subprocess.Popen(aligner_cmd,
                                                  stdout=subprocess.PIPE,
                                                  stderr=open("/dev/null", 'w'))
                                                  
                                                  
                    align_proc.communicate()
                    #print "Alignment Complete"                   
                    
                    allele_file[allele_det] = open("./Genotyping/{0}/Output".format(allele_det),"r")
                    allele_file[allele_det].readline()
                    allele_file[allele_det].readline()
                    allele_file[allele_det].readline()
            
            
                read_count = 0
                readsleft = True
            
                #Determine whether a given read matches to one or both alleles.
                #May need to check 'M' sites for mismatches.
                def best_read(read_data, allele_list):
                    best = [allele_list[0]]
                    b_matches = 0
                    
                    if len(allele_list) == 1:
                        return best
                    cigar_re = re.compile('\d+\w')
                    
                    setup = True
                    for allele in allele_list:
                        if(len(read_data[allele]) < 6):
                            continue
                        
                        cigar_data = cigar_re.findall(read_data[allele][5])
                        c_matches = 0
                        for x in cigar_data:
                            if 'M' in x:
                                c_matches += int(x.split('M')[0])
		                    
                            if 'D' in x:
                                c_matches -= int(x.split('D')[0])
                        #print allele
                        #print "Score: {0}".format(c_matches)
                        if setup:
                            setup = False
                            b_matches = c_matches
                        elif (c_matches == b_matches):
                            best.append(allele)
                        elif (c_matches > b_matches):
                            b_matches = c_matches
                            best = [allele]
                            
                    return best
                #Update reference with read data.    
                def updateRef(reference_allele, pos, seq):
                    if(pos >= len(reference_allele)):
                        return False
                    if reference_allele[pos] == None:
                        reference_allele[pos] = {}
                    
                    if(not seq in reference_allele[pos]):
                        reference_allele[pos][seq] = 1
                    else:
                        reference_allele[pos][seq] += 1
                    return True
                    
                    
                #Extract consensus sequence from read alignment data 
                def extract_sequence(read_seq, cigars, pos, reference_allele):
                    #print pos
                    current_pos = pos - 1 
                    count = 0
                    #while( current < reference_allele && count < read_length):
                    for cigar in cigars:
                        num, score = filter(None,re.split('(\d+)',cigar))
                        num = int(num)
                        if score == 'M' or 'S' or 'H' or '=' or 'X':
                            for i in range(0,num):
                                if(count + i < len(read_seq)):
                                    updateRef(reference_allele, current_pos + i, read_seq[count + i])
                            current_pos += num
                            count += num
                        elif score == 'D':
                            for i in range(0,num):
                                updateRef(reference_allele, current_pos + i, 'D')
                            current_pos += num
                            count += num
                        elif score == 'I':
                            updateRef(reference_allele, current_pos, read_seq[count:count+num])
                            count += num
                            current_pos += num
                        else:
                            count += num
                            current_pos += num
                
                #Scan through read alignments to closest alleles and output novel sequences (if any are found)
                while(readsleft):
                    allele_reads={}
                    for allele_det in alleles_found:
                        read = allele_file[allele_det].readline()
                        #print read
                        if(read == ''):
                            readsleft = False
                            break
                        allele_reads[allele_det] = read.split()
                        
                    if(readsleft):
                        alleles_sel = best_read(allele_reads,alleles_found)
                        for allele_sel in alleles_sel:
                            if(len(allele_reads[allele_sel]) >= 6):
                                cigar_data = cigar_re.findall(allele_reads[allele_sel][5])
                                extract_sequence(allele_reads[allele_sel][9],cigar_data,int(allele_reads[allele_sel][3]), allele_alignment_data[allele_sel])
                                
                            
                #Output novel allele sequence with variant report.
                def reportNewAllele(allele_data, allele_sequence):
                    i = 0
                    match = True
                    sequence = ""
                    for c in allele_sequence:
                        #if(i < 5):
                        #    i+=1
                        #    continue
                        
                        #if i >= len(allele_data):
                        #    break
                        if allele_data[i] == None or len(allele_data[i]) == 0:
                            print "Missing Coverage at Position: {0}".format(i+1)
                            match = False
                            i += 1
                            continue
                        #print allele_data[i]
                        variants = sorted(allele_data[i], key=lambda k: int(allele_data[i][k]), reverse=True)
                        sequence += variants[0]
                        if(variants[0]!= c):
                            if(variants[0] == 'D'):
                                print "Nucleotide Deletion at Position: {0}".format(i+1)
                            elif(len(variants[0]) > 1):
                                print "Insertion Sequence {0} at Position: {1}".format(variants[0],i+1)
                            else:
                                print "Nucleotide Change to {0} at Position: {1}".format(variants[0],i+1)
                                
                            match = False
                            
                        #match = match and (variants[0] == c)

                        i += 1
                    #Enable printing of novel sequence
                    print_output_sequence = True
                    if (not match and print_output_sequence):
                        #print allele_data
                        print "Assembled Sequence For New Allele:"
                        for s in range(0, len(sequence), 60):
                            print sequence[s:s+60]
                    return match
                    
                for allele_det in alleles_found:
                    gene_name = allele_det.split('*')[0]
                    allele_sequence = HLAs[gene_name][allele_det]
                    #print allele_alignment_data[allele_sel]
                    if(not reportNewAllele(allele_alignment_data[allele_det],allele_sequence)):
                        print "New Allele Detected. Highest match to existing allele: {0}".format(allele_det)
                    else:
                        print "Complete match to allele {0}".format(allele_det) 
                
                
            success = [False for i in range(len(test_HLA_names))]
            found_list = [False for i in range(len(test_HLA_names))]
            for prob_i in range(len(HLA_prob)):
                prob = HLA_prob[prob_i]
                found = False
                if simulation:
                    for name_i in range(len(test_HLA_names)):
                        test_HLA_name = test_HLA_names[name_i]
                        if prob[0] == test_HLA_name:
                            rank_i = prob_i
                            while rank_i > 0:
                                if prob == HLA_prob[rank_i - 1][1]:
                                    rank_i -= 1
                                else:
                                    break
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, test_HLA_name, prob[1] * 100.0)
                            if rank_i < len(success):
                                success[rank_i] = True
                            found_list[name_i] = True
                            found = True                        
                    if not False in found_list:
                        break
                if not found:
                    print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, prob[0], prob[1] * 100.0)
                    if best_alleles and prob_i < 2:
                        print >> sys.stdout, "SingleModel %s (abundance: %.2f%%)" % (prob[0], prob[1] * 100.0)
                if not simulation and prob_i >= 9:
                    break
            print >> sys.stderr

            if len(test_HLA_names) == 2 or not simulation:
                HLA_prob = joint_abundance(HLA_cmpt, HLA_lengths[gene])
                if len(HLA_prob) <= 0:
                    continue
                success = [False]
                for prob_i in range(len(HLA_prob)):
                    allele_pair, prob = HLA_prob[prob_i]
                    allele1, allele2 = allele_pair.split('-')
                    if best_alleles and prob_i < 1:
                        print >> sys.stdout, "PairModel %s (abundance: %.2f%%)" % (allele_pair, prob * 100.0)
                    if simulation:
                        if allele1 in test_HLA_names and allele2 in test_HLA_names:
                            rank_i = prob_i
                            while rank_i > 0:
                                if HLA_prob[rank_i-1][1] == prob:                                        
                                    rank_i -= 1
                                else:
                                    break
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, allele_pair, prob * 100.0)
                            if rank_i == 0:
                                success[0] = True
                            break
                    print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, allele_pair, prob * 100.0)
                    if not simulation and prob_i >= 9:
                        break
                print >> sys.stderr

                # Li's method
                """
                li_hla = os.path.join(ex_path, "li_hla/hla")
                if os.path.exists(li_hla):
                    li_hla_cmd = [li_hla,
                                  "hla",
                                  "hla_input.bam",
                                  "-b", "%s*BACKBONE" % gene]
                    li_hla_proc = subprocess.Popen(li_hla_cmd,
                                                   stdout=subprocess.PIPE,
                                                   stderr=open("/dev/null", 'w'))

                    # read in the result of Li's hla
                    for line in li_hla_proc.stdout:
                        allele1, allele2, score = line.strip().split()
                        score = float(score)
                        if simulation:
                            if allele1 in test_HLA_names and allele2 in test_HLA_names:
                                print >> sys.stderr, "\t\t\t*** 1 ranked %s-%s (score: %.2f)" % (allele1, allele2, score)
                                success[0] = True
                            else:
                                print >> sys.stderr, "\t\t\tLiModel fails"
                        if best_alleles:
                            print >> sys.stdout, "LiModel %s-%s (score: %.2f)" % (allele1, allele2, score)
                    li_hla_proc.communicate()
                """

            if simulation and not False in success:
                aligner_type = "%s %s" % (aligner, index_type)
                if not aligner_type in test_passed:
                    test_passed[aligner_type] = 1
                else:
                    test_passed[aligner_type] += 1

    #print HLA_counts
    if simulation:
        return test_passed


def read_HLA_alleles(fname, HLAs):
    for line in open(fname):
        if line.startswith(">"):
            HLA_name = line.strip().split()[0][1:]
            HLA_gene = HLA_name.split('*')[0]
            if not HLA_gene in HLAs:
                HLAs[HLA_gene] = {}
            if not HLA_name in HLAs[HLA_gene]:
                HLAs[HLA_gene][HLA_name] = ""
        else:
            HLAs[HLA_gene][HLA_name] += line.strip()
    return HLAs

"""
"""
def test_HLA_genotyping(base_fname,
                        reference_type,
                        hla_list,
                        partial,
                        aligners,
                        read_fname,
                        alignment_fname,
                        threads,
                        simulate_interval,
                        enable_coverage,
                        best_alleles,
                        exclude_allele_list,
                        default_allele_list,
                        num_mismatch,
                        verbose,
                        detect_allele,
                        error_rate,
                        daehwan_debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")

    simulation = (read_fname == [] and alignment_fname == "")

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if os.path.exists("hla.ref"):
        left = 0
        HLA_genes = set()
        BACKBONE = False
        for line in open("hla.ref"):
            HLA_name = line.strip().split()[0]
            if HLA_name.find("BACKBONE") != -1:
                BACKBONE = True
            HLA_gene = HLA_name.split('*')[0]
            HLA_genes.add(HLA_gene)
        delete_hla_files = False
        if reference_type == "gene":
            if not BACKBONE:
                delete_hla_files = True
        elif reference_type in ["chromosome", "genome"]:
            if BACKBONE:
                delete_hla_files = True
        else:
            assert False
        if not set(hla_list).issubset(HLA_genes):
            delete_hla_files = True
        if delete_hla_files:
            os.system("rm hla*")
    
    # Extract HLA variants, backbone sequence, and other sequeces  
    if len(base_fname) > 0:
        base_fname = "_" + base_fname
    base_fname = "hla" + base_fname
    
    HLA_fnames = [base_fname+"_backbone.fa",
                  base_fname+"_sequences.fa",
                  base_fname+".ref",
                  base_fname+".snp",
                  base_fname+".haplotype",
                  base_fname+".link",
                  base_fname+"_alleles_excluded.txt"]

    
    # Check if excluded alleles in current files match
    excluded_alleles_match = False
    if(os.path.exists(HLA_fnames[6])):
        afile = open(HLA_fnames[6],'r')
        afile.readline()
        lines = afile.read().split()
        excluded_alleles_match = (set(exclude_allele_list) == set(lines))
        afile.close()
    elif len(exclude_allele_list) == 0:
        excluded_alleles_match = True
        try:
            temp_name = HLA_fnames[6]
            HLA_fnames.remove(HLA_fnames[6])
            os.remove(temp_name)
        except OSError:
            pass
        
    if not excluded_alleles_match:
        print("Creating Allele Exclusion File. Alleles excluded: {0}\n".format(exclude_allele_list))
        afile = open(HLA_fnames[6],'w')
        afile.write("Alleles excluded:\n")
        afile.write("\n".join(exclude_allele_list))
        afile.close()
        
    #print HLA_fnames
    
    if (not check_files(HLA_fnames)) or (not excluded_alleles_match) :
        extract_hla_script = os.path.join(ex_path, "hisat2_extract_HLA_vars.py")
        extract_cmd = [extract_hla_script,
                       "--reference-type", reference_type,
                       "--hla-list", ','.join(hla_list)]

        if len(exclude_allele_list) > 0:
            #print exclude_allele_list
            extract_cmd += ["--exclude-allele-list", ",".join(exclude_allele_list)]

        if len(base_fname) > 3:
            extract_cmd += ["--base", base_fname]

        if partial:
            extract_cmd += ["--partial"]
        extract_cmd += ["--inter-gap", "30",
                        "--intra-gap", "50"]
        if verbose:
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)
            
    print "Base files built\n"

    # Build HISAT2 graph indexes based on the above information
    HLA_hisat2_graph_index_fnames = ["hla.graph.%d.ht2" % (i+1) for i in range(8)]
    if not check_files(HLA_hisat2_graph_index_fnames) or (not excluded_alleles_match):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     "-p", str(threads),
                     "--snp", HLA_fnames[3],
                     "--haplotype", HLA_fnames[4] ,
                     HLA_fnames[0],
                     "hla.graph"]
        if verbose:
            print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_hisat2_graph_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!  Perhaps, you may have forgotten to build hisat2 executables?"
            sys.exit(1)
    print "Step 1 Complete\n"
    # Build HISAT2 linear indexes based on the above information
    HLA_hisat2_linear_index_fnames = ["hla.linear.%d.ht2" % (i+1) for i in range(8)]
    if reference_type == "gene" and (not check_files(HLA_hisat2_linear_index_fnames) or (not excluded_alleles_match)):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     "%s,%s"%(HLA_fnames[0],HLA_fnames[1]),
                     "hla.linear"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_hisat2_linear_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)
            
    print "Step 2 Complete\n"
    # Build Bowtie2 indexes based on the above information
    HLA_bowtie2_index_fnames = ["hla.%d.bt2" % (i+1) for i in range(4)]
    HLA_bowtie2_index_fnames += ["hla.rev.%d.bt2" % (i+1) for i in range(2)]
    if reference_type == "gene" and (not check_files(HLA_bowtie2_index_fnames) or (not excluded_alleles_match)):
        build_cmd = ["bowtie2-build",
                     "%s,%s"%(HLA_fnames[0],HLA_fnames[1]),
                     "hla"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_bowtie2_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)

    print "Step 3 Complete\n"
    # Read partial alleles from hla.data (temporary)
    partial_alleles = set()
    for line in open("IMGTHLA/hla.dat"):
        if not line.startswith("DE"):
            continue
        allele_name = line.split()[1][4:-1]
        gene = allele_name.split('*')[0]
        if line.find("partial") != -1:
            partial_alleles.add(allele_name)

    if len(default_allele_list) != 0:
        #print os.getcwd()
        if not os.path.exists("./Default-HLA/hla_backbone.fa"):
            #current_path = os.getcwd()
            try:
                os.mkdir("./Default-HLA")
            except:
                pass
            #os.chdir(current_path + "/Default-HLA")
            
            extract_hla_script = os.path.join(ex_path, "hisat2_extract_HLA_vars.py")
            extract_cmd = [extract_hla_script,
                           "--reference-type", reference_type,
                           "--hla-list", ','.join(hla_list),
                           "--base", "./Default-HLA/hla"]

            if partial:
                extract_cmd += ["--partial"]
            extract_cmd += ["--inter-gap", "30",
                            "--intra-gap", "50"]
            if verbose:
                print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
            proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            
            if not os.path.exists("./Default-HLA/hla_backbone.fa"):
                print >> sys.stderr, "Error: extract_HLA_vars (Default) failed!"
                sys.exit(1)
    
    # Read HLA alleles (names and sequences)
    refHLAs, refHLA_loci = {}, {}
    for line in open("hla.ref"):
        HLA_name, chr, left, right, length, exon_str = line.strip().split()
        HLA_gene = HLA_name.split('*')[0]
        assert not HLA_gene in refHLAs
        refHLAs[HLA_gene] = HLA_name
        left, right = int(left), int(right)
        exons = []
        for exon in exon_str.split(','):
            exon_left, exon_right = exon.split('-')
            exons.append([int(exon_left), int(exon_right)])
        refHLA_loci[HLA_gene] = [HLA_name, chr, left, right, exons]
    HLAs = {}



    if reference_type == "gene":
        read_HLA_alleles(HLA_fnames[0], HLAs)
    read_HLA_alleles(HLA_fnames[1], HLAs)
    
    # HLA gene alleles
    HLA_names = {}
    for HLA_gene, data in HLAs.items():
        HLA_names[HLA_gene] = list(data.keys())

    # HLA gene allele lengths
    HLA_lengths = {}
    for HLA_gene, HLA_alleles in HLAs.items():
        HLA_lengths[HLA_gene] = {}
        for allele_name, seq in HLA_alleles.items():
            HLA_lengths[HLA_gene][allele_name] = len(seq)

    # Construct excluded alleles (Via default backbone data)
    custom_allele_check = False
    if len(default_allele_list) > 0:
        custom_allele_check = True
        HLAs_default = {}
        read_HLA_alleles("./Default-HLA/hla_backbone.fa",HLAs_default)
        read_HLA_alleles("./Default-HLA/hla_sequences.fa",HLAs_default)
        #HLA_lengths_default = {}
        
        
        for HLA_gene, HLA_alleles in HLAs_default.items():
            for allele_name, seq in HLA_alleles.items():
                if allele_name in default_allele_list:
                    HLA_lengths[HLA_gene][allele_name] = len(seq)
        
        #for allele_name, seq in HLAs_default.items():
         #   if allele_name in default_allele_list:
          #      HLA_lengths[allele_name] = len(seq)
            #if (allele_name in default_allele_list):
            #    HLA_lengths_default[allele_name] = len(seq)


    # Read HLA variants, and link information
    Vars, Var_list = {}, {}
    for line in open(HLA_fnames[3]):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        pos = int(pos)
        if reference_type != "gene":
            allele, dist = None, 0
            for tmp_gene, values in refHLA_loci.items():
                allele_name, chr, left, right, exons = values
                if allele == None:
                    allele = allele_name
                    dist = abs(pos - left)
                else:
                    if dist > abs(pos - left):
                        allele = allele_name
                        dist = abs(pos - left)
            
        gene = allele.split('*')[0]
        if not gene in Vars:
            Vars[gene] = {}
            assert not gene in Var_list
            Var_list[gene] = []
            
        assert not var_id in Vars[gene]
        left = 0
        if reference_type != "gene":
            _, _, left, _, _ = refHLA_loci[gene]
        Vars[gene][var_id] = [var_type, pos - left, data]
        Var_list[gene].append([pos - left, var_id])
        
    for gene, in_var_list in Var_list.items():
        Var_list[gene] = sorted(in_var_list)
        
    Links = {}
    for line in open(HLA_fnames[5]):
        var_id, alleles = line.strip().split('\t')
        alleles = alleles.split()
        assert not var_id in Links
        Links[var_id] = alleles

    # Scoring schemes from Sangtae Kim (Illumina)'s implementation
    # Currently not used.
    """
    max_qual_value = 100
    match_score, mismatch_score = [0] * max_qual_value, [0] * max_qual_value
    for qual in range(max_qual_value):
        error_rate = 0.1 ** (qual / 10.0)
        match_score[qual] = math.log(1.000000000001 - error_rate);
        mismatch_score[qual] = math.log(error_rate / 3.0);
    """
# Test HLA typing
    test_list = []
    if simulation:
        basic_test, pair_test = True, False
        if daehwan_debug:
            if "basic_test" in daehwan_debug:
                basic_test, pair_test = True, False
            else:
                basic_test, pair_test = False, True

        test_passed = {}
        test_list = []
        genes = list(set(hla_list) & set(HLA_names.keys()))
        if basic_test:
            if custom_allele_check:
                for allele in default_allele_list:
                    test_list.append([[allele]])
            else:
                for gene in genes:
                    HLA_gene_alleles = HLA_names[gene]
                    for HLA_name in HLA_gene_alleles:
                        if HLA_name.find("BACKBONE") != -1:
                            continue
                        test_list.append([[HLA_name]])
        if pair_test:
            test_size = 500
            allele_count = 2
            if custom_allele_check:
                if (default_allele_list) < allele_count:
                    print >> sys.stderr, "# of default alleles (%d) is at least %d" % (len(defeault_allele_list), allele_count)
                    sys.exit(1)
                    
                for test_i in range(1):
                    random.shuffle(default_allele_list)
                    test_pair = [default_allele_list[:allele_count]]
                    test_list.append(test_pair)
            else:
                for test_i in range(test_size):
                    test_pairs = []
                    for gene in genes:
                        HLA_gene_alleles = []

                        for allele in HLA_names[gene]:
                            if allele.find("BACKBONE") != -1:
                                continue
                            HLA_gene_alleles.append(allele)
                        nums = [i for i in range(len(HLA_gene_alleles))]
                        random.shuffle(nums)
                        test_pairs.append(sorted([HLA_gene_alleles[nums[i]] for i in range(allele_count)]))
                    test_list.append(test_pairs)

        for test_i in range(len(test_list)):
            if "test_id" in daehwan_debug:
                daehwan_test_ids = daehwan_debug["test_id"].split('-')
                if str(test_i + 1) not in daehwan_test_ids:
                    continue

            print >> sys.stderr, "Test %d" % (test_i + 1)
            test_HLA_list = test_list[test_i]
           
            # DK - for debugging purposes
            # test_HLA_list = [["A*11:50Q", "A*11:01:01:01", "A*01:01:01:01"]]
            for test_HLA_names in test_HLA_list:
                for test_HLA_name in test_HLA_names:
                    if custom_allele_check:
                        gene = test_HLA_name.split('*')[0]
                        test_HLA_seq = HLAs_default[gene][test_HLA_name]
                        seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                        print >> sys.stderr, "\t%s - %d bp (%s sequence)" % (test_HLA_name, len(test_HLA_seq), seq_type)
                        continue
                    gene = test_HLA_name.split('*')[0]
                    test_HLA_seq = HLAs[gene][test_HLA_name]
                    seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                    print >> sys.stderr, "\t%s - %d bp (%s sequence)" % (test_HLA_name, len(test_HLA_seq), seq_type)
                    
            if custom_allele_check:
                simulate_reads(HLAs_default, test_HLA_list, simulate_interval, error_rate)
            else:
                simulate_reads(HLAs, test_HLA_list, simulate_interval, error_rate)

            if "single-end" in daehwan_debug:
                read_fname = ["hla_input_1.fa"]
            else:
                read_fname = ["hla_input_1.fa", "hla_input_2.fa"]

            fastq = False
            tmp_test_passed = HLA_typing(ex_path,
                                         simulation,
                                         reference_type,
                                         test_HLA_list,
                                         partial,
                                         refHLAs,
                                         HLAs,                       
                                         HLA_names,
                                         HLA_lengths,
                                         refHLA_loci,
                                         Vars,
                                         Var_list,
                                         Links,
                                         exclude_allele_list,
                                         aligners,
                                         num_mismatch,
                                         fastq,
                                         read_fname,
                                         alignment_fname,
                                         threads,
                                         enable_coverage,
                                         best_alleles,
                                         detect_allele,
                                         HLA_fnames,
                                         verbose)

            for aligner_type, passed in tmp_test_passed.items():
                if aligner_type in test_passed:
                    test_passed[aligner_type] += passed
                else:
                    test_passed[aligner_type] = passed

                print >> sys.stderr, "\t\tPassed so far: %d/%d (%.2f%%)" % (test_passed[aligner_type], test_i + 1, (test_passed[aligner_type] * 100.0 / (test_i + 1)))


        for aligner_type, passed in test_passed.items():
            print >> sys.stderr, "%s:\t%d/%d passed (%.2f%%)" % (aligner_type, passed, len(test_list), passed * 100.0 / len(test_list))
    
    else: # With real reads or BAMs
        print >> sys.stderr, "\t", ' '.join(hla_list)
        fastq = True
        HLA_typing(ex_path,
                   simulation,
                   reference_type,
                   hla_list,
                   partial,
                   refHLAs,
                   HLAs,                       
                   HLA_names,
                   HLA_lengths,
                   refHLA_loci,
                   Vars,
                   Var_list,
                   Links,
                   exclude_allele_list,
                   aligners,
                   num_mismatch,
                   fastq,
                   read_fname,
                   alignment_fname,
                   threads,
                   enable_coverage,
                   best_alleles,
                   detect_allele,
                   HLA_fnames,
                   verbose)

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping')
    parser.add_argument("--base",
                        dest="base_fname",
                        type=str,
                        default="",
                        help="base filename for backbone HLA sequence, HLA variants, and HLA linking info")
    parser.add_argument("--default-list",
                        dest = "default_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of HLA alleles to be tested. Alleles are retrieved from default backbone data (all alleles included in backbone).")
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="A,B,C,DQA1,DQB1,DRB1",
                        help="A comma-separated list of HLA genes (default: A,B,C,DQA1,DQB1,DRB1)")
    parser.add_argument('--partial',
                        dest='partial',
                        action='store_true',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph,hisat2.linear,bowtie2.linear",
                        help="A comma-separated list of aligners (default: hisat2.graph,hisat2.linear,bowtie2.linear)")
    parser.add_argument("--reads",
                        dest="read_fname",
                        type=str,
                        default="",
                        help="Fastq read file name")
    parser.add_argument("--alignment",
                        dest="alignment_fname",
                        type=str,
                        default="",
                        help="BAM file name")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=1,
                        help="Reads simulated at every these base pairs (default: 1)")
    parser.add_argument("--coverage",
                        dest="coverage",
                        action='store_true',
                        help="Experimental purpose (assign reads based on coverage)")
    parser.add_argument("--best-alleles",
                        dest="best_alleles",
                        action='store_true',
                        help="")
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of alleles to be excluded. Enter a number N to randomly select N alleles for exclusion and N non-excluded alleles for testing (2N tested in total).")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument("--debug",
                        dest="debug",
                        type=str,
                        default="",
                        help="e.g., test_id:10,read_id:10000,basic_test")
    parser.add_argument("--detect-allele",
                        dest="detect_allele",
                        action='store_true',
                        help="Change test to detection of new alleles. Allow output of new allele sequences.")
    parser.add_argument("--error-rate",
                        dest="error_rate",
                        type=float,
                        default = -1.0,
                        help="Set read error rate (decimal).")                        


    args = parser.parse_args()
    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    args.hla_list = args.hla_list.split(',')
    if args.aligners == "":
        print >> sys.stderr, "Error: --aligners must be non-empty."
        sys.exit(1)    
    args.aligners = args.aligners.split(',')
    for i in range(len(args.aligners)):
        args.aligners[i] = args.aligners[i].split('.')
    if args.read_fname:
        args.read_fname = args.read_fname.split(',')
    else:
        args.read_fname = []
    if args.alignment_fname != "" and \
            not os.path.exists(args.alignment_fname):
        print >> sys.stderr, "Error: %s doesn't exist." % args.alignment_fname
        sys.exit(1)
    
    if len(args.default_allele_list) > 0:
        args.default_allele_list = args.default_allele_list.split(',')
        
    if len(args.exclude_allele_list) > 0:
        if args.exclude_allele_list.strip().isdigit():
            num_alleles = int(args.exclude_allele_list)
       
            HLAs_default = {}
            fname = "./IMGTHLA/hla_gen.fasta"
            for line in open(fname):
                if line.startswith(">"):
                    #print line
                    HLA_name = line.strip().split()[1]
                    HLA_gene = HLA_name.split('*')[0]
                    if not HLA_gene in HLAs_default:
                        HLAs_default[HLA_gene] = {}
                    if not HLA_name in HLAs_default[HLA_gene]:
                        HLAs_default[HLA_gene][HLA_name] = ""
                else:
                    HLAs_default[HLA_gene][HLA_name] += line.strip()

            allele_names = list(HLAs_default['A'].keys())
            random.shuffle(allele_names)
            args.exclude_allele_list = allele_names[0:num_alleles]
            args.default_allele_list = allele_names[num_alleles:2*num_alleles]
            
            args.default_allele_list = args.default_allele_list + args.exclude_allele_list
        else:
            args.exclude_allele_list = args.exclude_allele_list.split(',')
        
    debug = {}
    if args.debug != "":
        for item in args.debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                debug[key] = value
            else:
                debug[item] = 1

    random.seed(1)
    test_HLA_genotyping(args.base_fname,
                        args.reference_type,
                        args.hla_list,
                        args.partial,
                        args.aligners,
                        args.read_fname,
                        args.alignment_fname,
                        args.threads,
                        args.simulate_interval,
                        args.coverage,
                        args.best_alleles,
                        args.exclude_allele_list,
                        args.default_allele_list,
                        args.num_mismatch,
                        args.verbose,
                        args.detect_allele,
                        args.error_rate,
                        debug)
