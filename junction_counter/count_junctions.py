#!/usr/bin/env python

"""
basic module to use as an example.
"""
import argparse
import pandas as pd
import pysam

from tqdm import trange
from collections import defaultdict

def get_offset_m_basedon_n(
        cigartuples, n,
        include_jxn_span=True,
        include_insertions=False,
        include_deletions=True,
        verbose=False
):
    """
    for junction # n (0 based), return the left and right offsets m.
    So for example, if we have a cigar string such as 10M100N15M:
    this function will return:
    (if no flags): 10, 15
    (if include_jxn_span): 110, 115
    (if include_insertions):
    """
    all_counter = 0
    counter = 0
    left_accumulated_m = 0
    current_left_offset = 0
    current_right_offset = 0
    for t in cigartuples:
        if verbose:
            print(t)
        if t[0] == 0:
            left_accumulated_m += t[1]
        elif t[0] == 1 and include_insertions == True:  # insertion code
            left_accumulated_m += t[1]
        elif t[0] == 2 and include_deletions == True:
            left_accumulated_m += t[1]
        elif t[0] == 3:
            if include_jxn_span:
                current_left_offset = left_accumulated_m + t[1]
            else:
                current_left_offset = left_accumulated_m
            counter += 1

        if counter >= n:
            if verbose:
                print("ALL COUNTER", all_counter)
            for tr in cigartuples[all_counter:]:
                if verbose:
                    print("TR", tr)
                if tr[0] == 0:
                    current_right_offset += tr[1]
                elif tr[0] == 1 and include_insertions == True:
                    current_right_offset += tr[1]
                elif tr[0] == 2:  # and include_deletions == True:
                    current_right_offset += tr[1]
                elif tr[0] == 3 and include_jxn_span == True:
                    current_right_offset += tr[1]

            return [current_left_offset, current_right_offset]
        all_counter += 1
    return [current_left_offset, current_right_offset]


def get_jxc_counts_from_cigartuples(cigartuples):
    """
    From a list of cigar tuples, return the number of junctions found (N's)
    """
    count = 0
    for t in cigartuples:
        if t[0] == 3:
            count += 1
    return count


def parse_jxn_string(jxn_string):
    """
    Parses a line like: chr:start-stop:strand to return
    these proper fields.
    """
    chrom, pos, strand = jxn_string.split(':')
    start, stop = pos.split('-')
    return chrom, int(start), int(stop), strand


def return_spliced_junction_counts(jxn_string, bam_file, min_overlap):
    aligned_file = pysam.AlignmentFile(bam_file, "rb")
    chrom, start, stop, strand = parse_jxn_string(jxn_string)

    if start + 1 != stop:
        print(
        "WARNING: junction site is greater than 1 nt: {}".format(jxn_string))

    depth = 0  # initial number of reads supporting junction
    names = ''  # string of 'comma delimited' names of all reads supporting junction
    badflags = []  # flags of all the skipped reads, if failed QC (for debugging mostly)

    # jxn = '{}:{}-{}:{}'.format(row['chrom'],row['low_exon_end'],int(row['low_exon_end'])+1, row['strand'])
    for read in aligned_file.fetch(chrom, start,
                                   stop):  # get the reads which span this position
        ### cigar_tuple codes:
        # 0: match
        # 3: skip

        skip = False  # initially treat this read as supporting jxn
        left_span = False  # assume that this read does not support any junction yet
        right_span = False  # assume that this read does not support any junction yet

        jxc_count = 0  # number of junctions in the read
        if 'N' in read.cigarstring:  # only consider reads with jxns
            """ We need this?
            if read.get_tag('XS') != row['strand']:
                skip = True
            """
            if (
            not read.is_proper_pair) or read.is_qcfail or read.is_secondary:
                skip = True
                badflags.append(read.flag)

            jxc_count = get_jxc_counts_from_cigartuples(read.cigartuples)
            ###
            # Look at each junction that this read overlaps and see if any of them match
            # the one in focus.
            ###
            for j in range(1, jxc_count + 1):
                left, right = get_offset_m_basedon_n(
                    read.cigartuples, j, True, True, True
                )
                left_wo, right_wo = get_offset_m_basedon_n(
                    read.cigartuples, j, False, False, False
                )
                if left_wo < min_overlap or right_wo < min_overlap:
                    skip = True
                # if read.query_name == 'D80KHJN1:241:C5HF7ACXX:6:1105:6615:43587':
                #     print("RIGHT", read.reference_end, right, int(row['low_exon_end']), read.flag, read.cigartuples)
                #     print("LEFT", read.reference_start, read.query_alignment_start, left, int(row['low_exon_end']), read.flag)
                ### need to filter out read1 whose: 1) same strand
                # left span refers to this junction 'L', while right_span refers to R
                #   [L======R]-----[L=========R]----...
                if read.reference_start + read.query_alignment_start + left == start:
                    left_span = True  # read supports a left jxc
                if read.reference_end - right == start:
                    right_span = True  # read supports a right jxc

            if left_span == False and right_span == False:
                skip = True
            if not skip:
                names += read.query_name + ','
                depth += 1
    return depth, names[:-1]


def get_junction_sites(jxn_list, bam_file, min_overlap):
    jxn_dict = defaultdict(dict)
    progress = trange(len(jxn_list))
    for jxn in jxn_list:
        depth, names = return_spliced_junction_counts(jxn, bam_file, min_overlap)
        jxn_dict[jxn] = {'depth': depth, 'names': names}
        progress.update(1)
    return pd.DataFrame(jxn_dict).T


def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--jxn_file",
        required=True,
    )
    parser.add_argument(
        "--bam",
        required=True,
    )
    parser.add_argument(
        "--outfile",
        required=True,
    )
    parser.add_argument(
        "--min_overlap",
        required=False,
        default=0,
        type=int
    )
    args = parser.parse_args()

    jxn_file = args.jxn_file
    out_file = args.outfile
    bam = args.bam
    min_overlap = args.min_overlap
    jxn_list = []
    with open(jxn_file, 'r') as f:
        for line in f:
            jxn_list.append(line.rstrip())

    jxc_df = get_junction_sites(jxn_list, bam, min_overlap)
    jxc_df.to_csv(out_file, sep='\t')

if __name__ == "__main__":
    main()
