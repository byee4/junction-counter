#!/usr/bin/env python

"""
basic module to use as an example.
"""
import argparse
import pandas as pd
import pysam
from tqdm import trange
from collections import defaultdict, OrderedDict


def get_offset_m_basedon_n(
        cigartuples, n, include_jxn_span=True, include_insertions=False,
        include_deletions=True, verbose=False
):
    """
    for junction # n (0 based), return the left and right offsets m.
    So for example, if we have a cigar string such as 10M100N15M:
    this function will return:
    (if no flags): 10, 15
    (if include_jxn_span): 110, 115
    ...
    
    :param cigartuples: list
    :param n: int
    :param include_jxn_span: boolean
        if True, append the jxn span to either side
    :param include_insertions: boolean
        if True, append the insertion lengths to either side
    :param include_deletions: boolean
        if True, append the deletion lengths to either side
    :param verbose: 
        if True, print out CIGAR info and debugging stuff
    :return: 
    """

    all_counter = 0
    counter = 0
    left_accumulated_m = 0
    current_left_offset = 0
    current_right_offset = 0
    # print(cigartuples, n)
    for t in cigartuples:
        if verbose:
            print('T: ', t)
        if t[0] == 0:
            left_accumulated_m += t[1]
        elif t[0] == 1 and include_insertions == True:  # insertion code
            left_accumulated_m += t[1]
        elif t[0] == 2 and include_deletions == True:
            left_accumulated_m += t[1]
        elif t[0] == 3:
            if include_jxn_span:
                # print('adding {} to {}, {}'.format(
                #   current_left_offset, left_accumulated_m, t[1])
                # )
                current_left_offset += left_accumulated_m + t[1]
            else:
                current_left_offset += left_accumulated_m
            left_accumulated_m = 0
            counter += 1
            # print('counter = {}'.format(counter))
            # print('current offset: {}'.format(current_left_offset))
        if verbose:
            print('current offset: {}'.format(current_left_offset))
        if counter >= n:
            # print("COUNTER greater than : {}, {}".format(counter, n))
            # print("CURRENT OFFSET LEFT: {}".format(current_left_offset))
            # print("ALL COUNTER", all_counter)
            for tr in cigartuples[all_counter:]:
                # print("TR", tr)
                if tr[0] == 0:
                    current_right_offset += tr[1]
                elif tr[0] == 1 and include_insertions == True:
                    current_right_offset += tr[1]
                elif tr[0] == 2:  # and include_deletions == True:
                    current_right_offset += tr[1]
                elif tr[0] == 3 and include_jxn_span == True:
                    current_right_offset += tr[1]
            # print("RETURNING: {}, {}".format(
            #   current_left_offset, current_right_offset)
            # )
            return [current_left_offset, current_right_offset]
        all_counter += 1

    return [current_left_offset, current_right_offset]


def get_jxc_counts_from_cigartuples(cigartuples):
    """
    From a list of cigar tuples, return the number of introns found (N's) that
    the reads supports/spans.

    :param cigartuples: list
        list of cigar tuples (usually from pysam.AlignedRead.cigartuples
    :return count: int
        number of N's found.
    """
    count = 0
    for t in cigartuples:
        if t[0] == 3:
            count += 1
    return count


def parse_interval_to_jxn_string(chrom, start, end, strand):
    """
    Parses a genomic interval into a jxn-string format.
    :param interval: pybedtools.Interval
    :return jxn_string: string
    """
    # TODO: deprecate this format
    return '{}:{}-{}:{}'.format(
        chrom,
        start,
        end,
        strand
    )


def parse_jxn_string(jxn_string):
    """
    Parses a line like: chr:5'-3':strand to return
    these proper fields.
    :param jxn_string: string
    """
    chrom, pos, strand = jxn_string.split(':')
    five, three = pos.split('-')
    return chrom, int(five), int(three), strand


def passes_strand_filter(strand, read, library='reverse_pe'):
    """
    Given a library (TruSeq rstrand PE only for now), return True if
    the reads are on the proper strand.

    :param strand: string
        either '+' or '-' indicating the strand of the gene
    :param read: pysam.AlignedSegment
        see: http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
    :param library: string
        one of 'reverse_pe', 'forward_pe', 'reverse_se', forward_se'
    :return:
    """
    if library == 'reverse_pe':
        if strand == '+':
            if read.is_reverse and read.is_read2:
                return False
            if not read.is_reverse and read.is_read1:
                return False
        elif strand == '-':
            if read.is_reverse and read.is_read1:
                return False
            if not read.is_reverse and read.is_read2:
                return False
        else:
            print("wrong strand")
            return 1
        return True
    elif library == 'forward_pe':
        if strand == '+':
            if read.is_reverse and read.is_read1:
                return False
            if not read.is_reverse and read.is_read2:
                return False
        elif strand == '-':
            if read.is_reverse and read.is_read2:
                return False
            if not read.is_reverse and read.is_read1:
                return False
        else:
            print("wrong strand")
            return 1
        return True
    elif library == 'reverse_se':
        if strand == '+':
            if not read.is_reverse:
                return False
        elif strand == '-':
            if read.is_reverse:
                return False
        else:
            print("wrong strand")
            return 1
        return True
    elif library == 'forward_se':
        if strand == '+':
            if read.is_reverse:
                return False
        elif strand == '-':
            if not strand.is_reverse:
                return False
        else:
            print("wrong strand")
            return 1
        return True
    else:
        print("wrong library")
        return 1


def right_span(read, five, min_overlap):
    """
    Returns whether or a not a read that overlaps the 'right side' of the exon
    supports an inclusion or an exclusion event.

    :param read: pysam.AlignedSegment
        pysam read that spans the
    :param five: int
        the lower genomic coordinate of an intron
    :param min_overlap: int
    :return:
    """
    supports_skipping = True  # initially treat this read as supporting jxn
    supports_inclusion = True

    ### LOOK AT JUNCTIONS DENOTED BY CIGAR ###
    if 'N' in read.cigarstring:  # only consider reads with jxns
        right_span = False  # assume that this read does not support any junction yet
        jxc_count = get_jxc_counts_from_cigartuples(read.cigartuples)
        ###
        # Look at each junction that this read overlaps and see if any of them
        # line up or stop at the one in focus. 'left span' refers to a read spanning
        # the left side of an exon
        ###
        for j in range(1, jxc_count + 1):
            # print("RIGHT SPAN", j, read.query_name)
            left, right = get_offset_m_basedon_n(
                read.cigartuples, j, True, True, True
            )

            left_wo, right_wo = get_offset_m_basedon_n(
                read.cigartuples, j, False, False, False
            )
            if left_wo < min_overlap or right_wo < min_overlap:
                right_span = False
            elif read.reference_end - right == five:
                right_span = True  # read supports a right jxc

        if not right_span:
            supports_skipping = False
            # If no junctions line up, check the exon/intron boundaries
            # We assume that this read overlaps the exon by at least one base,
            # let's see if it overlaps an intron too
            # if read.query_name == 'D80KHJN1:241:C5HF7ACXX:6:2113:2908:58516':
            #     print('right', read.get_overlap(five+min_overlap, five+min_overlap+1))
            if read.get_overlap(five+min_overlap, five+min_overlap+1) == 0:
                supports_inclusion = False
            else:
                supports_inclusion = True
        else:
            supports_skipping = True
            supports_inclusion = False

    else:  # no N's found, just need to determine whether this read overlaps
        # if a read doesn't span an intron, it doesn't support skipping.
        supports_skipping = False
        if read.get_overlap(five, five+min_overlap) > 0:
            supports_inclusion = True
        else:
            supports_inclusion = False

    return supports_skipping, supports_inclusion


def left_span(read, three, min_overlap):
    """
    Returns whether or a not a read that overlaps the 'left side' of the exon
    supports an inclusion or an exclusion event.

    :param read:
    :param three:
    :param min_overlap:
    :return:
    """
    supports_skipping = True  # initially treat this read as supporting jxn
    supports_inclusion = True

    ### LOOK AT JUNCTIONS DENOTED BY CIGAR ###
    if 'N' in read.cigarstring:  # only consider reads with jxns
        left_span = False  # assume that this read does not support any jxn yet
        jxc_count = get_jxc_counts_from_cigartuples(read.cigartuples)
        ###
        # Look at each junction that this read overlaps and see
        # if any of them line up or stop at the one in focus.
        # 'left span' refers to a read spanning the left side of an exon
        ###
        for j in range(1, jxc_count + 1):
            left, right = get_offset_m_basedon_n(
                read.cigartuples, j, True, True, True
            )
            # This block handles overlaps
            # (need minimum x offset to call a site included/excluded)
            left_wo, right_wo = get_offset_m_basedon_n(
                read.cigartuples, j, False, False, False
            )
            if left_wo < min_overlap or right_wo < min_overlap:
                left_span = False

            # print(j, read.reference_start, read.query_alignment_start, left,
            # read.reference_start + read.query_alignment_start + left, three)
            elif read.reference_start + read.query_alignment_start + left == three:
                left_span = True  # read supports a left jxc

        if not left_span:
            supports_skipping = False
            # If no junctions line up, check the exon/intron boundaries
            # We assume that this read overlaps the exon by at least one base,
            # let's see if it overlaps an intron too
            if read.get_overlap(three-min_overlap, three) == 0:
                supports_inclusion = False
            else:
                supports_inclusion = True
        else:
            supports_skipping = True
            supports_inclusion = False

    else:
        # no N's found, so we just need to determine whether or not
        # this read overlaps enough. If a read doesn't span an intron,
        # it doesn't support skipping.
        supports_skipping = False
        # if read.query_name == 'D80KHJN1:241:C5HF7ACXX:6:2113:2908:58516':
        #     print('left', read.get_overlap(three - min_overlap - 1, three - min_overlap))
        if read.get_overlap(three - min_overlap - 1, three - min_overlap) > 0:
            supports_inclusion = True
        else:
            supports_inclusion = False
    return supports_skipping, supports_inclusion


def total_inclusion(read, five, three):
    """
    Returns True if the read is completely within an intron boundary.

    :param read: pysam.AlignedSegment
        read
    :param five: int
        'lower' position of an intron
    :param three: int
        'upper' position of an intron
    :return:
    """
    if read.reference_start >= five and read.reference_end <= three:
        return True
    else:
        return False


def return_spliced_junction_counts(jxn_string, bam_file, min_overlap, library):
    """
    Returns the junction inclusion and exclusion counts for a given junction.
    
    :param jxn_string: string
        in the format "chr:start-end:strand"
    :param bam_file: string
        bam filename
    :param min_overlap: int
        minimum overlap required to support site (1nt default)
    :param library: string
        default 'reverse_pe'
    :return: 
    """
    aligned_file = pysam.AlignmentFile(bam_file, "rb")
    chrom, five, three, strand = parse_jxn_string(jxn_string)

    skip_names_list = []
    incl_names_list = []
    skip_names = ''  # 'comma delimited' skip_names of reads supporting jxn
    incl_names = ''

    for read in aligned_file.fetch(
            chrom, five - min_overlap, three + min_overlap
    ):  # five denotes 0-based start of intron

        ### WARNING THIS ASSUMES TRUSEQ PAIRED END ###
        ### INITIAL READ QC CHECKS ###
        if (not passes_strand_filter(strand, read, library)) or \
                (not read.is_proper_pair) or \
                read.is_qcfail or \
                read.is_secondary:
            pass
        else:
            # does this read support inclusion or skipping from the right side
            # of the exon/left side of the intron?
            supports_skipping_right, supports_inclusion_right = right_span(
                read, five, min_overlap
            )
            # does this read support inclusion or skipping from the left side
            # of the exon/right side of the intron?
            supports_skipping_left, supports_inclusion_left = left_span(
                read, three, min_overlap
            )
            # does this read support inclusion/is located between two exons?
            supports_total_inclusion = total_inclusion(
                read, five, three
            )

            # read must support both the left and right side of a jxc
            if supports_skipping_right and supports_skipping_left:
                skip_names_list.append(read.query_name)
            # read originating from one exon need only to be from one side.
            elif supports_inclusion_right or supports_inclusion_left:
                incl_names_list.append(read.query_name)
            # read supports inclusion if it's between two junctions
            elif supports_total_inclusion:
                incl_names_list.append(read.query_name)

    # Removing duplicated skip_names will filter
    # double counting R1/R2 if they span the same region.
    for name in set(skip_names_list):
        skip_names += name + ','
    for name in set(incl_names_list):
        incl_names += name + ','

    return len(set(skip_names_list)), len(set(incl_names_list)), skip_names[:-1], incl_names[:-1]


def get_junction_sites(jxn_list, bam_file, min_overlap):
    """
    returns the depth (number of reads supporting) a skipped event,
    the depth of an inclusion event, and the names of both events.

    :param jxn_list:
    :param bam_file:
    :param min_overlap:
    :return:
    """
    jxn_dict = defaultdict(dict)
    progress = trange(len(jxn_list))
    for jxn in jxn_list:
        skip_depth, incl_depth, skip_names, incl_names = return_spliced_junction_counts(
            jxn, bam_file, min_overlap
        )
        jxn_dict[jxn] = OrderedDict(
            {
                'skip_depth': skip_depth, 'incl_depth': incl_depth,
                'skip_names': skip_names, 'incl_names': incl_names
            }
        )


        progress.update(1)
    return pd.DataFrame(jxn_dict).T


def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bed",
        required=True,
        help="a bed file containing intron/spliced exon positions to check"
    )
    parser.add_argument(
        "--bam",
        required=True,
        help='a bam file with an index (.bai) in the same location'
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help='outputs a tab separated file containing depth and read names for'
             'each splice event'
    )
    parser.add_argument(
        "--min_overlap",
        required=False,
        default=1,
        type=int
    )
    parser.add_argument(
        "--library",
        required=False,
        help='default "reverse_pe". Either reverse_pe, forward_pe, '
             'reverse_se, forward_se'
    )
    args = parser.parse_args()

    jxn_file = args.bed
    out_file = args.outfile
    bam = args.bam
    min_overlap = args.min_overlap
    library = args.library

    jxn_list = []
    with open(jxn_file, 'r') as f:
        for line in f:
            chrom, start, end, _, _, strand = line.rstrip().split('\t')
            jxn_string = parse_interval_to_jxn_string(
                chrom, start, end, strand
            )
            jxn_list.append(jxn_string)

    jxc_df = get_junction_sites(jxn_list, bam, min_overlap, library)
    jxc_df.to_csv(out_file, sep='\t')

if __name__ == "__main__":
    main()
