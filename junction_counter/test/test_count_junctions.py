from junction_counter import count_junctions as j
import pysam

"""
To create a single read bam file:

samtools view data/ENCFF468DJE.bam | grep 'DF8F08P1:365:C60G9ACXX:5:1208:17397:64846' > \
data/tmp.sam;
cat data/header.txt data/tmp.sam | samtools view \
-bS > junction_counter/test/bams/multi_junction_skip_2.bam;
samtools index junction_counter/test/bams/multi_junction_skip_2.bam
"""
def single_jxn_skip_overlapping_mate_1():
    """
    D80KHJN1:241:C5HF7ACXX:6:1314:9778:90604
    :return:
    """
    return 'bams/single_junction_skip_overlapping_mate_1.bam'

def single_jxn_skip_overlapping_mate_2():
    """
    D80KHJN1:241:C5HF7ACXX:6:1310:8490:70959
    :return:
    """
    return 'bams/single_junction_skip_overlapping_mate_2.bam'

def single_jxn_skip_1():
    """
    D80KHJN1:241:C5HF7ACXX:6:1307:9279:87799
    :return:
    """
    return 'bams/single_junction_skip_1.bam'

def single_jxn_skip_2():
    """
    D80KHJN1:241:C5HF7ACXX:6:2207:4811:73625
    :return:
    """
    return 'bams/single_junction_skip_2.bam'

def single_jxn_skip_3():
    """
    D80KHJN1:241:C5HF7ACXX:6:1309:15137:12983
    :return:
    """
    return 'bams/single_junction_skip_3.bam'

def single_jxn_skip_4():
    """
    D80KHJN1:241:C5HF7ACXX:6:2113:6890:67612
    :return:
    """
    return 'bams/single_junction_skip_4.bam'


def single_jxn_no_skip_1():
    """
    D80KHJN1:241:C5HF7ACXX:6:1110:6416:79053
    :return:
    """
    return 'bams/single_junction_no_skip_1.bam'

def single_jxn_no_skip_2():
    """
    D80KHJN1:241:C5HF7ACXX:6:1314:11862:13068
    :return:
    """
    return 'bams/single_junction_no_skip_2.bam'

def multi_jxn_skip_1():
    """
    D80KHJN1:241:C5HF7ACXX:6:2307:8849:50473
    :return:
    """
    return 'bams/multi_junction_skip_1.bam'

def multi_jxn_skip_2():
    """
    DF8F08P1:365:C60G9ACXX:5:1208:17397:64846
    :return:
    """
    return 'bams/multi_junction_skip_2.bam'

def inclusion_1():
    """
    D80KHJN1:241:C5HF7ACXX:6:1213:6350:10569
    :return:
    """
    return 'bams/inclusion_1.bam'

def inclusion_2():
    """
    D80KHJN1:241:C5HF7ACXX:6:2103:13975:89570
    :return:
    """
    return 'bams/inclusion_2.bam'

def inclusion_3():
    """
    D80KHJN1:241:C5HF7ACXX:6:2109:20023:77241
    :return:
    """
    return 'bams/inclusion_3.bam'

def inclusion_4():
    """
    D80KHJN1:241:C5HF7ACXX:6:1301:3675:25338
    :return:
    """
    return 'bams/inclusion_4.bam'

def inclusion_5():
    """
    D80KHJN1:241:C5HF7ACXX:6:2114:8225:47562
    :return:
    """
    return 'bams/inclusion_5.bam'

def inclusion_6():
    """
    D80KHJN1:241:C5HF7ACXX:6:2309:19086:46598
    :return:
    """
    return 'bams/inclusion_6.bam'

def inclusion_7():
    """
    D80KHJN1:241:C5HF7ACXX:6:2113:2908:58516
    :return:
    """
    return 'bams/inclusion_7.bam'

def indel_1():
    """
    DF8F08P1:365:C60G9ACXX:5:2313:8791:67629
    (From ENCFF468DJE.bam)
    :return:
    """
    return 'bams/indel_1.bam'

def test_right_span_1():
    print("Tests a single read with a single junction N")
    print("Region: chr10:135209281-135209666:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1307:9279:87799")
    print("Should support skipping (read spans this region(intron)")
    bam = single_jxn_skip_1()
    jxc_coord = 'chr10:135209281-135209666:+'
    min_overlap = 1
    aligned_file = pysam.AlignmentFile(bam, "rb")
    chrom, five, three, strand = j.parse_jxn_string(jxc_coord)
    for read in aligned_file.fetch(chrom, five-min_overlap, five):
        skip, incl = j.right_span(read, five, min_overlap)
        assert skip == True
        assert incl == False

def test_left_span_1():
    print("Tests a single read with a single junction N")
    print("Region: chr10:135209281-135209666:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1307:9279:87799")
    print("Should support skipping (read spans this region(intron)")
    bam = single_jxn_skip_1()
    jxc_coord = 'chr10:135209281-135209666:+'
    min_overlap = 1
    aligned_file = pysam.AlignmentFile(bam, "rb")
    chrom, five, three, strand = j.parse_jxn_string(jxc_coord)
    for read in aligned_file.fetch(chrom, five - min_overlap, five):
        skip, incl = j.left_span(read, three, min_overlap)
        assert skip == True
        assert incl == False

def test_single_junction_skip_1():
    print("Tests a single read with a single junction N")
    print("Region: chr10:135209281-135209666:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1307:9279:87799")
    print("Should support skipping (read spans this region(intron)")
    bam = single_jxn_skip_1()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr10:135209281-135209666:+'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only
    )
    assert sk == 1
    assert skn == 'D80KHJN1:241:C5HF7ACXX:6:1307:9279:87799'
    assert inc == 0


def test_right_span_2():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2207:4811:73625")
    print("Should support skipping (read spans this region(intron)")
    bam = single_jxn_skip_2()
    jxc_coord = 'chr10:122666367-122668067:+'
    min_overlap = 1
    aligned_file = pysam.AlignmentFile(bam, "rb")
    chrom, five, three, strand = j.parse_jxn_string(jxc_coord)
    for read in aligned_file.fetch(chrom, five-min_overlap, five):
        skip, incl = j.right_span(read, five, min_overlap)
        assert skip == True
        assert incl == False


def test_left_span_2():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2207:4811:73625")
    print("Should support skipping (read spans this region(intron)")
    bam = single_jxn_skip_2()
    jxc_coord = 'chr10:122666367-122668067:+'
    min_overlap = 1
    aligned_file = pysam.AlignmentFile(bam, "rb")
    chrom, five, three, strand = j.parse_jxn_string(jxc_coord)
    for read in aligned_file.fetch(chrom, five - min_overlap, five):
        skip, incl = j.left_span(read, three, min_overlap)
        assert skip == True
        assert incl == False


def test_single_junction_skip_2():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2207:4811:73625")
    print("Should not support skipping "
          "(read span is less than min_overlap (11))")
    bam = single_jxn_skip_2()
    min_overlap = 11
    jxc_only = False
    jxc_coord = 'chr10:122666367-122668067:+'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 0

def test_single_junction_skip_2_2():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2207:4811:73625")
    print("Should support skipping "
          "(read span covers min_overlap (10))")
    bam = single_jxn_skip_2()
    min_overlap = 10
    jxc_only = False
    jxc_coord = 'chr10:122666367-122668067:+'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 1
    assert inc == 0


def test_single_junction_skip_2_3():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2207:4811:73625")
    print("Should support skipping "
          "(read span covers min_overlap (10))"
          "(jxn_only = True should not change skipping result.")
    bam = single_jxn_skip_2()
    min_overlap = 10
    jxc_only = True
    jxc_coord = 'chr10:122666367-122668067:+'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 1
    assert inc == 0


def test_right_span_4():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2113:6890:67612")
    print("Should support skipping (read spans this region(intron)")
    bam = single_jxn_skip_4()
    jxc_coord = 'chr10:122666367-122668067:+'
    min_overlap = 1
    aligned_file = pysam.AlignmentFile(bam, "rb")
    chrom, five, three, strand = j.parse_jxn_string(jxc_coord)
    for read in aligned_file.fetch(chrom, five-min_overlap, five):
        skip, incl = j.right_span(read, five, min_overlap)
        assert skip == True
        assert incl == False


def test_left_span_4():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2113:6890:67612")
    print("Should support skipping (read spans this region(intron)")
    bam = single_jxn_skip_4()
    jxc_coord = 'chr10:122666367-122668067:+'
    min_overlap = 8
    aligned_file = pysam.AlignmentFile(bam, "rb")
    chrom, five, three, strand = j.parse_jxn_string(jxc_coord)
    for read in aligned_file.fetch(chrom, five - min_overlap, five):
        skip, incl = j.left_span(read, three, min_overlap)
        assert skip == False
        assert incl == False


def test_single_junction_skip_4():
    print("Tests a single read with a single junction N")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2113:6890:67612")
    print("Should not support skipping (read span is less than min_overlap (8))")
    bam = single_jxn_skip_4()
    min_overlap = 8
    jxc_only = False
    jxc_coord = 'chr10:122666367-122668067:+'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 0

def test_single_junction_no_skip_1():
    print("Tests a single read with a single junction N")
    print("Region: chr16:83842351-83842909:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1110:6416:79053")
    print("Should be filtered (flags 161/81 indicate "
          "they are not properly paired")
    bam = single_jxn_no_skip_1()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr16:83842351-83842909:+'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 0


def test_single_junction_no_skip_2():
    print("Tests a single read with a single junction N")
    print("Region: chr1:40319740-40322948:-")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1314:11862:13068")
    print("Should be filtered as reads support junctions but not this one")
    bam = single_jxn_no_skip_2()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr1:40319740-40322948:-'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 0


def test_single_junction_no_skip_3():
    pass


def test_single_junction_overlapping_mate_1():
    print("Tests read with a single junction whose r1 and r2 mate overlap")
    print("Region: chr16:83842351-83842909:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1314:9778:90604")
    print("Should return just 1 read supporting this skipping, "
          "we dont want to double count overlapping mates")
    bam = single_jxn_skip_overlapping_mate_1()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr16:83842351-83842909:+' # ENST00000433866.2_intron_1_0_chr16_83842352_f
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 1
    assert skn == 'D80KHJN1:241:C5HF7ACXX:6:1314:9778:90604'
    assert inc == 0

def test_single_junction_overlapping_mate_2():
    print("Tests a single read with a single junction N")
    print("Region: chr1:40319740-40322948")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1310:8490:70959")
    print("Should return just 1 read supporting this skipping, "
          "we dont want to double count overlapping mates")
    bam = single_jxn_skip_overlapping_mate_2()
    min_overlap = 1
    jxc_only = False
    library = 'reverse_pe'
    jxc_coord = 'chr1:40319740-40322948:-'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 1
    assert skn == 'D80KHJN1:241:C5HF7ACXX:6:1310:8490:70959'
    assert inc == 0

def test_multi_junction_skip_1():
    print("Tests whether we can properly count "
          "reads containing multiple junctions (N)")
    print("Region: chr10:135212730-135213032:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2307:8849:50473")
    print("This read supports this junction but also others, "
          "so should call this as +1 for skipping")
    bam = multi_jxn_skip_1()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr10:135212730-135213032:+'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )

    assert sk == 1
    assert skn == 'D80KHJN1:241:C5HF7ACXX:6:2307:8849:50473'
    assert inc == 0

def test_inclusion_1():
    print("Tests whether a read that extends past the junction is "
          "correctly called as +1 inclusion")
    print("Region: chr1:40319740-40322948:-")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1213:6350:10569")
    print("Read should be called as included and not excluded.")
    bam = inclusion_1()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr1:40319740-40322948:-'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 1
    assert incn == 'D80KHJN1:241:C5HF7ACXX:6:1213:6350:10569'

def test_inclusion_2():
    print("Tests whether a read that is completely contained within "
          "a jxn region is correctly called as +1 inclusion")
    print("Region: chr1:40319740-40322948:-")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2103:13975:89570")
    print("Read should be called as included and not excluded.")
    bam = inclusion_2()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr1:40319740-40322948:-'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 1
    assert incn == 'D80KHJN1:241:C5HF7ACXX:6:2103:13975:89570'

def test_inclusion_4():
    print("Tests whether a read that is completely contained within "
          "a jxn region is correctly called as +1 inclusion")
    print("Region: chr1:40319740-40322948:-")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:1301:3675:25338")
    print("Read should be called as included and not excluded.")
    bam = inclusion_4()
    min_overlap = 1
    jxc_only = False
    jxc_coord = 'chr1:40319740-40322948:-'
    library = 'reverse_pe'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 1
    assert incn == 'D80KHJN1:241:C5HF7ACXX:6:1301:3675:25338'

def test_inclusion_5():
    print("Tests whether a read that is completely within the junction span"
          "is called as inclusion (if we are not just counting junction "
          "spanning reads (jxc_only is false))")
    print("Region: chr10:14883239-14884108:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2114:8225:47562")
    print("Read should be called as included and not excluded.")
    bam = inclusion_5()
    min_overlap = 1
    jxc_only = False
    library = 'reverse_pe'
    jxc_coord = 'chr10:14883239-14884108:+'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 1
    assert incn == 'D80KHJN1:241:C5HF7ACXX:6:2114:8225:47562'

def test_inclusion_6():
    print("Tests whether a read that is completely within the junction span"
          "is called as nothing (if we are JUST counting junction "
          "spanning reads (jxc_only is True))")
    print("Region: chr10:14882156-14884108:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2114:8225:47562")
    print("Read should be called as nothing.")
    bam = inclusion_5()
    min_overlap = 1
    jxc_only = True
    library = 'reverse_pe'
    jxc_coord = 'chr10:14882156-14884108:+'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 0

def test_right_span_inclusion_7():
    print("Tests a inclusion read that doesn't quite overlap over "
          "the min_overlap reqs")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2113:2908:58516")
    print("Should support nothing (would support inclusion "
          "if min overhang wasn't 34")
    bam = inclusion_7()
    min_overlap = 34
    jxc_coord = 'chr10:122666367-122668067:+'
    aligned_file = pysam.AlignmentFile(bam, "rb")
    chrom, five, three, strand = j.parse_jxn_string(jxc_coord)
    for read in aligned_file.fetch(chrom, five-min_overlap, five):
        skip, incl = j.right_span(read, five, min_overlap)
        assert skip == False
        assert incl == False

def test_inclusion_7():
    print("Tests whether a read that is an inclusion read "
          "(r2 does not meet the min overlap but r1 does) "
          " is still called an inclusion read")
    print("Region: chr10:122666367-122668067:+")
    print("Read: D80KHJN1:241:C5HF7ACXX:6:2113:2908:58516")
    print("Read should be called as included and not excluded.")
    bam = inclusion_7()
    min_overlap = 34
    jxc_only = False
    library = 'reverse_pe'
    jxc_coord = 'chr10:122666367-122668067:+'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 0
    assert inc == 1

### TEST INDEL ###

def test_indel_1():
    print("Tests whether a read with indels (2I) is correctly"
          "counted as an exclusion read")
    print("Region: chr9:108151389-108151467:+")
    print("Read: DF8F08P1:365:C60G9ACXX:5:2313:8791:67629")
    print("Read should be called as excluded and not included.")
    bam = indel_1()
    min_overlap = 10
    jxc_only = False
    library = 'reverse_pe'
    jxc_coord = 'chr9:108151389-108151467:+'
    sk, inc, skn, incn, _, _ = j.return_spliced_junction_counts(
        jxc_coord, bam, min_overlap, jxc_only, library
    )
    assert sk == 1
    assert inc == 0


### TEST GET_OFFSET_M_BASEDON_N ###

def test_get_offset_m_baseon_n_1():
    print("Ignoring in/del for now, looks at the appropriate junction "
          "(n=1) and determines the left/right offsets.")
    print("cigartuple: [(0, 16), (3, 654), (0, 57), (3, 302), (0, 27)]")
    cigartuples = [(0, 16), (3, 654), (0, 57), (3, 302), (0, 27)]
    n = 1
    include_jxn_span = True
    include_insertions = True
    include_deletions = True
    verbose = False
    left, right = j.get_offset_m_basedon_n(
        cigartuples, n, include_jxn_span,
        include_insertions, include_deletions, verbose
    )
    assert left == 670
    assert right == 1040

def test_get_offset_m_baseon_n_2():
    print("Ignoring in/del for now, looks at the appropriate junction (n=2) "
          "and determines the left/right offsets.")
    print("cigartuple: [(0, 16), (3, 654), (0, 57), (3, 302), (0, 27)]")
    cigartuples = [(0, 16), (3, 654), (0, 57), (3, 302), (0, 27)]
    n = 2
    include_jxn_span = True
    include_insertions = True
    include_deletions = True
    verbose = False
    left, right = j.get_offset_m_basedon_n(
        cigartuples, n, include_jxn_span,
        include_insertions, include_deletions, verbose
    )
    assert left == 1029
    assert right == 329

def test_get_offset_m_baseon_n_3():
    print("Ignoring in/del for now, looks at the appropriate junction (n=1) "
          "and determines the left/right offsets.")
    print("cigartuple: [(0, 10), (3, 15), (0, 20)]")
    cigartuples = [(0, 10), (3, 15), (0, 20)]
    n = 1
    include_jxn_span = True
    include_insertions = True
    include_deletions = True
    verbose = False
    left, right = j.get_offset_m_basedon_n(
        cigartuples, n, include_jxn_span,
        include_insertions, include_deletions, verbose
    )
    assert left == 25
    assert right == 35

def test_get_offset_m_baseon_n_4():
    print("Ignoring in/del for now, looks at the appropriate junction (n=1) "
          "and determines the left/right offsets.")
    print("cigartuple: [(0, 1), (3, 3708), (0, 99)],")
    cigartuples = [(0, 1), (3, 3708), (0, 99)]
    n = 1
    include_jxn_span = True
    include_insertions = True
    include_deletions = True
    verbose = True
    left, right = j.get_offset_m_basedon_n(
        cigartuples, n, include_jxn_span,
        include_insertions, include_deletions, verbose
    )
    assert left == 3709
    assert right == 3807

def test_get_offset_m_baseon_n_5():
    print("Ignoring in/del for now, looks at the appropriate junction (n=1) "
          "and determines the left/right offsets.")
    print("cigartuple: [(0, 74), (3, 475), (0, 29)],")
    cigartuples = [(0, 74), (3, 475), (0, 29)]
    n = 1
    include_jxn_span = True
    include_insertions = True
    include_deletions = True
    verbose = True
    left, right = j.get_offset_m_basedon_n(
        cigartuples, n, include_jxn_span,
        include_insertions, include_deletions, verbose
    )
    assert left == 549
    assert right == 504