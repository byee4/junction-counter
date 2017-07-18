# junction-counter

Basically a few steps upstream of any splicing algorithm.
If you just want to count reads supporting splice junctions
given an RNA-SEQ bam file and a spliced region,
you can use this tool to do it quickly:

![Alt text](images/chr10-100182270-100183359.png)

Requirements:

```
samtools
pysam
pandas
tqdm
```

(there's a weird shared libraries problem with samtools, but
assuming you can get around that, these packages are all you need
to run).

Usage:

Refer to the bash script run_junction_counter.sh

```
python junction_counter/count_junctions.py \
--bed data/intron_jxnlist.bed \
--bam ENCFF756RDZ.bam \
--outfile data/output.txt \
--min_overlap 1
```

# More stuff
jxnlist.txt contains a bedfile of introns we want to count junctions for

output.txt contains the junction list, the number of reads supporting that junction,
and the read names supporting that junction.

You can download the BAM file [here](https://www.encodeproject.org/files/ENCFF756RDZ/@@download/ENCFF756RDZ.bam)
 (make sure to index this as well!)
