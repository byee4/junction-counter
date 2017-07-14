# junction-counter

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
--jxn_file data/jxnlist.txt \
--bam ENCFF756RDZ.bam \
--outfile data/output.txt \
--min_overlap 0
```

# More stuff

jxnlist.txt contains a line-delimited list of junctions
that follow the format chr:start-stop:strand

output.txt contains the junction list, the number of reads supporting that junction,
and the read names supporting that junction.

You can download the BAM file [here](https://www.encodeproject.org/files/ENCFF756RDZ/@@download/ENCFF756RDZ.bam)
 (make sure to index this as well!)
