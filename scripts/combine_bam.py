import pysam
import random
import itertools
import math
import re

min_mapq = snakemake.params.min_mapq
sam_list = [pysam.AlignmentFile(b, 'rb', check_sq=False) for b in snakemake.input.bam]
sam_comb = pysam.AlignmentFile(snakemake.output[0], "wb", template=sam_list[0])
for line_tuple in itertools.zip_longest(*sam_list):
    match_list = [-math.inf if line.is_unmapped else line.get_tag('AS')
                  for line in line_tuple]
    ## select the line(s) with the highest alignment score (nucleotide matches)
    line_list = [line_tuple[i] for i in range(0,2) if match_list[i]==max(match_list)]
    if len(line_list) == 1:
        ## if there is one alignment with better alignment, pick that one
        out_line = line_list[0]
    else:
        ## else pick lowest mapping quality, this will best represent real mapping
        ## quality.
        mapq_list = [line.mapping_quality for line in line_tuple]
        line_list = [line_list[i] for i in range(0,2) if
                     mapq_list[i]==max(mapq_list)]
        if len(line_list) == 1:
            out_line = line_list[0]
        else:
            out_line = line_list[round(random.random())]
    ## if quality is equal or higher than the minimum mapping quality, write it down
    if out_line.mapping_quality >= min_mapq:
        sam_comb.write(out_line)
[b.close() for b in sam_list]
sam_comb.close()
