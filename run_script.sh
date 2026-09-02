#!/bin/bash

# ml purge
# ml gatk/4.6.2.0
# python3 requires: pip install pysam

CROSS_CHECK='src/crosscheck_fingerprints.py'
WORKDIR='workdir'
PANEL='/path/to/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt'

python3 ${CROSS_CHECK} \
  --input-list ${WORKDIR}/sample_full_path.txt \
  --haplotype-map ${PANEL} \
  --output cross_all_sample.csv \
  --group-output inferred_group.csv \
  --workdir ${WORKDIR} \
  --crosscheck-by FILE \
  --id-suffix "sorted.bam" \
  --no-docker \
  --num-threads 2
