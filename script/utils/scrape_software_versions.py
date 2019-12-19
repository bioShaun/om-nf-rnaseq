#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
import os

regexes = {
    'nf-core/rnaseq': ['v_ngi_rnaseq.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'fastp': ['v_fastp.txt', r'fastp (\S+)'],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'STAR': ['v_star.txt', r"STAR_(\S+)"],
    'HISAT2': ['v_hisat2.txt', r"version (\S+)"],
    'Picard': ['v_markduplicates.txt', r"([\d\.]+)-SNAPSHOT"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'featureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'deepTools': ['v_deeptools.txt', r"bamCoverage (\S+)"],
    'StringTie': ['v_stringtie.txt', r"(\S+)"],
    'GffCompare': ['v_gffcompare.txt', r"gffcompare v(\S+)"],
    'Preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'RSeQC': ['v_rseqc.txt', r"read_duplication.py ([\d\.]+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()

# Search each file using its regex
for k, v in regexes.items():
    if not os.path.exists(v[0]):
        continue
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Strip STAR or HiSAT2
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ("software\tversion")
for k, v in results.items():
    print("{}\t{}".format(k, v))
