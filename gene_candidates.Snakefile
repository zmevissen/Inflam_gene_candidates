#!/usr/bin/env python3
import pandas
bbmap = "quay.io/biocontainers/bbmap:38.98--h5c4e2a8_1"


#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########
sample_table_loc = "data/sample_table/sample_table.csv"
reads_dir = 'data/fastq_repaired'

#########
# MAIN #
#########

sample_table = pandas.read_csv(
    sample_table_loc,
    index_col="sample_spm_name")

paired_sample_names = sorted(set(sample_table[sample_table.LibraryLayout == 'PAIRED'].index))

#########
# RULES #
#########

rule target:
    input:
        expand('output/trim/{sample}_{r}.fastq.gz',
               sample=paired_sample_names, r=[1, 2])

rule trim:
    input:
        r1 = Path(reads_dir, "{sample}_1.fastq.gz").resolve(),
        r2 = Path(reads_dir, "{sample}_2.fastq.gz").resolve()
        
    output:
        r1 = 'output/trim/{sample}_1.fastq.gz',
        r2 = 'output/trim/{sample}_2.fastq.gz'
        
    log:
        'output/logs/trim.{sample}.log'
    threads:
        10
    resources:
        time = 59,
        mem_mb = 10 * 1000
    container:
        bbmap
    shell:
        'bbduk.sh '
        '-Xmx{resources.mem_mb}m '
        'zl=9 '
        'in={input.r1} '
        'in2={input.r1} '
        'out={output.r1} '
        'out2={output.r2} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'
