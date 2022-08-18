#!/usr/bin/env python3

#TODO check with Tom that this is the correct link 
bbmap = quay.io/biocontainers/bbmap:38.98--h5c4e2a8_1


#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    sample = wildcards.sample
    my_filename = sample_table.loc[sample]['filename']
    return(Path(reads_dir, my_filename).resolve().as_posix())


#########
# RULES #
#########





rule trim:
    input:
        get_reads
    output:
        'output/trim/{sample}.fastq.gz'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/trim.{sample}.log'
    threads:
        1
    resources:
        time = 59,
        mem_mb = 10 * 1000
    container:
        bbmap
    shell:
        'bbduk.sh '
        '-Xmx{resources.mem_mb}m '
        'zl=9 '
        'in={input} '
        'out={output.r1} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'
