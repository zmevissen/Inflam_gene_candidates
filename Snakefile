#!/usr/bin/env python3

bbmap = "docker://quay.io/biocontainers/bbmap:38.96--h5c4e2a8_1"
rule trim:
    input:
        'data/1_control_18S_2019_minq7.fastq'
    output:
        r1 = 'output/010_process/1_control_18S_2019_minq7.r1.fastq.gz'
    #params:
    #    adapters = '/adapters.fa'
    log:
        'output/logs/trim.1_control_18S_2019_minq7.log'
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
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'
