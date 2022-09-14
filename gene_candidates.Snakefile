#!/usr/bin/env python3
import pandas
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

#############
# CONTAINERS #
#############

bbmap = "docker://quay.io/biocontainers/bbmap:38.98--h5c4e2a8_1"


#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########
sample_table_loc = "data/sample_table/sample_table.csv"
reads_dir = 'data/fastq_repaired'
# TODO check with Tom this is correct: genomes downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/
ref_gff_url = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.gff.gz')
ref_gff = 'GCA_000001405.29_GRCh38.p14_genomic.gff.gz'

#########
# MAIN #
#########
FTP = FTPRemoteProvider()

sample_table = pandas.read_csv(
    sample_table_loc,
    index_col="sample_spm_name")

paired_sample_names = sorted(set(sample_table[sample_table.LibraryLayout == 'PAIRED'].index))

#########
# RULES #
#########

# rule target:
#     input:
#         expand('output/star/pass1/{sample}_{r}.fastq.gz',
#                sample=paired_sample_names, r=[1, 2])
        
        
        
        
        
# rule star_second_pass:
#     input:
#         r1 = 'output/010_process/{sample}.r1.fastq.gz',
#         star_reference = 'output/007_star-index/SA',
#         junctions = expand('output/025_star/pass1/{sample}.SJ.out.tab',
#                            sample=all_samples)
#     output:
#         bam = 'output/025_star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
#         counts = 'output/025_star/pass2/{sample}.ReadsPerGene.out.tab'
#     threads:
#         10
#     params:
#         genome_dir = 'output/007_star-index',
#         prefix = 'output/025_star/pass2/{sample}.'
#     log:
#         'output/logs/star_second_pass.{sample}.log'
#     resources:
#         time = 59,
#         mem_mb = 32 * 1000
#     container:
#         star
#     shell:
#         'STAR '
#         '--runThreadN {threads} '
#         '--genomeDir {params.genome_dir} '
#         '--sjdbFileChrStartEnd {input.junctions} '
#         '--outSAMtype BAM SortedByCoordinate '
#         '--outReadsUnmapped Fastx '
#         '--quantMode GeneCounts '
#         '--readFilesCommand zcat '
#         '--readFilesIn {input.r1} '
#         '--outFileNamePrefix {params.prefix} '
#         '--outTmpDir ' + dontmaketempdir() + ' '
#         '&> {log}'
        
        
        
        
# rule star_first_pass:
#     input:
#         r1 = 'output/trim/{sample}_1.fastq.gz',
#         r2 = 'output/trim/{sample}_2.fastq.gz'
#         star_reference = 'output/007_star-index/SA'
#     output:
#         sjdb = 'output/025_star/pass1/{sample}.SJ.out.tab'
#     threads:
#         10
#     params:
#         genome_dir = 'output/007_star-index',
#         prefix = 'output/025_star/pass1/{sample}.'
#     log:
#         'output/logs/star_first_pass.{sample}.log'
#     resources:
#         time = 59,
#         mem_mb = 32 * 1000
#     container:
#         star
#     shell:
#         'STAR '
#         '--runThreadN {threads} '
#         '--genomeDir {params.genome_dir} '
#         '--outSJfilterReads Unique '
#         '--outSAMtype None '          # troubleshoot gtf
#         # '--outSAMtype SAM '               # troubleshoot gtf
#         # '--quantMode GeneCounts '       # troubleshoot gtf
#         '--readFilesCommand zcat '
#         '--readFilesIn {input.r1} '
#         '--outFileNamePrefix {params.prefix} '
#         '&> {log}'


        
        
# rule star_index:
#     input:
#         fasta = ref,
#         gff = gff
#     output:
#         'output/007_star-index/SA'
#     params:
#         outdir = 'output/007_star-index'
#     log:
#         'output/logs/star_index.log'
#     threads:
#         10
#     resources:
#         time = 30,
#         mem_mb = 32 * 1000
#     container:
#         star
#     shell:
#         'STAR '
#         '--runThreadN {threads} '
#         '--runMode genomeGenerate '
#         '--genomeDir {params.outdir} '
#         '--genomeFastaFiles {input.fasta} '
#         '--sjdbGTFfile {input.gff} '
#         '--genomeSAindexNbases 12 '
#         '--outTmpDir ' + dontmaketempdir() + ' '
#         '--sjdbGTFtagExonParentTranscript Parent '
#         '--sjdbGTFtagExonParentGene gene '
#         # '--sjdbGTFtagExonParentGeneName Name '
#         '&> {log}'
   
rule target:
    input:
        f'output/ref/{ref_gff}'
    
rule download_ref_gff:
    input:
        FTP.remote(ref_gff_url, keep_local=True)
    output:
        f'output/ref/{ref_gff}'
    log:
        'output/logs/ref_gff.log'
    shell:
        'mv {input} {output} '

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
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'
