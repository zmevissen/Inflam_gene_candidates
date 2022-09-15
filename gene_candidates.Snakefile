#!/usr/bin/env python3
import pandas
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from tempfile import mkdtemp

#############
# CONTAINERS #
#############

bbmap = "docker://quay.io/biocontainers/bbmap:38.98--h5c4e2a8_1"
#TODO is this the correct container? just got the first one
star = "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"

#############
# FUNCTIONS #
#############

#TODO ask what this does?
def dontmaketempdir():
    return Path(mkdtemp(), 'tmp').resolve().as_posix()

###########
# GLOBALS #
###########
sample_table_loc = "data/sample_table/sample_table.csv"
reads_dir = 'data/fastq_repaired'

# TODO check with Tom this is correct: genomes downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/
ref_gff_url = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.gff.gz')
ref_gff = 'GCA_000001405.29_GRCh38.p14_genomic.gff.gz'

ref_fna_url = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz')  
ref_fna = 'GCA_000001405.29_GRCh38.p14_genomic.fna.gz'  



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

rule target:
    input:
        expand('output/star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
               sample=paired_sample_names), 
        expand('output/star/pass2/{sample}.ReadsPerGene.out.tab',
               sample=paired_sample_names)
        
        
rule star_second_pass:
    input:
        r1 = 'output/trim/{sample}_1.fastq.gz',
        r2 = 'output/trim/{sample}_2.fastq.gz',
        star_reference = 'output/star-index/SA',
        junctions = expand('output/star/pass1/{sample}.SJ.out.tab',
                           sample=paired_sample_names)
    output:
        bam = 'output/star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
        counts = 'output/star/pass2/{sample}.ReadsPerGene.out.tab'
    threads:
        10
    params:
        genome_dir = 'output/star-index',
        prefix = 'output/star/pass2/{sample}.'
    log:
        'output/logs/star_second_pass.{sample}.log'
    resources:
        time = 59,
        mem_mb = 32 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} ' #Do we need to use a scratch directory?
        '--sjdbFileChrStartEnd {input.junctions} '  #WHy use this instead of just the GTF file?
        '--outSAMtype BAM SortedByCoordinate '  #Explain sorted coordinate command?
        '--outReadsUnmapped Fastx '   #Explain? Is this a file that contains the reads that were unmapped?
        '--quantMode GeneCounts '   #This tells it to output counts for reads per gene but explain?
        '--readFilesCommand zcat '
        '--readFilesIn {input.r1} {input.r2}'
        '--outFileNamePrefix {params.prefix} '
        '--outTmpDir ' + dontmaketempdir() + ' ' #explain? 
        '&> {log}'
        

        

 #TODO is the purpose of the first pass to identify splice junctions? How does this work?       
 #Can STAR handle multiple inputs? Will it produce 1 or 2 junction files?
rule star_first_pass:
    input:
        r1 = 'output/trim/{sample}_1.fastq.gz',
        r2 = 'output/trim/{sample}_2.fastq.gz',
        star_reference = 'output/star/star-index/SA' #What do the numbers mean? 007
    output:
        sjdb = 'output/star/pass1/{sample}.SJ.out.tab'
    threads:
        10
    params:
        index = 'output/star/star-index',
        prefix = 'output/star/pass1/{sample}.'  #Can I remove the 025 prefix before star?
    log:
        'output/logs/star_first_pass.{sample}.log'
    resources:
        time = 59,
        mem_mb = 32 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} ' #is 10 good for us?
        '--genomeDir {params.index} '
        '--outSJfilterReads Unique ' #does this mean only uniquely mapped reads are kept?
        '--outSAMtype None '     #So we are only interested in the splice junctions in this pass?
        '--readFilesCommand zcat '  #Will zcat work or do I need to use gunzip because they are .gz?
        '--readFilesIn {input.r1} {input.r2}' #Will this work?
        '--outFileNamePrefix {params.prefix} ' #What does this do? is it just the file name before the SJ.out.tab'
        '&> {log}'
        
        
rule star_index:
    input:
        gff = f'output/ref/{ref_gff}', 
        fasta = f'output/ref/{ref_fna}'
    output:
        'output/star/star-index/SA'
    params:
        outdir = 'output/star/star-index'
    log:
        'output/logs/star_index.log'
    threads:
        10
    resources:
        time = 30,
        mem_mb = 32 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.outdir} ' #does this need to be a scratch directory?
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gff} ' 
        '--genomeSAindexNbases 12 '  #explain? why 12 instead of default 14?
        '--outTmpDir ' + dontmaketempdir() + ' ' #explain?
        '--sjdbGTFtagExonParentTranscript Parent ' #explain?
        '--sjdbGTFtagExonParentGene gene ' #explain?
        '&> {log}'
        
rule download_ref_fna:
    input:
        FTP.remote(ref_fna_url, keep_local=True)
    output:
        f'output/ref/{ref_fna}'
    log:
        'output/logs/ref_fna.log'
    shell:
        'cp {input} {output} '
    
rule download_ref_gff:
    input:
        FTP.remote(ref_gff_url, keep_local=True)
    output:
        f'output/ref/{ref_gff}'
    log:
        'output/logs/ref_gff.log'
    shell:
        'cp {input} {output} '

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
