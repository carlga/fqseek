import os
import re
import glob

from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d")

# URLs for data retrieval
SORTMERNADB_LINK = 'https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz'
KRAKEN2_STANDARD16_LINK = 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20221209.tar.gz'
KRAKEN2_PLUSPFP16_LINK = 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz'

# Set targets
RUN_RULES = []
INIT_HISAT2 = INIT_SORTMERNADB = INIT_KALLISTO = INIT_KRAKEN2 = None

if config["init"]["hisat2_ref"] and config["init"]["hisat2_index_base"]:
    INIT_HISAT2 = 'hisat2_index/' + config['init']['hisat2_index_base'] + '.1.ht2'
    RUN_RULES.append(INIT_HISAT2)

if config["init"]["sortmerna_ref"]:
    INIT_SORTMERNADB = 'sortmerna_rRNA_ref/rRNA_ref_complete'
    RUN_RULES.append(INIT_SORTMERNADB)

if config["init"]["kallisto_ref"] and config["init"]["kallisto_index"]:
    INIT_KALLISTO = 'kallisto_index/' + config['init']['kallisto_index']
    RUN_RULES.append(INIT_KALLISTO)

if config["init"]["kraken2_db_type"] and config["init"]["bracken_read_len"]:
    INIT_KRAKEN2 = (
        'kraken2_db/' + 
        config['init']['kraken2_db_type'] + 
        '/database' +
        str(config["init"]["bracken_read_len"]) +
        'mers.kmer_distrib'
    )
    RUN_RULES.append(INIT_KRAKEN2)

# Distribute allocated jobs/cores
AVAIL_THREADS = config['init']['jobs']
INIT_HISAT2_THREADS = INIT_SORTMERNADB_THREADS = INIT_KALLISTO_THREADS = INIT_KRAKEN2_THREADS = 0

if INIT_SORTMERNADB:
    INIT_SORTMERNADB_THREADS = 1
    AVAIL_THREADS -= 1

if INIT_KALLISTO:
    INIT_KALLISTO_THREADS = 1
    AVAIL_THREADS -= 1

if INIT_HISAT2:
    INIT_HISAT2_THREADS = round(AVAIL_THREADS/2) if INIT_KRAKEN2 else AVAIL_THREADS
    AVAIL_THREADS -= INIT_HISAT2_THREADS

if INIT_KRAKEN2:
    if config["init"]["kraken2_db_type"] in ['std16_prebuilt', 'plusPFP16_prebuilt']:
        INIT_KRAKEN2_THREADS = 1
        AVAIL_THREADS -= 1
    else:
        INIT_KRAKEN2_THREADS = AVAIL_THREADS
        AVAIL_THREADS = 0
    INIT_HISAT2_THREADS += AVAIL_THREADS


rule all:
    input:
        RUN_RULES


#
# Build HISAT2 index files
#
if INIT_HISAT2 and config['init']['hisat2_gtf']:

    # Include splice sites and exons
    rule hisat2_index:
        """
        Builds index for HISAT2.
        """
        input:
            ref = config['init']['hisat2_ref'],
            gtf = config['init']['hisat2_gtf']
        output:
            INIT_HISAT2
        params:
            ref = os.path.basename(config['init']['hisat2_ref']),
            gtf = os.path.basename(config['init']['hisat2_gtf']),
            base = config['init']['hisat2_index_base'],
            outdir = "hisat2_index"
        log: "logs/" + TIMESTAMP + "_hisat2_index.log"
        threads: INIT_HISAT2_THREADS
        # resources: cpus=4, mem_mb=4000, time_min=200
        shell:
            '''
            echo 'Extracting splicesites from '{params.gtf}'...' > {log}

            hisat2_extract_splice_sites.py \
                {input.gtf} \
                > {params.outdir}/{params.gtf}.ss.tsv \
                2>> {log}

            echo 'Finished extracting splicesites' >> {log}


            echo 'Extracting exons from '{params.gtf}'...' >> {log}

            hisat2_extract_exons.py \
                {input.gtf} \
                > {params.outdir}/{params.gtf}.exons.tsv \
                2>> {log}

            echo 'Finished extracting exons' >> {log}


            echo 'Building HISAT2 index using '{params.ref}'...' >> {log}

            hisat2-build \
                -p {threads} \
                --ss {params.outdir}/{params.gtf}.ss.tsv \
                --exon {params.outdir}/{params.gtf}.exons.tsv \
                {input.ref} {params.outdir}/{params.base} \
                &>> {log}

            echo 'Finished HISAT2 index' >> {log}
            '''

elif INIT_HISAT2:

    rule hisat2_index:
        """
        Builds index for HISAT2.
        """
        input:
            ref = config['init']['hisat2_ref']
        output:
            INIT_HISAT2
        params:
            ref = os.path.basename(config['init']['hisat2_ref']),
            base = config['init']['hisat2_index_base'],
            outdir = "hisat2_index"
        log: "logs/" + TIMESTAMP + "_hisat2_index.log"
        threads: INIT_HISAT2_THREADS
        # resources: cpus=4, mem_mb=4000, time_min=200
        shell:
            '''
            echo 'Building HISAT2 index using '{params.ref}'...' > {log}

            hisat2-build \
                -p {threads} \
                {input.ref} \
                {params.outdir}/{params.base} \
                &>> {log}

            echo 'Finished HISAT2 index' >> {log}
            '''


#
# Download SortMeRNA rRNA reference files
#
if INIT_SORTMERNADB:

    rule download_sortmerna_db:
            """
            Retrieves sortmerna rRNA reference database.
            """
            output:
                complete = temp(INIT_SORTMERNADB)
            params:
                link = SORTMERNADB_LINK,
                filename = os.path.basename(SORTMERNADB_LINK),
                outdir = 'sortmerna_rRNA_ref'
            log: "logs/" + TIMESTAMP + "_download_sortmerna_db.log"
            threads: INIT_SORTMERNADB_THREADS
            # resources: cpus=4, mem_mb=4000, time_min=200
            shell:
                '''
                echo 'Downloading sortmerna rRNA reference database...' > {log}

                curl -L {params.link} --remote-name --silent
                tar -zxf {params.filename} -C {params.outdir}
                rm {params.filename}

                echo 'Finished downloading rRNA reference' >> {log}

                touch {output.complete}
                '''


#
# Builds kallisto index file
#
if INIT_KALLISTO:

    rule kallisto_index:
            """
            Builds index for kallisto.
            """
            input:
                ref = config['init']['kallisto_ref']
            output:
                idx = INIT_KALLISTO
            params:
                ref = os.path.basename(config['init']['kallisto_ref'])
            log: "logs/" + TIMESTAMP + "_kallisto_index.log"
            threads: INIT_KALLISTO_THREADS
            # resources: cpus=4, mem_mb=4000, time_min=200
            shell:
                '''
                echo 'Building kallisto index using '{params.ref}'...' > {log}

                kallisto index --index {output.idx} {input.ref} &>> {log}

                echo 'Finished kallisto index' >> {log}
                '''


#
# Builds database for kraken2/bracken analysis
#
if INIT_KRAKEN2:

    rule kraken2_db:
        """
        Builds kraken2/bracken database.
        """
        output:
            db = INIT_KRAKEN2
        params:
            db_type = config['init']['kraken2_db_type'],
            read_len = config['init']['bracken_read_len'],
            std16_link = KRAKEN2_STANDARD16_LINK,
            std16_filename = os.path.basename(KRAKEN2_STANDARD16_LINK),
            plusPFP16_link = KRAKEN2_PLUSPFP16_LINK,
            plusPFP16_filename = os.path.basename(KRAKEN2_PLUSPFP16_LINK),
            outdir = 'kraken2_db/' + config['init']['kraken2_db_type']
        log: "logs/" + TIMESTAMP + "_kraken2_db.log"
        threads: INIT_KRAKEN2_THREADS
        # resources: cpus=4, mem_mb=4000, time_min=200
        run:
            shell(" echo 'Building kraken2/bracken database for '{params.db_type}'...' > {log} ")

            if params['db_type'] == 'standard':
                shell(
                    '''
                    kraken2-build \
                        --threads {threads} \
                        --use-ftp \
                        --standard \
                        --db {params.outdir} \
                        &>> {log}
                    '''
                )
            elif params['db_type'] in ['greengenes','rdp','silva']:
                shell(
                    '''
                    kraken2-build \
                        --threads {threads} \
                        --use-ftp \
                        --special {params.db_type} \
                        --db {params.outdir} \
                        &>> {log}
                    '''
                )
            elif params['db_type'] == 'nt':
                shell(
                    '''
                    kraken2-build \
                        --use-ftp \
                        --download-taxonomy \
                        --db {params.outdir} \
                        &>> {log}

                    kraken2-build \
                        --use-ftp \
                        --download-library {params.db_type} \
                        --db {params.outdir} \
                        &>> {log}

                    kraken2-build \
                        --threads {threads} \
                        --build \
                        --fast-build \
                        --db {params.outdir} \
                        &>> {log}
                    '''
                )
            elif params['db_type'] == 'std16_prebuilt':
                shell(
                    '''
                    curl -L {params.std16_link} --remote-name --silent
                    tar -zxf {params.std16_filename} -C {params.outdir}
                    rm {params.std16_filename}

                    echo 'Finished kraken2/bracken database' >> {log}
                    '''
                )
            elif params['db_type'] == 'plusPFP16_prebuilt':
                shell(
                    '''
                    curl -L {params.plusPFP16_link} --remote-name --silent
                    tar -zxf {params.plusPFP16_filename} -C {params.outdir}
                    rm {params.plusPFP16_filename}

                    echo 'Finished kraken2/bracken database' >> {log}
                    '''
                )

            if params['db_type'] not in ['std16_prebuilt', 'plusPFP16_prebuilt']:
                shell(
                    '''
                    bracken-build \
                        -t {threads} \
                        -k 35 \
                        -l {params.read_len} \
                        -d {params.outdir} \
                        &>> {log}

                    echo 'Finished kraken2/bracken database' >> {log}
                    '''
                )


onsuccess:
    shell(" echo Initialization successful! ")
onerror:
    shell(" echo Initialization halted! ")
