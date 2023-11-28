rule kraken_build_db:
    input:
        rules.common_download_db.output
    output:
        #name_table = os.path.join(DBPATH,"kraken/taxonomy/names.dmp"),
        #tax_tree = os.path.join(DBPATH,"kraken/taxonomy/nodes.dmp"),
        #conversion_table = os.path.join(DBPATH,"kraken/seqid2taxid.map"),
        #ref_seqs = os.path.join(DBPATH,"kraken/data/SILVA_138.1_SSURef_NR99_tax_silva.fasta")
        out = os.path.join(DBPATH, "kraken/kraken_done")
    threads:
        config["kraken"]["dbthreads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["kraken"]["dbmemory"]
    params:
        db = os.path.join(DBPATH,"common"),
        db_type = config["kraken"]["dbtype"]
    conda:
        os.path.join(ENVDIR,config["kraken"]["environment"])
    log:
        "logs/kraken_build_db.log"
    benchmark:
        "benchmarks/kraken_build_db.txt"
    shell:
        """
        # --db {params.db} --special {params.db_type} \
        #  --threads {threads} > {log} 2>&1

        kraken2-build --build --fast-build --threads {threads} --db {params.db}

        touch {output}
        """

rule kraken_classify:
    input:
        rules.kraken_build_db.output,
        fastq = get_seqfiletype
    output:
        report = temp("classifications/{run}/kraken/{sample}.kraken.report"),
        out = temp("classifications/{run}/kraken/{sample}.kraken.out")
    threads:
        config["kraken"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["kraken"]["memory"]
    params:
        db_dir = os.path.join(DBPATH,"common"),
        confidence_score = config["kraken"]["confidence_score"]
    conda:
        os.path.join(ENVDIR,config["kraken"]["environment"])
    log:
        "logs/{run}/kraken_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/kraken_classify_{sample}.txt"
    shell:
        """
        kraken2 --db {params.db_dir} \
            --confidence {params.confidence_score} \
            --output {output.out} \
            --report {output.report} \
            --gzip-compressed \
            --threads {threads} {input.fastq} 2> {log}

        """


rule kraken_tomat:
    input:
        kraken_out = "classifications/{run}/kraken/{sample}.kraken.out",
        silva_seqs = os.path.join(DBPATH,"common/data/SILVA_138.1_SSURef_NR99_tax_silva.fasta"),
        kraken_map = os.path.join(DBPATH,"common/seqid2taxid.map")
    output:
        taxlist = "classifications/{run}/kraken/{sample}.kraken.taxlist",
        taxmat = "classifications/{run}/kraken/{sample}.kraken.taxmat",
        otumat = "classifications/{run}/kraken/{sample}.kraken.otumat"
    threads: 1
    conda:
        os.path.join(ENVDIR,config["kraken"]["environment"])
    log:
        "logs/{run}/kraken_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/kraken_tomat_{sample}.txt"
    shell:
        """
        python3 {SRCDIR}/tomat.py -k {input.kraken_out} -f {input.silva_seqs} \
          -m {input.kraken_map} 2> {log}
        """

