# rule centrifuge_get_db:
#     input:
#         seq = os.path.join(DBPATH,"common/ref-seqs.fna"),
#         tax = os.path.join(DBPATH,"common/ref-taxonomy.txt"),
#     output:
#         name_table = os.path.join(DBPATH,"centrifuge/taxonomy/names.dmp"),
#         tax_tree = os.path.join(DBPATH,"centrifuge/taxonomy/nodes.dmp"),
#         tax_map = os.path.join(DBPATH,"centrifuge/ref-tax.map"),
#         ref_seqs = os.path.join(DBPATH,"centrifuge/ref-seqs.fna"),
#         ref_tax = os.path.join(DBPATH,"centrifuge/ref-taxonomy.txt"),
#         tax_folder = directory(os.path.join(DBPATH,"centrifuge/taxonomy"))
#     threads: 1
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * config["centrifuge"]["dbmemory"]
#     params:
#         map_url = config["centrifuge"]["taxmapurl"]
#     conda:
#         os.path.join(ENVDIR,config["centrifuge"]["environment"])
#     log:
#         "logs/centrifuge_get_db.log"
#     benchmark:
#         "benchmarks/centrifuge_get_db.txt"
#     shell:
#         """
#         centrifuge-download -o {output.tax_folder} taxonomy > {log} 2>&1
        
#         wget {params.map_url} -q -O - | gzip -d -c - | \
#             awk '{{print $1\".\"$2\".\"$3\"\t\"$(NF)}}' \
#             > {output.tax_map} 2>> {log}
        
#         {SRCDIR}/todb.py -s {input.seq} -t {input.tax} -m centrifuge \
#             -S {output.ref_seqs} -T {output.ref_tax} 2>> {log}
#         """

rule centrifuge_build_db:
    input:
        name_table = os.path.join(DBPATH,"common/taxonomy/names.dmp"),
        tax_tree = os.path.join(DBPATH,"common/taxonomy/nodes.dmp"),
        conversion_table = os.path.join(DBPATH,"common/seqid2taxid.map"),
        ref_seqs = os.path.join(DBPATH,"common/ref-seqs.fna")
    output:
        touch(os.path.join(DBPATH,"centrifuge/CENTRIFUGE_DB_BUILD"))
    threads:
        config["centrifuge"]["dbthreads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["centrifuge"]["dbmemory"]
    params:
        prefix = os.path.join(DBPATH,"centrifuge/ref-db")
    conda:
        os.path.join(ENVDIR,config["centrifuge"]["environment"])
    log:
        "logs/centrifuge_build_db.log"
    benchmark:
        "benchmarks/centrifuge_build_db.txt"
    shell:
        """
        centrifuge-build \
          --threads {threads} \
          --conversion-table {input.conversion_table} \
          --taxonomy-tree {input.tax_tree} \
          --name-table {input.name_table} \
          {input.ref_seqs} {params.prefix} > {log} 2>&1
        """


rule centrifuge_classify:
    input:
        rules.centrifuge_build_db.output,
        fastq = get_seqfiletype,
        ref_seqs = os.path.join(DBPATH,"centrifuge/ref-seqs.fna")
        #ref_seqs = "db/kraken/data/SILVA_132_SSURef_Nr99_tax_silva.fasta"
    output:
        report = temp("classifications/{run}/centrifuge/{sample}.report.tsv"),
        classification = "classifications/{run}/centrifuge/{sample}.centrifuge.out",
    threads:
        config["centrifuge"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["centrifuge"]["memory"]
    params:
        index_prefix = os.path.join(DBPATH,"centrifuge/ref-db")
    conda:
        os.path.join(ENVDIR,config["centrifuge"]["environment"])
    log:
        "logs/{run}/centrifuge_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/centrifuge_classify_{sample}.txt"
    shell:
        """
        centrifuge -x {params.index_prefix} \
          -U {input.fastq} \
          --threads {threads} \
          --report-file {output.report} \
          -S  {output.classification} \
          --met-stderr > {log} 2>&1
        """


rule centrifuge_tomat:
    input:
        out = "classifications/{run}/centrifuge/{sample}.centrifuge.out",
        ref_seqs = os.path.join(DBPATH,"centrifuge/ref-seqs.fna")
    output:
        taxlist = "classifications/{run}/centrifuge/{sample}.centrifuge.taxlist",
        taxmat = "classifications/{run}/centrifuge/{sample}.centrifuge.taxmat",
        otumat = "classifications/{run}/centrifuge/{sample}.centrifuge.otumat"
    threads: 1
    conda:
        os.path.join(ENVDIR,config["centrifuge"]["environment"])
    log:
        "logs/{run}/centrifuge_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/centrifuge_tomat_{sample}.txt"
    shell:
        "{SRCDIR}/tomat.py -c {input.out} -f {input.ref_seqs} 2> {log}"
