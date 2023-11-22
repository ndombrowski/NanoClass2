rule rdp_build_db:
    input:
        seq = os.path.join(DBPATH,"common/ref-seqs.fna"),
        tax = os.path.join(DBPATH,"common/ref-taxonomy.txt")
    output:
        seq = os.path.join(DBPATH,"rdp/ref-seqs.fna.gz"),
        tax = os.path.join(DBPATH,"rdp/ref-taxonomy.txt")
    threads: 1
    log:
        "logs/rdp_db.log"
    benchmark:
        "benchmarks/rdp_db.txt"
    conda:
        os.path.join(ENVDIR,config["rdp"]["environment"])
    shell:
        """
        {SRCDIR}/todb.py -s {input.seq} -t {input.tax} -m rdp -S tmp.seq -T {output.tax} 2> {log}
        gzip -c tmp.seq > {output.seq} && rm tmp.seq 2>> {log}
        """


rule rdp_classify:
    input:
        db = os.path.join(DBPATH,"rdp/ref-seqs.fna.gz"),
        query = rules.prep_fasta_query.output
    output:
        "classifications/{run}/rdp/{sample}.rdp.taxlist"
    params:
        config["rdp"]["pctthreshold"]
    threads:
        config["rdp"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["rdp"]["memory"]
    conda:
        os.path.join(ENVDIR,config["rdp"]["environment"])
    log:
        "logs/{run}/rdp_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/rdp_classify_{sample}.txt"
    shell:
        "Rscript {SRCDIR}/assigntaxonomy.R {input.db} {input.query} {output} {threads} {params} 2> {log}"


rule rdp_tomat:
    input:
        list = "classifications/{run}/rdp/{sample}.rdp.taxlist",
    output:
        taxmat = "classifications/{run}/rdp/{sample}.rdp.taxmat",
        otumat = "classifications/{run}/rdp/{sample}.rdp.otumat"
    threads: 1
    conda:
        os.path.join(ENVDIR,config["rdp"]["environment"])
    log:
        "logs/{run}/rdp_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/rdp_tomat_{sample}.txt"
    shell:
        "{SRCDIR}/tomat.py -l {input.list} 2> {log}"
