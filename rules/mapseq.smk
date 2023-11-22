rule mapseq_classify:
    input:
        query = rules.prep_fasta_query.output,
        ref_tax = "db/common/ref-taxonomy.txt",
        ref_seqs = "db/common/ref-seqs.fna"
    output:
        temp("classifications/{run}/mapseq/{sample}.mapseq.out")
    threads:
        config["mapseq"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["mapseq"]["memory"]
    singularity:
        config["mapseq"]["container"]
    log:
        "logs/{run}/mapseq_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/mapseq_classify_{sample}.txt"
    shell:
        """
        mapseq -nthreads {threads} {input.query} {input.ref_seqs} \
           {input.ref_tax}  > {output} 2> {log}
        """

rule mapseq_tomat:
    input:
        out = "classifications/{run}/mapseq/{sample}.mapseq.out",
        db = "db/common/ref-taxonomy.txt"
    output:
        taxlist = "classifications/{run}/mapseq/{sample}.mapseq.taxlist",
        taxmat = "classifications/{run}/mapseq/{sample}.mapseq.taxmat",
        otumat = "classifications/{run}/mapseq/{sample}.mapseq.otumat"
    threads: 1
    log:
        "logs/{run}/mapseq_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/mapseq_tomat_{sample}.txt"
    shell:
        "scripts/tomat.py -b {input.out} -t {input.db} 2> {log}"

