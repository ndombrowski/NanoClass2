rule idtaxa_build_db:
    input:
        ref_seqs = "db/common/ref-seqs.fna",
        ref_tax = "db/common/ref-taxonomy.txt"
    output:
        "db/idtaxa/ref-db.Rdata"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["idtaxa"]["dbmemory"]
    conda:
        config["idtaxa"]["environment"]
    log:
        "logs/idtaxa_learn_taxa.log"
    benchmark:
        "benchmarks/idtaxa_learn_taxa.txt"
    shell:
        "Rscript scripts/learntaxa.R {input.ref_seqs} {input.ref_tax} {output} 2> {log}"


rule idtaxa_classify:
    input:
        db = "db/idtaxa/ref-db.Rdata",
        query = rules.prep_fasta_query.output
    output:
        tmp = temp("classifications/{run}/idtaxa/{sample}.idtaxa.tmp"),
        out = "classifications/{run}/idtaxa/{sample}.idtaxa.taxlist"
    threads:
        config["idtaxa"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["idtaxa"]["memory"]
    params:
        config["idtaxa"]["pctthreshold"]
    conda:
        config["idtaxa"]["environment"]
    log:
        "logs/{run}/idtaxa_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/idtaxa_classify_{sample}.txt"
    shell:
        """
        Rscript scripts/idtaxa.R {input.db} {input.query} {output.tmp} {threads} {params} 2> {log}
        awk 'BEGIN{{FS=OFS="\\t"}} {{for (i=1; i<=7; i++) if ($i ~ /^ *$/) $i="NA"}} 1' {output.tmp} > {output.out}
        """


rule idtaxa_tomat:
    input:
        list = "classifications/{run}/idtaxa/{sample}.idtaxa.taxlist",
    output:
        taxmat = "classifications/{run}/idtaxa/{sample}.idtaxa.taxmat",
        otumat = "classifications/{run}/idtaxa/{sample}.idtaxa.otumat"
    threads: 1
    conda:
        config["idtaxa"]["environment"]
    log:
        "logs/{run}/idtaxa_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/idtaxa_tomat_{sample}.txt"
    shell:
        "scripts/tomat.py -l {input.list} 2> {log}"

