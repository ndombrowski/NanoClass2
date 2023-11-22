rule spingo_build_db:
    input:
        seq = "db/common/ref-seqs.fna",
        tax = "db/common/ref-taxonomy.txt"
    output:
        seq = "db/spingo/ref-seqs.fna",
        tax = "db/spingo/ref-taxonomy.txt"
    threads:
        config["spingo"]["dbthreads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["spingo"]["dbmemory"]
    log:
        "logs/spingo_db.log"
    benchmark:
        "benchmarks/spingo_db.txt"
    conda:
        config["spingo"]["environment"]
    shell:
        """
        scripts/todb.py -s {input.seq} -t {input.tax} -m spingo \
            -S {output.seq} -T {output.tax} 2> {log}
        spindex -k 8 -p {threads} -d {output.seq} 2>> {log}
        """


rule spingo_classify:
    input:
        db = "db/spingo/ref-seqs.fna",
        fasta = rules.prep_fasta_query.output
    output:
        out = "classifications/{run}/spingo/{sample}.spingo.taxlist",
        tmp = temp("classifications/{run}/spingo/{sample}.spingo.tmp")
    threads:
        config["spingo"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["spingo"]["memory"]
    conda:
        config["spingo"]["environment"]
    log:
        "logs/{run}/spingo_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/spingo_classify_{sample}.txt"
    shell:
        """
        spingo -d {input.db} -k 8 -a \
          -p {threads} -i {input.fasta} \
          > {output.tmp} 2> {log}
        sed 's/ /\\t/' {output.tmp} | \
        awk -F '\\t' -v OFS='\\t' '{{
            if (NR == 1)
                print "#readid","Domain","Phylum","Class","Order","Family","Genus";
            else
                print $1, $4, $6, $8, $10, $12, $14
        }}' > {output.out} 2> {log}
        """


rule spingo_tomat:
    input:
        list = "classifications/{run}/spingo/{sample}.spingo.taxlist",
    output:
        taxmat = "classifications/{run}/spingo/{sample}.spingo.taxmat",
        otumat = "classifications/{run}/spingo/{sample}.spingo.otumat"
    threads: 1
    conda:
        config["spingo"]["environment"]
    log:
        "logs/{run}/spingo_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/spingo_tomat_{sample}.txt"
    shell:
        "scripts/tomat.py -l {input.list} 2> {log}"
