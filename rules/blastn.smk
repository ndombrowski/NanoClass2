rule blastn_build_db:
    input:
        "db/common/ref-seqs.fna"
    output:
        touch("db/blastn/DB_BUILD")
    threads: 1
    log:
        "logs/blastn_build_db.log"
    benchmark:
        "benchmarks/blastn_build_db.txt"
    conda:
        config["blastn"]["environment"]
    shell:
        """
        makeblastdb -in {input} -parse_seqids \
            -dbtype nucl 2>&1 | tee -a {log}
        """


rule blastn_chunk:
    input:
        fasta = rules.prep_fasta_query.output
    output:
        temp(expand("classifications/{{run}}/blastn/{{sample}}/{{sample}}.part-{chunk}.fasta",
            chunk = range(1, (config["blastn"]["threads"]+1))
        ))
    params:
        n_chunks = config["blastn"]["threads"],
        out_dir = "classifications/{run}/blastn/{sample}"
    conda:
        config["blastn"]["environment"]
    shell:
        """
        fasta-splitter --n-parts {params.n_chunks} {input} \
            --nopad --out-dir {params.out_dir}
        """


rule blastn_classify:
    input:
        rules.blastn_build_db.output,
        db = "db/common/ref-seqs.fna",
        fasta = "classifications/{run}/blastn/{sample}/{sample}.part-{chunk}.fasta"
    output:
        all = temp("classifications/{run}/blastn/{sample}/{chunk}.blastn.out")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["blastn"]["memory"]
    params:
        eval = config["blastn"]["evalue"],
        pident = config["blastn"]["pctidentity"],
        nseqs = config["blastn"]["ntargetseqs"]
    conda:
        config["blastn"]["environment"]
    log:
        "logs/{run}/blastn_classify_{sample}_{chunk}.log"
    benchmark:
        "benchmarks/{run}/blastn_classify_{sample}_{chunk}.txt"
    shell:
        """
        blastn -task 'blastn' -db {input.db} -evalue {params.eval} \
            -perc_identity {params.pident} -max_target_seqs {params.nseqs} \
            -query {input.fasta} -out {output.all} -outfmt 6 \
            2> {log}
        """

rule blastn_aggregate:
    input:
        out = expand("classifications/{{run}}/blastn/{{sample}}/{chunk}.blastn.out",
            chunk = range(1, (config["blastn"]["threads"]+1))
        ),
        benchm = expand("benchmarks/{{run}}/blastn_classify_{{sample}}_{chunk}.txt",
            chunk = range(1, (config["blastn"]["threads"]+1))
        )
    output:
        out = temp("classifications/{run}/blastn/{sample}.blastn.out"),
        benchm = "benchmarks/{run}/blastn_classify_{sample}.txt"
    threads: 1
    params:
        threads = config["blastn"]["threads"],
        alnlen = config["blastn"]["alnlength"]
    conda:
        config["blastn"]["environment"]
    log:
        "logs/{run}/blastn_aggregate_{sample}.log"
    benchmark:
        "benchmarks/{run}/blastn_aggregate_{sample}.txt"
    shell:
        """
        cat {input.out} | awk -F '\\t' -v l={params.alnlen} \
            '{{if ($4 > l) print $0}}' \
            > {output.out}
        cat {input.benchm} | grep -P "^[0-9]" | \
          awk '{{sum+=$1}} END {{print "s"; print sum/{params.threads}}}' \
          > {output.benchm} 2> {log}
        """

rule blast_tolca:
    input:
        blast = "classifications/{run}/blastn/{sample}.blastn.out",
        db = "db/common/ref-taxonomy.txt"
    output:
        "classifications/{run}/blastn/{sample}.blastn.taxlist"
    threads: 1
    params:
        lcacons = config["blastn"]["lcaconsensus"]
    conda:
        config["blastn"]["environment"]
    log:
        "logs/{run}/blastn_tolca_{sample}.log"
    benchmark:
        "benchmarks/{run}/blastn_tolca_{sample}.txt"
    shell:
        """
        scripts/tolca.py -b {input.blast} -t {input.db} \
             -l {output} -c {params.lcacons} > {log}
        """

rule blastn_tomat:
    input:
        "classifications/{run}/blastn/{sample}.blastn.taxlist"
    output:
        taxmat = "classifications/{run}/blastn/{sample}.blastn.taxmat",
        otumat = "classifications/{run}/blastn/{sample}.blastn.otumat"
    threads: 1
    conda:
        config["blastn"]["environment"]
    log:
        "logs/{run}/blastn_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/blastn_tomat_{sample}.txt"
    shell:
        "scripts/tomat.py -l {input} 2> {log}"
