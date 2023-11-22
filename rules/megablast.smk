rule megablast_build_db:
    input:
        "db/blastn/DB_BUILD"
    output:
        touch("db/megablast/DB_BUILD")


rule megablast_chunk:
    input:
        fasta = rules.prep_fasta_query.output
    output:
        temp(expand("classifications/{{run}}/megablast/{{sample}}/{{sample}}.part-{chunk}.fasta",
            chunk = range(1, (config["megablast"]["threads"]+1))
        ))
    params:
        n_chunks = config["megablast"]["threads"],
        out_dir = "classifications/{run}/megablast/{sample}"
    conda:
        config["megablast"]["environment"]
    shell:
        """
        fasta-splitter --n-parts {params.n_chunks} {input} \
            --nopad --out-dir {params.out_dir}
        """

rule megablast_classify:
    input:
        rules.megablast_build_db.output,
        db = "db/common/ref-seqs.fna",
        fasta = "classifications/{run}/megablast/{sample}/{sample}.part-{chunk}.fasta"
    output:
        all = temp("classifications/{run}/megablast/{sample}/{chunk}.megablast.out")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["megablast"]["memory"]
    params:
       	eval = config["blastn"]["evalue"],
        pident = config["megablast"]["pctidentity"],
        nseqs = config["blastn"]["ntargetseqs"]
    conda:
        config["megablast"]["environment"]
    log:
        "logs/{run}/megablast_classify_{sample}_{chunk}.log"
    benchmark:
        "benchmarks/{run}/megablast_classify_{sample}_{chunk}.txt"
    shell:
        """
        blastn -task 'megablast' -db {input.db} -evalue {params.eval} \
            -perc_identity {params.pident} -max_target_seqs {params.nseqs} \
            -query {input.fasta} -out {output.all} -outfmt 6 \
            2> {log} 
        """
   

rule megablast_aggregate:
    input:
        out = expand("classifications/{{run}}/megablast/{{sample}}/{chunk}.megablast.out",
            chunk = range(1, (config["megablast"]["threads"]+1))
        ),
        benchm = expand("benchmarks/{{run}}/megablast_classify_{{sample}}_{chunk}.txt",
            chunk = range(1, (config["megablast"]["threads"]+1))
        )
    output:
        out = temp("classifications/{run}/megablast/{sample}.megablast.out"),
        benchm = "benchmarks/{run}/megablast_classify_{sample}.txt"
    threads: 1
    params:
        threads = config["megablast"]["threads"],
        alnlen = config["megablast"]["alnlength"]
    conda:
        config["megablast"]["environment"]
    log:
        "logs/{run}/megablast_aggregate_{sample}.log"
    benchmark:
        "benchmarks/{run}/megablast_aggregate_{sample}.txt"
    shell:
        """
	cat {input.out} | awk -F '\\t' -v l={params.alnlen} \
            '{{if ($4 > l) print $0}}' \
            > {output.out}
        cat {input.benchm} | grep -P "^[0-9]" | \
          awk '{{sum+=$1}} END {{print "s"; print sum / {params.threads}}}' \
          > {output.benchm} 2> {log}
        """


rule megablast_tolca:
    input:
        blast = "classifications/{run}/megablast/{sample}.megablast.out",
        db = "db/common/ref-taxonomy.txt"
    output:
        "classifications/{run}/megablast/{sample}.megablast.taxlist",
    threads: 1
    params:
        lcacons = config["megablast"]["lcaconsensus"]
    conda:
        config["megablast"]["environment"]
    log:
        "logs/{run}/megablast_tolca_{sample}.log"
    benchmark:
        "benchmarks/{run}/megablast_tolca_{sample}.txt"
    shell:
      	"""
        scripts/tolca.py -b {input.blast} -t {input.db} \
             -l {output} -c {params.lcacons} > {log}
        """

rule megablast_tomat:
    input:
        "classifications/{run}/megablast/{sample}.megablast.taxlist"
    output: 
        taxmat = "classifications/{run}/megablast/{sample}.megablast.taxmat",
        otumat = "classifications/{run}/megablast/{sample}.megablast.otumat"
    threads: 1
    conda:
        config["megablast"]["environment"]
    log:
        "logs/{run}/megablast_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/megablast_tomat_{sample}.txt"
    shell:
        "scripts/tomat.py -l {input} 2> {log}"
