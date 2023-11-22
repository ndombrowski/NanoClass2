rule dcmegablast_build_db:
    input:
        "db/blastn/DB_BUILD"
    output:
        touch("db/dcmegablast/DB_BUILD")


rule dcmegablast_chunk:
    input:
        fasta = rules.prep_fasta_query.output
    output:
        temp(expand("classifications/{{run}}/dcmegablast/{{sample}}/{{sample}}.part-{chunk}.fasta",
            chunk = range(1, (config["dcmegablast"]["threads"]+1))
        ))
    params:
        n_chunks = config["dcmegablast"]["threads"],
        out_dir = "classifications/{run}/dcmegablast/{sample}"
    conda:
        config["dcmegablast"]["environment"]
    shell:
        """
        fasta-splitter --n-parts {params.n_chunks} {input} \
            --nopad --out-dir {params.out_dir}
        """

rule dcmegablast_classify:
    input:
        rules.dcmegablast_build_db.output,
        db = "db/common/ref-seqs.fna",
        fasta = "classifications/{run}/dcmegablast/{sample}/{sample}.part-{chunk}.fasta"
    output:
        all = temp("classifications/{run}/dcmegablast/{sample}/{chunk}.dcmegablast.out")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["dcmegablast"]["memory"]
    params:
       	eval = config["blastn"]["evalue"],
        pident = config["dcmegablast"]["pctidentity"],
        nseqs = config["blastn"]["ntargetseqs"]
    conda:
        config["dcmegablast"]["environment"]
    log:
        "logs/{run}/dcmegablast_classify_{sample}_{chunk}.log"
    benchmark:
        "benchmarks/{run}/dcmegablast_classify_{sample}_{chunk}.txt"
    shell:
        """
        blastn -task 'dc-megablast' -db {input.db} -evalue {params.eval} \
            -perc_identity {params.pident} -max_target_seqs {params.nseqs} \
            -query {input.fasta} -out {output.all} -outfmt 6 \
            2> {log} 
        """
   

rule dcmegablast_aggregate:
    input:
        out = expand("classifications/{{run}}/dcmegablast/{{sample}}/{chunk}.dcmegablast.out",
            chunk = range(1, (config["dcmegablast"]["threads"]+1))
        ),
        benchm = expand("benchmarks/{{run}}/dcmegablast_classify_{{sample}}_{chunk}.txt",
            chunk = range(1, (config["dcmegablast"]["threads"]+1))
        )
    output:
        out = temp("classifications/{run}/dcmegablast/{sample}.dcmegablast.out"),
        benchm = "benchmarks/{run}/dcmegablast_classify_{sample}.txt"
    threads: 1
    params:
        threads = config["dcmegablast"]["threads"],
        alnlen = config["dcmegablast"]["alnlength"]
    conda:
        config["dcmegablast"]["environment"]
    log:
        "logs/{run}/dcmegablast_aggregate_{sample}.log"
    benchmark:
        "benchmarks/{run}/dcmegablast_aggregate_{sample}.txt"
    shell:
        """
	cat {input.out} | awk -F '\\t' -v l={params.alnlen} \
            '{{if ($4 > l) print $0}}' \
            > {output.out}
        cat {input.benchm} | grep -P "^[0-9]" | \
          awk '{{sum+=$1}} END {{print "s"; print sum / {params.threads}}}' \
          > {output.benchm} 2> {log}
        """


rule dcmegablast_tolca:
    input:
        blast = "classifications/{run}/dcmegablast/{sample}.dcmegablast.out",
        db = "db/common/ref-taxonomy.txt"
    output:
        "classifications/{run}/dcmegablast/{sample}.dcmegablast.taxlist",
    threads: 1
    params:
        lcacons = config["dcmegablast"]["lcaconsensus"]
    conda:
        config["dcmegablast"]["environment"]
    log:
        "logs/{run}/dcmegablast_tolca_{sample}.log"
    benchmark:
        "benchmarks/{run}/dcmegablast_tolca_{sample}.txt"
    shell:
      	"""
        scripts/tolca.py -b {input.blast} -t {input.db} \
             -l {output} -c {params.lcacons} > {log}
        """

rule dcmegablast_tomat:
    input:
        "classifications/{run}/dcmegablast/{sample}.dcmegablast.taxlist"
    output: 
        taxmat = "classifications/{run}/dcmegablast/{sample}.dcmegablast.taxmat",
        otumat = "classifications/{run}/dcmegablast/{sample}.dcmegablast.otumat"
    threads: 1
    conda:
        config["dcmegablast"]["environment"]
    log:
        "logs/{run}/dcmegablast_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/dcmegablast_tomat_{sample}.txt"
    shell:
        "scripts/tomat.py -l {input} 2> {log}"
