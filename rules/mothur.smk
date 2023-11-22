rule mothur_build_db:
    input:
        aln = "db/common/ref-seqs.aln",
        tax = "db/common/ref-taxonomy.txt"
    output:
        aln = "db/mothur/ref-seqs.aln",
        tax = "db/mothur/ref-taxonomy.txt"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["mothur"]["dbmemory"]
    log:
        "logs/mothur_build_db.log"
    benchmark:
        "benchmarks/mothur_build_db.txt"
    conda:
        config["mothur"]["environment"]
    shell:
        """
        scripts/todb.py -s {input.aln} -t {input.tax} -m mothur \
            -S {output.aln} -T {output.tax} 2> {log}
        """

rule mothur_classify:
    input:
        query = rules.prep_fasta_query.output,
        aln = "db/mothur/ref-seqs.aln"
    output:
        dir = temp(directory("classifications/{run}/mothur/{sample}/")),
        out = temp("classifications/{run}/mothur/{sample}.mothur.out")
    params:
        file = "{sample}.align.report"
    threads:
        config["mothur"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["mothur"]["memory"]
    conda:
        config["mothur"]["environment"]
    log:
        "logs/{run}/mothur_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/mothur_classify_{sample}.txt"
    shell:
        """
        mkdir -p {output.dir}
        cp {input} {output.dir}
        queryf=$(basename {input.query})
        alnf=$(basename {input.aln})
        echo \"align.seqs(candidate={output.dir}/$queryf, template={output.dir}/$alnf, \
            processors={threads}, ksize=6, align=needleman); -q\" > {output.dir}/mothur.cmd
        mothur {output.dir}/mothur.cmd 2> {log}
        awk -F '\\t' -v OFS='\\t' '{{if (NR==1) printf "%s","#"; print $1, $3}}' \
          {output.dir}/{params.file} > {output.out} 2>> {log} 
        """


rule mothur_tomat:
    input:
        tax = "db/mothur/ref-taxonomy.txt",
        out = "classifications/{run}/mothur/{sample}.mothur.out"
    output:
        taxlist = "classifications/{run}/mothur/{sample}.mothur.taxlist",
        taxmat = "classifications/{run}/mothur/{sample}.mothur.taxmat",
        otumat = "classifications/{run}/mothur/{sample}.mothur.otumat"
    threads: 1
    log:
        "logs/{run}/mothur_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/mothur_tomat_{sample}.txt"
    conda:
        config["mothur"]["environment"]
    shell:
        "scripts/tomat.py -b {input.out} -t {input.tax} 2> {log}"
