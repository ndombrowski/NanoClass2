def get_fastq(wildcards):
    return smpls.loc[(wildcards.run, wildcards.sample), ["path"]].dropna()

def get_seqfiletype(wildcards):
    if config["subsample"]["skip"] is True:
        return "data/{run}/chopper/{sample}.filtered.fastq.gz"
    else:
        return "data/{run}/chopper/{sample}.subsampled.fastq.gz"


rule prep_porechop:
    input:
        get_fastq
    output:
        "data/{run}/porechopped/{sample}.trimmed.fastq.gz"
    threads:
        config["porechop"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["porechop"]["memory"]
    priority: 50
    conda:
        os.path.join(ENVDIR,config["porechop"]["environment"])
    params:
        check_reads = config["porechop"]["checkreads"]
    log:
        "logs/{run}/prep_porechop_{sample}.log"
    benchmark:
        "benchmarks/{run}/prep_porechop_{sample}.txt"
    shell:
        """
        porechop --input {input} \
          --output {output} \
          --threads {threads} \
          --check_reads {params.check_reads} \
          --discard_middle > {log}
        """

rule run_chopper:
    input:
        "data/{run}/porechopped/{sample}.trimmed.fastq.gz"
    output:
        "data/{run}/chopper/{sample}.filtered.fastq.gz"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["qual_filtering"]["memory"]
    priority: 50
    params:
        threads = config["qual_filtering"]["threads"],
        min_len = config["qual_filtering"]["minlen"],
        max_len = config["qual_filtering"]["maxlen"],
        quality = config["qual_filtering"]["quality"],
        headcrop = config["qual_filtering"]["headcrop"],
        tailcrop = config["qual_filtering"]["tailcrop"]
    log:
        "logs/{run}/prep_chopper_{sample}.log"
    benchmark:
        "benchmarks/{run}/prep_chopper_{sample}.txt"
    conda:
        os.path.join(ENVDIR,config["qual_filtering"]["environment"])
    shell:
        """
        gunzip -c {input} |\
            chopper -q {params.quality} \
            --headcrop {params.headcrop} --tailcrop {params.tailcrop}   \
            -l {params.min_len} --maxlength {params.max_len}  --threads {params.threads} |\
            gzip > {output}  2> {log}
        """


rule prep_subsample:
    input:
        fastq = "data/{run}/chopper/{sample}.filtered.fastq.gz",
    output:
        "data/{run}/chopper/{sample}.subsampled.fastq.gz"
    threads: 1
    priority: 50
    params:
        n = config["subsample"]["samplesize"],
        seed = 12345
    conda:
        os.path.join(ENVDIR,config["subsample"]["environment"])
    log:
        "logs/{run}/prep_subsample_{sample}.log"
    benchmark:
        "benchmarks/{run}/prep_subsample_{sample}.txt"
    shell:
        """
        seqtk sample -s {params.seed} {input} {params.n} | \
            gzip -c > {output}
        """


rule prep_fasta_query:
    input:
        get_seqfiletype
    output:
        "data/{run}/chopper/{sample}.fasta"
    threads: 1
    priority: 50
    conda:
        os.path.join(ENVDIR,config["subsample"]["environment"])
    shell:
        "zcat < {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}"



rule prep_chopper_plot:
    input:
        "data/{run}/chopper/{sample}.filtered.fastq.gz"
    output:
        report("plots/{run}/chopper/{sample}.filtered.pdf", 
               caption="../report/fig-chopper.rst", 
               category="Read-processing"
              )
    log:
        "logs/{run}/prep_chopper_plot_{sample}.log"
    benchmark:
        "benchmarks/{run}/prep_chopper_plot_{sample}.txt"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["nanoplot"]["memory"]
    params:
        downsample = config["nanoplot"]["downsample"]
    conda:
        os.path.join(ENVDIR,config["subsample"]["environment"])
    shell:
        """
        pistis --fastq {input} --output {output} \
          --downsample {params.downsample} 2> {log}
        """


rule prep_chopper_stats:
    input:
        "data/{run}/chopper/{sample}.filtered.fastq.gz"
    output:
        "stats/{run}/chopper/{sample}.filtered.txt"
    log:
        "logs/{run}/prep_chopper_stats_{sample}.log"
    benchmark:
        "benchmarks/{run}/prep_chopper_stats_{sample}.txt"
    threads:
        config["nanostats"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["nanostats"]["memory"]
    conda:
        os.path.join(ENVDIR,config["subsample"]["environment"])
    shell:
        """
        NanoStat --fastq {input} --name {output} \
          --threads {threads} 2> {log}
        """
