rule minimap_classify:
    input:
        target = os.path.join(DBPATH,"common/ref-seqs.fna"),
        query = get_seqfiletype
    output:
        temp("classifications/{run}/minimap/{sample}.minimap.bam")
    threads: 
        config["minimap"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["minimap"]["memory"]
    conda:
        os.path.join(ENVDIR,config["minimap"]["environment"])
    params:
        extra = "-K 25M --no-kalloc --print-qname -aLx map-ont",
        nseqs = config["minimap"]["ntargetseqs"]
    log:
        "logs/{run}/minimap_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/minimap_classify_{sample}.txt"
    shell:
        """
        minimap2 {params.extra} -t {threads} -N {params.nseqs} \
          {input.target} {input.query} -o {output} 2> {log}
        """

rule minimap_bam2out:
    input:
        "classifications/{run}/minimap/{sample}.minimap.bam"
    output:
        "classifications/{run}/minimap/{sample}.minimap.out"
    threads:
        config["minimap"]["threads"]
    conda:
        os.path.join(ENVDIR,config["minimap"]["environment"])
    log:
        "logs/{run}/minimap_sortbam_{sample}.log"
    benchmark:
        "benchmarks/{run}/minimap_sortbam_{sample}.txt"
    shell:
        """
        samtools sort -@ {threads} {input} | \
        samtools view -@ {threads} | \
        cut -f 1,3 > {output} 2> {log}
        """

rule minimap_tolca:
    input:
        mm = "classifications/{run}/minimap/{sample}.minimap.out",
        db = os.path.join(DBPATH,"common/ref-taxonomy.txt")
    output:
        "classifications/{run}/minimap/{sample}.minimap.taxlist"
    threads: 1
    params:
        lcacons = config["minimap"]["lcaconsensus"]
    conda:
        os.path.join(ENVDIR,config["blastn"]["environment"])
    log:
        "logs/{run}/minimap_tolca_{sample}.log"
    benchmark:
        "benchmarks/{run}/minimap_tolca_{sample}.txt"
    shell:
        """
	{SRCDIR}/tolca.py -b {input.mm} -t {input.db} \
             -l {output} -c {params.lcacons} > {log}
        """


rule minimap_tomat:
    input:
        "classifications/{run}/minimap/{sample}.minimap.taxlist",
    output:
        taxmat = "classifications/{run}/minimap/{sample}.minimap.taxmat",
        otumat = "classifications/{run}/minimap/{sample}.minimap.otumat"
    threads: 1
    conda:
        os.path.join(ENVDIR,config["minimap"]["environment"])
    log:
        "logs/{run}/minimap_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/minimap_tomat_{sample}.txt"
    shell:
        "{SRCDIR}/tomat.py -l {input} 2> {log}"
