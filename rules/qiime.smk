rule qiime_import_query:
    input:
        fasta = rules.prep_fasta_query.output
    output:
        temp("classifications/{run}/qiime/{sample}.seqs.qza")
    threads:
        1
    conda:
        config["qiime"]["environment"]
    log:
        "logs/{run}/qiime_import_query_{sample}.log"
    benchmark:
        "benchmarks/{run}/qiime_import_query_{sample}.txt"
    shell:
        """
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {input.fasta} \
            --output-path {output} \
            2>&1 | tee -a {log}
        """

rule qiime_import_refseq:
    input:
        "db/common/ref-seqs.fna"
    output:
        "db/qiime/ref-seqs.qza"
    threads: 1
    conda:
        config["qiime"]["environment"]
    log:
        "logs/qiime_import_refseq.log"
    benchmark:
        "benchmarks/qiime_import_refseq.txt"
    shell:
        """
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {input} \
            --output-path {output} \
            2>&1 | tee -a {log}
        """

rule qiime_import_taxonomy:
    input:
        "db/common/ref-taxonomy.txt"
    output:
        "db/qiime/ref-taxonomy.qza"
    threads: 1
    conda:
        config["qiime"]["environment"]
    log:
        "logs/qiime_import_taxonomy.log"
    benchmark:   
        "benchmarks/qiime_import_taxonomy.txt"
    shell:
        """
        qiime tools import \
            --type FeatureData[Taxonomy] \
            --input-format HeaderlessTSVTaxonomyFormat \
            --input-path {input} \
            --output-path {output} \
            2>&1 | tee -a {log}
        """

rule qiime_classify:
    input:
        query="classifications/{run}/qiime/{sample}.seqs.qza",
        taxo="db/qiime/ref-taxonomy.qza",
        ref_seqs="db/qiime/ref-seqs.qza"
    output:
        qza = temp("classifications/{run}/qiime/{sample}.qiime.qza"),
        path = temp(directory("classifications/{run}/qiime/{sample}/")),
        out = "classifications/{run}/qiime/{sample}.qiime.out"
    params:
        pident = config["qiime"]["pctidentity"],
        nseqs = config["qiime"]["ntargetseqs"],
        lcacons = config["qiime"]["lcaconsensus"]
    threads:
        config["qiime"]["threads"]
    conda:
        config["qiime"]["environment"]
    log:
        "logs/{run}/qiime_classify_{sample}.log"
    benchmark:
        "benchmarks/{run}/qiime_classify_{sample}.txt"
    shell:
        """
        qiime feature-classifier classify-consensus-vsearch \
            --i-query {input.query} \
            --i-reference-reads {input.ref_seqs} \
            --i-reference-taxonomy {input.taxo} \
            --o-classification {output.qza} \
            --p-threads {threads} \
            --p-perc-identity {params.pident} \
            --p-maxaccepts {params.nseqs} \
            --p-min-consensus {params.lcacons} \
            2>&1 | tee -a {log}
        qiime tools export \
            --input-path {output.qza} \
            --output-path {output.path} 2>&1 | tee -a {log}
        mv {output.path}/taxonomy.tsv {output.out} 2>> {log}
        """


rule qiime_taxlist:
    input:
        "classifications/{run}/qiime/{sample}.qiime.out"
    output:
        "classifications/{run}/qiime/{sample}.qiime.taxlist"
    threads: 1
    log:
        "logs/{run}/qiime_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/qiime_tomat_{sample}.txt"
    conda:
        config["qiime"]["environment"]
    shell:
        """
        cut -f 1,2 {input} | grep -P -v '\\tUnassigned$' | \
        awk -F '\\t|;' -v OFS='\\t' '{{NF=7}}1' | \
        awk 'BEGIN {{ FS = OFS = "\\t" }} {{ for(i=1; i<=NF; i++) if($i == "") $i = "NA" }}; 1' | \
        sed "s/D_[[:digit:]]__//g" | \
        sed "s/^Feature ID.*/\#readid\\tDomain\\tPhylum\\tClass\\tOrder\\tFamily\\tGenus/" \
        > {output} 2> {log}
        """

rule qiime_tomat:
    input:
        list = "classifications/{run}/qiime/{sample}.qiime.taxlist",
    output:
        taxmat = "classifications/{run}/qiime/{sample}.qiime.taxmat",
        otumat = "classifications/{run}/qiime/{sample}.qiime.otumat"
    threads: 1
    log:
        "logs/{run}/qiime_tomat_{sample}.log"
    benchmark:
        "benchmarks/{run}/qiime_tomat_{sample}.txt"
    shell:
        "scripts/tomat.py -l {input.list} 2> {log}"

