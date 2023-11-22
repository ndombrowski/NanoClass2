rule common_download_db:
    output:
        ref_tax = os.path.join(DBPATH,"common/ref-taxonomy.txt"),
        ref_seqs = os.path.join(DBPATH,"common/ref-seqs.fna"),
        ref_aln = os.path.join(DBPATH,"common/ref-seqs.aln")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["common"]["dbmemory"]
    params:
        url = config["common"]["dburl"],
        ssu = config["common"]["ssu"]
    log:
        "logs/common_download_db.log"
    benchmark:
        "benchmarks/common_download_db.txt"
    conda:
        os.path.join(ENVDIR,config["common"]["environment"])
    shell:
        """
        wget -q -O {DBPATH}/common/db.zip {params.url}
        
        unzip -q -p -j {DBPATH}/common/db.zip \
            */taxonomy/{params.ssu}_only/99/majority_taxonomy_7_levels.txt \
            > {output.ref_tax} 2> {log}
        
        unzip -q -p -j {DBPATH}/common/db.zip \
            */rep_set/rep_set_{params.ssu}_only/99/silva_*_99_{params.ssu}.fna \
            > {output.ref_seqs} 2>> {log}
        
        unzip -q -p -j {DBPATH}/common/db.zip \
            */rep_set_aligned/99/99_alignment.fna.zip > tmp.aln.zip
        
        unzip -q -p tmp.aln.zip > {output.ref_aln} 2>> {log}
        
        rm {DBPATH}/common/db.zip tmp.aln.zip 2>> {log}
        """

rule common_plot_tax:
    input:
        expand("classifications/{smpls.run}/{method}/{smpls.sample}.{method}.taxmat",
            method = config["methods"], smpls =  smpls.itertuples()
        )
    output:
        report(expand("plots/{absrel}-Phylum-by-{grouper}.pdf", 
            absrel = ["aabund","rabund"], grouper = config["common"]["group-by"]),
            caption="../report/fig-phylum.rst", category="Classification"
        ),
        report(expand("plots/{absrel}-Class-by-{grouper}.pdf", 
            absrel = ["aabund","rabund"], grouper = config["common"]["group-by"]),
            caption="../report/fig-class.rst", category="Classification"
        ),
        report(expand("plots/{absrel}-Order-by-{grouper}.pdf",
            absrel = ["aabund","rabund"], grouper = config["common"]["group-by"]),
            caption="../report/fig-order.rst", category="Classification"
        ),
        report(expand("plots/{absrel}-Family-by-{grouper}.pdf",
            absrel = ["aabund","rabund"], grouper = config["common"]["group-by"]),
            caption="../report/fig-family.rst", category="Classification"
        ),
        report(expand("plots/{absrel}-Genus-by-{grouper}.pdf",
            absrel = ["aabund","rabund"], grouper = config["common"]["group-by"]),
            caption="../report/fig-genus.rst", category="Classification"
        )
    threads: 1
    log:
        "logs/common_plot_tax.log"
    benchmark:
        "benchmarks/common_plot_tax.txt"
    conda:
        os.path.join(ENVDIR,config["common"]["environment"])
    shell:
        """
        mkdir -p ./plots
        Rscript {SRCDIR}/barplot.R {input} 2> {log}
        """

rule common_get_precision:
    input:
        expand("classifications/{smpls.run}/{method}/{smpls.sample}.{method}.taxlist",
            method = config["methods"], smpls =  smpls.itertuples()
        )
    output:
        expand("classifications/{smpls.run}/{method}/{smpls.sample}.{method}.precision",
            method = config["methods"], smpls =  smpls.itertuples()
        )
    threads: 1
    log:
        "logs/common_get_precision.log"
    benchmark:
        "benchmarks/common_get_precision.txt"
    conda:
        os.path.join(ENVDIR,config["common"]["environment"])
    shell:
        "{SRCDIR}/toconsensus.py -l {input} 2> {log}"


rule common_plot_precision:
    input:
        expand("classifications/{smpls.run}/{method}/{smpls.sample}.{method}.precision",
            method = config["methods"], smpls =  smpls.itertuples()
        )
    output:
        report("plots/precision.pdf", caption="../report/fig-precision.rst", category="Precision")
    threads: 1
    log:
        "logs/common_plot_precision.log"
    benchmark:
        "benchmarks/common_plot_precision.txt"
    conda:
        os.path.join(ENVDIR,config["common"]["environment"])
    shell:
        """
        mkdir -p ./plots
        Rscript {SRCDIR}/lineplot.R {input} 2> {log}
        """


rule common_plot_runtime:
    input:
        expand("benchmarks/{smpls.run}/{method}_classify_{smpls.sample}.txt",
            method = config["methods"], smpls =  smpls.itertuples()
        )
    output:
        report(expand("plots/runtime-by-{grouper}.pdf", grouper = config["common"]["group-by"]),
            caption="../report/fig-runtime.rst", category="Runtime"
        ),
        report(expand("plots/runtime_log-by-{grouper}.pdf", grouper = config["common"]["group-by"]),
            caption="../report/fig-runtime-log.rst", category="Runtime"
        )
    threads: 1
    log:
        "logs/common_plot_runtime.log"
    benchmark:
        "benchmarks/common_plot_runtime.txt"
    conda:
        os.path.join(ENVDIR,config["common"]["environment"])
    shell:
        """
        mkdir -p ./plots
        mkdir -p ./tables
        Rscript {SRCDIR}/timeplot.R {input} 2> {log}
        """
