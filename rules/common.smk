rule common_download_db:
    output:
        ref_tax = os.path.join(DBPATH,"common/ref-taxonomy.txt"),
        ref_seqs = os.path.join(DBPATH,"common/ref-seqs.fna"),
        aln = os.path.join(DBPATH,"common/ref-seqs.aln"),
        names = os.path.join(DBPATH,"common/taxonomy/names.dmp"),
        nodes = os.path.join(DBPATH,"common/seqid2taxid.map"),
        seq_map = os.path.join(DBPATH,"common/taxonomy/nodes.dmp"),
        silva_seqs = os.path.join(DBPATH,"common/data/SILVA_138.1_SSURef_NR99_tax_silva.fasta")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["common"]["dbmemory"]
    params:
        ssu = config["common"]["ssu"],
        ssu_uppercase = config["common"]["ssu"].upper()
    log:
        "logs/common_download_db.log"
    benchmark:
        "benchmarks/common_download_db.txt"
    conda:
        os.path.join(ENVDIR,config["common"]["environment"])
    shell:
        """
        #prep folders 
        mkdir -p {DBPATH}/common/data {DBPATH}/common/taxonomy {DBPATH}/common/library

        #download data
        wget -P {DBPATH}/common https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_{params.ssu_uppercase}Ref_NR99_tax_silva.fasta.gz
        wget -P {DBPATH}/common https://ftp.arb-silva.de/release_138.1/Exports/taxonomy/tax_slv_{params.ssu}_138.1.txt.gz
        wget -P {DBPATH}/common https://ftp.arb-silva.de/release_138.1/Exports/taxonomy/tax_slv_{params.ssu}_138.1.acc_taxid.gz
        wget -P {DBPATH}/common https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_{params.ssu_uppercase}Ref_NR99_tax_silva_full_align_trunc.fasta.gz

        gzip -d {DBPATH}/common/*gz

        #get a taxonomy file
        perl {SRCDIR}/build_silva_taxonomy.pl {DBPATH}/common/tax_slv_{params.ssu}_138.1.txt

        mv names.dmp nodes.dmp {DBPATH}/common/taxonomy
        mv {DBPATH}/common/SILVA_138.1_{params.ssu_uppercase}Ref_NR99_tax_silva.fasta {DBPATH}/common/tax_slv_{params.ssu}_138.1.txt {DBPATH}/common/data 
        mv {DBPATH}/common/tax_slv_{params.ssu}_138.1.acc_taxid {DBPATH}/common/seqid2taxid.map

        #cleanup aln
        cut -f1 -d " " {DBPATH}/common/SILVA_138.1_{params.ssu_uppercase}Ref_NR99_tax_silva_full_align_trunc.fasta > {DBPATH}/common/ref-seqs.aln
        rm {DBPATH}/common/SILVA_138.1_{params.ssu_uppercase}Ref_NR99_tax_silva_full_align_trunc.fasta

        #clean U 
        sed -e '/^>/!y/U/T/' {DBPATH}/common/data/SILVA_138.1_{params.ssu_uppercase}Ref_NR99_tax_silva.fasta | cut -f1 -d " " > {DBPATH}/common/ref-seqs.fna

        #create symlink for kraken
        ln -s {DBPATH}/common/ref-seqs.fna {DBPATH}/common/library/ref-seqs.fna

        taxonkit lineage --data-dir {DBPATH}/common/taxonomy <(awk '{{print $1}}' {DBPATH}/common/taxonomy/names.dmp | sort | uniq) | \
            taxonkit reformat -r NA -a -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" --data-dir {DBPATH}/common/taxonomy | \
            awk 'BEGIN {{ FS = OFS = "\t" }} {{print $1,"D_0__"$3";D_1__"$4";D_2__"$5";D_3__"$6";D_4__"$7";D_5__"$8";D_6__"$9}}' > {DBPATH}/common/taxid_to_tax.txt
        
        LC_ALL=C join -1 2 -2 1 -t $'\t' <(LC_ALL=C sort -k 2 {DBPATH}/common/seqid2taxid.map) <(LC_ALL=C sort -k1 {DBPATH}/common/taxid_to_tax.txt) \
            | awk -F "\t" -v OFS="\t" '{{print $2, $3}}' > {DBPATH}/common/ref-taxonomy.txt

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
