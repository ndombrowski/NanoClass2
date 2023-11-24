## General

NanoClass2 is a taxonomic meta-classifier for 16S/18S amplicon sequencing data generated with the Oxford Nanopore MinION. The first iteration of this workflow, NanoClass, was originally developed by Evelien Jongepier and the original version can be found [here](https://ejongepier.github.io/NanoClass/getting_started.html).

With a single command, you can run upto 11 classification tools on multiple samples in parallel, including BLASTN, Centrifuge, Kraken2, IDTAXA, MegaBLAST, dcMegaBLAST, Minimap2, Mothur, QIIME2, RDP and SPINGO. Read preparation steps, such as quality trimming, length filtering and sub-sampling, are an integral part of the pipeline.

NanoClass automatically installs all software packages and dependencies, downloads and builds required taxonomic databases and runs the analysis on your samples.

For detailed instructions on how to use NanoClass2, please visit [the manual](https://ndombrowski.github.io/NanoClass2/).

**Please notice:**

NanoClass is using Kraken2 as one of the classifiers but users have noticed that when using older Kraken versions there is a problem when downloading the Silva database.

We plan to make kraken2 available again but since this means updating the used Silva database throughout the workflow this is still in progress. For now, we recommend running this workflow without Kraken2 if you setup NanoClass2 for the first time.
