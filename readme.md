## General

NanoClass2 is a taxonomic meta-classifier for 16S/18S amplicon sequencing data generated with the Oxford Nanopore MinION. The first iteration of this workflow, NanoClass, was originally developed by Evelien Jongepier and the original version can be found [here](https://ejongepier.github.io/NanoClass/getting_started.html).

With a single command, you can run up to 11 classification tools on multiple samples in parallel, including BLASTN, Centrifuge, Kraken2, IDTAXA, MegaBLAST, dcMegaBLAST, Minimap2, Mothur, QIIME2, RDP and SPINGO. Read preparation steps, such as quality trimming, length filtering and sub-sampling, are an integral part of the pipeline.

NanoClass automatically installs all software packages and dependencies, downloads and builds required taxonomic databases and runs the analysis on your samples.

For detailed instructions on how to use NanoClass2, please visit [the manual](https://ndombrowski.github.io/NanoClass2/).


**Cite**

If you found this software useful and used it in your own research, please cite us using the following doi:

[![DOI](https://zenodo.org/badge/722016433.svg)](https://doi.org/10.5281/zenodo.14264106)
