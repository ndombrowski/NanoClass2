import re
import os
import yaml
from typing import List
from pathlib import Path
import pandas as pd
from snakemake.utils import min_version, validate

min_version("5.23.0")
report: "report/workflow.rst"

# ======================================================
# Config files
# ======================================================

#configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")

smpls = pd.read_csv(config["samples"], dtype=str).set_index(["run", "sample"], drop=False)
smpls.index = smpls.index.set_levels([i.astype(str) for i in smpls.index.levels])
validate(smpls, schema="schemas/samples.schema.yaml")

wildcard_constraints:
    sample = '[A-Za-z0-9]+',
    run = '[A-Za-z0-9]+'


# ======================================================
# directories
# ======================================================

ROOTDIR = srcdir(".")

DBPATH = os.path.expandvars(config['db_path'])
if not os.path.isabs(DBPATH):
    DBPATH = os.path.join(ROOTDIR, DBPATH)
if not os.path.exists(DBPATH):
    os.makedirs(DBPATH)

SRCDIR = srcdir("scripts")
ENVDIR = srcdir("envs")
OUTPUTDIR = os.path.abspath(os.path.expandvars(config['output_dir']))

workdir:
    OUTPUTDIR



# ======================================================
# Rules
# ======================================================

optrules = [] 
optrules.extend(["plots/precision.pdf"] if len(config["methods"]) > 2 else [])


rule all:
    input:
        expand("data/{samples.run}/porechopped/{samples.sample}.trimmed.fastq.gz",
            samples=smpls.itertuples()
        ),
        expand("plots/{samples.run}/nanofilt/{samples.sample}.filtered.pdf",
            samples=smpls.itertuples(), method=config["methods"]
        ),
        expand("stats/{samples.run}/nanofilt/{samples.sample}.filtered.txt",
            samples=smpls.itertuples(), method=config["methods"]
        ),
        expand("classifications/{samples.run}/{method}/{samples.sample}.{method}.taxmat",
            samples=smpls.itertuples(), method=config["methods"]
        ),
        expand("plots/{absrel}-Genus-by-{grouper}.pdf", 
            absrel = ["aabund","rabund"], grouper=config["common"]["group-by"]
        ), 
        expand("plots/runtime-by-{grouper}.pdf", grouper=config["common"]["group-by"]
        ),
        optrules


rule filter:
    input:
        expand("plots/{samples.run}/nanofilt/{samples.sample}.filtered.pdf",
            samples=smpls.itertuples(), method=config["methods"]
        ),
        expand("stats/{samples.run}/nanofilt/{samples.sample}.filtered.txt",
            samples=smpls.itertuples(), method=config["methods"]
        )



# ======================================================
# Functions and Classes
# ======================================================


onsuccess:
    print("NanoClass finished!")
    print("To generate a report run: snakemake --report report/NanoClass.zip")


onerror:
    print("Note the path to the log file for debugging.")
    print("Documentation is available at: https://ejongepier.github.io/NanoClass/")
    print("Issues can be raised at: https://github.com/ejongepier/NanoClass/issues")


# ======================================================
# Global variables
# ======================================================

#======================================================
# Include
#======================================================

include: "rules/preprocess.smk"
include: "rules/common.smk"
include: "rules/blastn.smk"
include: "rules/centrifuge.smk"
include: "rules/dcmegablast.smk"
include: "rules/idtaxa.smk"
include: "rules/kraken.smk"
include: "rules/mapseq.smk"
include: "rules/megablast.smk"
include: "rules/minimap.smk"
include: "rules/mothur.smk"
include: "rules/qiime.smk"
include: "rules/rdp.smk"
include: "rules/spingo.smk"
