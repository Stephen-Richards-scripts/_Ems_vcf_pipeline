# ================= Snakefile =====================
# Entry point for Em's VCF pipeline

import os
import pandas as pd

configfile: "workflow/config/config.yaml"

SAMPLES_DF = pd.read_csv("workflow/config/samples.tsv", sep="\t")
SAMPLES = list(SAMPLES_DF["sample"])

READ_SET     = config["read_set"]
REFERENCE_BN = config["reference"]                      # basename
REF_FASTA    = f"results/{READ_SET}/reference/{REFERENCE_BN}"

CALLER       = config.get("caller", "bcftools")

include: "workflow/rules/ref.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/trim.smk"
include: "workflow/rules/align.smk"

if CALLER == "bcftools":
    include: "workflow/rules/call_bcftools.smk"
elif CALLER == "gatk":
    include: "workflow/rules/call_gatk.smk"
else:
    raise ValueError(f"Unknown caller: {CALLER}")

# FINAL TARGETS
rule all:
    input:
        # reference indices
        REF_FASTA + ".fai",
        REF_FASTA.replace(".fa", ".dict"),

        # alignment outputs
        expand(f"results/{READ_SET}/alignments/{{sample}}.dedup.cram", sample=SAMPLES),
        expand(f"results/{READ_SET}/alignments/{{sample}}.dedup.cram.crai", sample=SAMPLES),

        # variant calling final output
        f"results/{READ_SET}/variants/{CALLER}/cohort.vcf.gz",
        f"results/{READ_SET}/variants/{CALLER}/cohort.vcf.gz.tbi"
