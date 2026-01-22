# ========== ref.smk ==========
# Create working copy + symlink reference + index

READ_SET = config["read_set"]

REF_SRC = os.path.join(config["project_root"], "data", "references", config["reference"])
REF_FASTA = f"results/{READ_SET}/reference/{config['reference']}"

rule symlink_reference:
    output:
        REF_FASTA
    shell:
        "mkdir -p results/{config[read_set]}/reference && "
        "ln -sf {REF_SRC} {output}"

rule faidx:
    input:
        REF_FASTA
    output:
        REF_FASTA + ".fai"
    shell:
        "samtools faidx {input}"

rule dict:
    input:
        REF_FASTA
    output:
        REF_FASTA.replace('.fa', '.dict')
    shell:
        """
        gatk CreateSequenceDictionary \
            -R {input} \
            -O {output}
        """

rule bwa_index:
    input:
        REF_FASTA
    output:
        REF_FASTA + ".pac"   # bwa‑mem2 index artifact
    shell:
        "bwa-mem2 index {input}"
