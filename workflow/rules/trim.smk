# ========== trim.smk ==========
READ_SET = config["read_set"]
ADAPTERS = os.path.join(config["project_root"], "data", "adapters", "adapters.fasta")

rule trim_reads:
    input:
        R1=lambda wc: str(SAMPLES_DF.set_index("sample").loc[wc.sample, "R1"]),
        R2=lambda wc: str(SAMPLES_DF.set_index("sample").loc[wc.sample, "R2"])
    output:
        R1t = f"results/{READ_SET}/trimmed/{{sample}}.R1.trimmed.fastq.gz",
        R2t = f"results/{READ_SET}/trimmed/{{sample}}.R2.trimmed.fastq.gz"
    threads: config["threads_per_sample"]
    shell:
        """
        mkdir -p results/{READ_SET}/trimmed

        cutadapt \
          -j {threads} \
          -q 15,15 \
          -m 36 \
          -g file:{ADAPTERS} \
          -G file:{ADAPTERS} \
          -a file:{ADAPTERS} \
          -A file:{ADAPTERS} \
          -o {output.R1t} \
          -p {output.R2t} \
          {input.R1} {input.R2}
        """
