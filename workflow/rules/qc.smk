# ========== qc.smk ==========

READ_SET = config["read_set"]
SAMPLES_DF = pd.read_csv("workflow/config/samples.tsv", sep="\t")

rule fastqc_raw:
    input:
        R1=lambda wildcards: str(SAMPLES_DF.set_index("sample").loc[wildcards.sample, "R1"]),
        R2=lambda wildcards: str(SAMPLES_DF.set_index("sample").loc[wildcards.sample, "R2"]),
    output:
        f"results/{READ_SET}/qc/fastqc_raw/{{sample}}_R1_fastqc.zip",
        f"results/{READ_SET}/qc/fastqc_raw/{{sample}}_R2_fastqc.zip",
    threads: 4
    shell:
        """
        mkdir -p results/{READ_SET}/qc/fastqc_raw
        fastqc -t {threads} -o results/{READ_SET}/qc/fastqc_raw {input.R1} {input.R2}
        """
