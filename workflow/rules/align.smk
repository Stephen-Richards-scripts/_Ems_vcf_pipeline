# ========== align.smk ==========
READ_SET = config["read_set"]
REF_FASTA = f"results/{READ_SET}/reference/{config['reference']}"

rule align_sort:
    input:
        R1=f"results/{READ_SET}/trimmed/{{sample}}.R1.trimmed.fastq.gz",
        R2=f"results/{READ_SET}/trimmed/{{sample}}.R2.trimmed.fastq.gz",
        REF_FASTA=REF_FASTA,
        FAI=REF_FASTA + ".fai"
    output:
        bam=f"results/{READ_SET}/alignments/{{sample}}.sorted.bam"
    threads: config["threads_per_sample"]
    shell:
        """
        mkdir -p results/{READ_SET}/alignments

        bwa-mem2 mem -M -t {threads} \
          -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:{READ_SET}\\tPU:{wildcards.sample}" \
          {input.REF_FASTA} {input.R1} {input.R2} \
        | samtools view -b - \
        | samtools sort -@ {threads} -m 1500M -T $TMPDIR/{wildcards.sample} \
          -o {output.bam}
        """


rule markdup_to_cram:
    input:
        bam=f"results/{READ_SET}/alignments/{{sample}}.sorted.bam",
        REF_FASTA=REF_FASTA
    output:
        cram=f"results/{READ_SET}/alignments/{{sample}}.dedup.cram",
        crai=f"results/{READ_SET}/alignments/{{sample}}.dedup.cram.crai",
        metrics=f"results/{READ_SET}/alignments/{{sample}}.markdup.metrics.txt"
    threads: 4
    shell:
        """
        # 1) Mark duplicates -> BAM (GATK cannot write CRAM)
        gatk MarkDuplicates \
            -I {input.bam} \
            -O results/{READ_SET}/alignments/{wildcards.sample}.dedup.bam \
            -M {output.metrics} \
            --CREATE_INDEX true \
            --TMP_DIR $TMPDIR

        # 2) Convert BAM -> CRAM
        samtools view -T {input.REF_FASTA} -C \
            -o {output.cram} \
            results/{READ_SET}/alignments/{wildcards.sample}.dedup.bam

        samtools index {output.cram}

        # 3) Cleanup BAM intermediates
        rm results/{READ_SET}/alignments/{wildcards.sample}.dedup.bam
        rm results/{READ_SET}/alignments/{wildcards.sample}.dedup.bam.bai
        """

