
# ========== call_gatk.smk ==========

READ_SET = config["read_set"]
REF_FASTA = f"results/{READ_SET}/reference/{config['reference']}"

CRAMS = expand(f"results/{READ_SET}/alignments/{{sample}}.dedup.cram",
               sample=SAMPLES)

rule hapcaller_gvcf:
    input:
        cram=f"results/{READ_SET}/alignments/{{sample}}.dedup.cram",
        crai=f"results/{READ_SET}/alignments/{{sample}}.dedup.cram.crai",
        REF_FASTA
    output:
        gvcf=f"results/{READ_SET}/gvcf/{{sample}}.g.vcf.gz"
    threads: 4
    shell:
        """
        mkdir -p results/{READ_SET}/gvcf

        gatk HaplotypeCaller \
            -R {input.REF_FASTA} \
            -I {input.cram} \
            -O {output.gvcf} \
            -ERC GVCF \
            --native-pair-hmm-threads 4

        bcftools index -t {output.gvcf}
        """

rule combine_gvcfs:
    input:
        gvcfs=expand(f"results/{READ_SET}/gvcf/{{sample}}.g.vcf.gz", sample=SAMPLES)
    output:
        combined=f"results/{READ_SET}/variants/gatk/combined.g.vcf.gz"
    shell:
        """
        mkdir -p results/{READ_SET}/variants/gatk

        gatk CombineGVCFs \
            -R {REF_FASTA} \
            {' '.join(['--variant ' + g for g in input.gvcfs])} \
            -O {output.combined}

        bcftools index -t {output.combined}
        """

rule genotype_gvcfs:
    input:
        combined=f"results/{READ_SET}/variants/gatk/combined.g.vcf.gz"
    output:
        vcf=f"results/{READ_SET}/variants/gatk/cohort.vcf.gz",
        tbi=f"results/{READ_SET}/variants/gatk/cohort.vcf.gz.tbi"
    shell:
        """
        gatk GenotypeGVCFs \
            -R {REF_FASTA} \
            -V {input.combined} \
            -O {output.vcf}

        bcftools index -t {output.vcf}
        """

