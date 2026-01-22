# ========== call_bcftools.smk ==========

READ_SET = config["read_set"]
REF_FASTA = f"results/{READ_SET}/reference/{config['reference']}"

CRAMS = expand(f"results/{READ_SET}/alignments/{{sample}}.dedup.cram",
               sample=SAMPLES)


rule bcftools_call:
    input:
        ref = REF_FASTA,
        crams = CRAMS
    output:
        vcf = f"results/{READ_SET}/variants/bcftools/cohort.vcf.gz",
        tbi = f"results/{READ_SET}/variants/bcftools/cohort.vcf.gz.tbi"
    threads: 16
    shell:
        """
        mkdir -p results/{READ_SET}/variants/bcftools

        bcftools mpileup -Ou \
            -f {input.ref} \
            -q 20 -Q 20 \
            -a FORMAT/DP,FORMAT/AD \
            {input.crams} \
        | bcftools call -mv -Oz -o {output.vcf}

        bcftools index -t {output.vcf}
        """
	
CRAMS = expand(f"results/{READ_SET}/alignments/{{sample}}.dedup.cram",
               sample=SAMPLES)
