#checkpoint touch_lead_snp_files:
#    input:
#        "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps.distance_clumped"
#    output:
#        directory("results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps")
#    params:
#        out_dir = "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps"
#    localrule: True
#    run:
#        daf = pd.read_csv(input[0], sep = '\t', header = 0)
#
#        shell(f"mkdir {params.out_dir}")
#
#        for x in daf['SNPID']:
#            shell(f"touch {params.out_dir}/{x.replace(':', '_')}.done")

rule subset_snps_for_lead_snp_window:
    input:
        multiext("results/1kG/hg38/eur/snps_only/qc/all/merged", ".pgen", ".pvar.zst", ".psam")
    output:
        temp(multiext("results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood", ".pgen", ".pvar.zst", ".psam")),
        range_file = temp("results/1kG/hg38/eur/snps_only/{snp_id}/range.txt")
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/qc/all/merged",
        out_stem = "results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood",
        window_width = 1e6
    threads: 8
    resources:
    group: "gwas"
    run:
        chrom = wildcards.snp_id.split('_')[0]
        start = int(wildcards.snp_id.split('_')[1])-int(params.window_width/2)
        stop = int(wildcards.snp_id.split('_')[1])+int(params.window_width/2)

        with open(output.range_file, 'w') as outfile:
            outfile.write(f"{chrom}\t{start}\t{stop}\t{wildcards.snp_id}\n")

        shell("plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --extract range {output.range_file} --make-pgen vzs --out {params.out_stem}")

# NB: --r2 is not yet implemented in plink2
rule convert_lead_snp_window_to_bfiles:
    input:
        multiext("results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood", ".pgen", ".pvar.zst", ".psam")
    output:
        temp(multiext("results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood", ".bed", ".bim", ".fam"))
    params:
        stem = "results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood"
    threads: 8
    resources:
    group: "gwas"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.stem} vzs --make-bed --out {params.stem}"

rule calculate_ld_in_lead_snp_window:
    input:
        multiext("results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood", ".bed", ".bim", ".fam")
    output:
        temp("results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood.ld.gz")
    params:
        stem = "results/1kG/hg38/eur/snps_only/{snp_id}/neighbourhood"
    threads: 8
    resources:
    group: "gwas"
    shell:
        "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.stem} --r2 inter-chr gz --out {params.stem}"
