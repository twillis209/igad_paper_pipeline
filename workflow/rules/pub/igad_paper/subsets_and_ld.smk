rule subset_sumstats_for_tnfaip3_locus_plot:
    input:
        "results/igad_meta/meta.tsv.gz"
    output:
        sumstats = "results/pub/igad_paper/misc/tnfaip3_locus_snps.tsv",
        ids = "results/pub/igad_paper/misc/tnfaip3_locus_snp_ids.txt"
    params:
        chrom = 6,
        start = 137600000,
        stop = 138000000
    localrule: True
    shell:"""
        zcat {input} | awk 'NR == 1 || $1 == {params.chrom} && $2 >= {params.start} && $2 <= {params.stop}' > {output.sumstats}
        zcat {input} | awk '$1 == {params.chrom} && $2 >= {params.start} && $2 <= {params.stop} {{print $17}}' > {output.ids}
    """

use rule subset_1kGP_data_for_ld_matrix as subset_1kGP_data_for_tnfaip3_locus with:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/pub/igad_paper/misc/tnfaip3_locus_snp_ids.txt"
    output:
        multiext("results/pub/igad_paper/misc/tnfaip3_locus_snps", ".pgen", ".pvar.zst", ".psam")
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/pub/igad_paper/misc/tnfaip3_locus_snps"

use rule calculate_ld_for_subset_about_variant as calculate_ld_for_tnfaip3_locus with:
    input:
        multiext("results/pub/igad_paper/misc/tnfaip3_locus_snps", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/pub/igad_paper/misc/tnfaip3_locus_snps.phased", ".vcor2", ".vcor2.vars")
    params:
        stem = "results/pub/igad_paper/misc/tnfaip3_locus_snps"

use rule merge_sumstats_with_r2 as merge_tnfaip3_sumstats_with_r2 with:
    input:
        gwas = "results/pub/igad_paper/misc/tnfaip3_locus_snps.tsv",
        ld = "results/pub/igad_paper/misc/tnfaip3_locus_snps.phased.vcor2",
        ld_vars = "results/pub/igad_paper/misc/tnfaip3_locus_snps.phased.vcor2.vars"
    output:
        "results/pub/igad_paper/misc/tnfaip3_locus_snps_with_r2_for_{variant_id}.tsv.gz"

use rule subset_sumstats_for_tnfaip3_locus_plot as subset_sumstats_for_cd86_locus_plot with:
    input:
        "results/igad_meta/meta.tsv.gz"
    output:
        sumstats = "results/pub/igad_paper/misc/cd86_locus_snps.tsv",
        ids = "results/pub/igad_paper/misc/cd86_locus_snp_ids.txt"
    params:
        chrom = 3,
        start = 122e6,
        stop = 122.2e6

use rule subset_1kGP_data_for_ld_matrix as subset_1kGP_data_for_cd86_locus with:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/pub/igad_paper/misc/cd86_locus_snp_ids.txt"
    output:
        multiext("results/pub/igad_paper/misc/cd86_locus_snps", ".pgen", ".pvar.zst", ".psam")
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/pub/igad_paper/misc/cd86_locus_snps"

use rule calculate_ld_for_subset_about_variant as calculate_ld_for_cd86_locus with:
    input:
        multiext("results/pub/igad_paper/misc/cd86_locus_snps", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/pub/igad_paper/misc/cd86_locus_snps.phased", ".vcor2", ".vcor2.vars")
    params:
        stem = "results/pub/igad_paper/misc/cd86_locus_snps"

use rule merge_sumstats_with_r2 as merge_cd86_sumstats_with_r2 with:
    input:
        gwas = "results/pub/igad_paper/misc/cd86_locus_snps.tsv",
        ld = "results/pub/igad_paper/misc/cd86_locus_snps.phased.vcor2",
        ld_vars = "results/pub/igad_paper/misc/cd86_locus_snps.phased.vcor2.vars"
    output:
        "results/pub/igad_paper/misc/cd86_locus_snps_with_r2_for_{variant_id}.tsv.gz"

rule subset_sumstats_for_ahi1_tnfaip3_locus_plot:
    input:
        "results/pub/igad_paper/misc/merged_meta_cfdr_and_aux_data.tsv.gz"
    output:
        sumstats = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps.tsv",
        ids = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snp_ids.txt"
    params:
        start = 133e6,
        stop = 140e6,
        chrom = 6
    localrule: True
    shell:"""
        zcat {input} | awk 'NR == 1 || $2 == {params.chrom} && $3 >= {params.start} && $3 <= {params.stop}' > {output.sumstats}
        zcat {input} | awk '$2 == {params.chrom} && $3 >= {params.start} && $3 <= {params.stop} {{print $1}}' > {output.ids}
    """

use rule subset_1kGP_data_for_ld_matrix as subset_1kGP_data_for_ahi1_tnfaip3_locus with:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snp_ids.txt"
    output:
        multiext("results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps", ".pgen", ".pvar.zst", ".psam")
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps"

use rule calculate_ld_for_subset_about_variant as calculate_ld_for_ahi1_tnfaip3_locus with:
    input:
        multiext("results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps.phased", ".vcor2", ".vcor2.vars")
    params:
        stem = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps"

use rule merge_sumstats_with_r2 as merge_ahi1_tnfaip3_sumstats_with_r2 with:
    input:
        gwas = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps.tsv",
        ld = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps.phased.vcor2",
        ld_vars = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps.phased.vcor2.vars"
    output:
        "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps_with_r2_for_{variant_id}.tsv.gz"

use rule subset_sumstats_for_ahi1_tnfaip3_locus_plot as subset_sumstats_for_cd68_locus_plot with:
    output:
        sumstats = "results/pub/igad_paper/misc/cd68_locus_snps.tsv",
        ids = "results/pub/igad_paper/misc/cd68_locus_snp_ids.txt"
    params:
        start = 745e4,
        stop = 77e5,
        chrom = 17
    localrule: True

use rule subset_1kGP_data_for_ld_matrix as subset_1kGP_data_for_cd68_locus with:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/pub/igad_paper/misc/cd68_locus_snp_ids.txt"
    output:
        multiext("results/pub/igad_paper/misc/cd68_locus_snps", ".pgen", ".pvar.zst", ".psam")
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/pub/igad_paper/misc/cd68_locus_snps"

use rule calculate_ld_for_subset_about_variant as calculate_ld_for_cd68_locus with:
    input:
        multiext("results/pub/igad_paper/misc/cd68_locus_snps", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/pub/igad_paper/misc/cd68_locus_snps.phased", ".vcor2", ".vcor2.vars")
    params:
        stem = "results/pub/igad_paper/misc/cd68_locus_snps"

use rule merge_sumstats_with_r2 as merge_cd68_sumstats_with_r2 with:
    input:
        gwas = "results/pub/igad_paper/misc/cd68_locus_snps.tsv",
        ld = "results/pub/igad_paper/misc/cd68_locus_snps.phased.vcor2",
        ld_vars = "results/pub/igad_paper/misc/cd68_locus_snps.phased.vcor2.vars"
    output:
        "results/pub/igad_paper/misc/cd68_locus_snps_with_r2_for_{variant_id}.tsv.gz"

# Reusing CD68 window SNPs
use rule merge_sumstats_with_r2 as merge_tnfsf13_sumstats_with_r2 with:
    input:
        gwas = "results/pub/igad_paper/misc/cd68_locus_snps.tsv",
        ld = "results/pub/igad_paper/misc/cd68_locus_snps.phased.vcor2",
        ld_vars = "results/pub/igad_paper/misc/cd68_locus_snps.phased.vcor2.vars"
    output:
        "results/pub/igad_paper/misc/tnfsf13_locus_snps_with_r2_for_{variant_id}.tsv.gz"

use rule subset_sumstats_for_ahi1_tnfaip3_locus_plot as subset_sumstats_for_inava_locus_plot with:
    output:
        sumstats = "results/pub/igad_paper/misc/inava_locus_snps.tsv",
        ids = "results/pub/igad_paper/misc/inava_locus_snp_ids.txt"
    params:
        start = 200.8e6,
        stop = 201.1e6,
        chrom = 1
    localrule: True

use rule subset_1kGP_data_for_ld_matrix as subset_1kGP_data_for_inava_locus with:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/pub/igad_paper/misc/inava_locus_snp_ids.txt"
    output:
        multiext("results/pub/igad_paper/misc/inava_locus_snps", ".pgen", ".pvar.zst", ".psam")
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/pub/igad_paper/misc/inava_locus_snps"

use rule calculate_ld_for_subset_about_variant as calculate_ld_for_inava_locus with:
    input:
        multiext("results/pub/igad_paper/misc/inava_locus_snps", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/pub/igad_paper/misc/inava_locus_snps.phased", ".vcor2", ".vcor2.vars")
    params:
        stem = "results/pub/igad_paper/misc/inava_locus_snps"

use rule merge_sumstats_with_r2 as merge_inava_sumstats_with_r2 with:
    input:
        gwas = "results/pub/igad_paper/misc/inava_locus_snps.tsv",
        ld = "results/pub/igad_paper/misc/inava_locus_snps.phased.vcor2",
        ld_vars = "results/pub/igad_paper/misc/inava_locus_snps.phased.vcor2.vars"
    output:
        "results/pub/igad_paper/misc/inava_locus_snps_with_r2_for_{variant_id}.tsv.gz"

rule look_up_gws_snps_in_tnfaip3_signal:
    input:
        "results/igad_meta/meta.tsv.gz"
    output:
        "results/pub/igad_paper/misc/tnfaip3_signal_snps.tsv"
    params:
        # NB: From Ensembl
        wakmar2_interval = [137823673, 137868233],
        tnfaip3_interval = [137867214, 137883314]
    threads: 8
    resources:
    localrule: True
    conda: env_path("locuszoomr.yaml")
    script: script_path("igad_meta/look_up_gws_snps_in_tnfaip3_signal.R")

rule fetch_gws_associations_for_tnfaip3_gws_snps:
    input:
        "results/pub/igad_paper/misc/tnfaip3_signal_snps.tsv"
    output:
        "results/pub/igad_paper/misc/tnfaip3_signal_snps_phewas_lookup.tsv"
    localrule: True
    script: script_path("pub/igad_paper/misc/lookup_tnfaip3_phewas.py")
