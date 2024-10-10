import requests
import json
import pandas as pd

def get_per_iteration_cfdr_files(w):
    checkpoint_output = checkpoints.write_out_per_iteration_cfdr.get(**w).output[0]

    return expand("results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/per_iteration/it_{it}.tsv.gz", prin_trait = w.prin_trait, aux_traits = w.aux_traits, variant_set = w.variant_set, variant_type = w.variant_type, window_size = w.window_size, r2 = w.r2, cfdr_threshold = w.cfdr_threshold, it = glob_wildcards(os.path.join(checkpoint_output, "it_{it}.tsv.gz")).it)

def get_per_iteration_cfdr_manhattans(w):
    checkpoint_output = checkpoints.write_out_per_iteration_cfdr.get(**w).output[0]

    return expand("results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration/it_{it}_annotated_manhattan.distance_clumped.png", prin_trait = w.prin_trait, aux_traits = w.aux_traits, variant_set = w.variant_set, variant_type = w.variant_type, window_size = w.window_size, r2 = w.r2, cfdr_threshold = w.cfdr_threshold, clump_threshold = w.clump_threshold, it = glob_wildcards(os.path.join(checkpoint_output, "it_{it}.tsv.gz")).it)

def get_cfdr_manhattan_title(w):
    pretty_aux_names = ', '.join([get_metadata_field(x, 'pretty_name') for x in w.aux_traits.split('_and_')])

    return '%s conditioned on %s, iteration %d' % (get_metadata_field(w.prin_trait, 'pretty_name'), pretty_aux_names, w.iteration)

rule run_iterative_cfdr:
    input:
        gwas_file = "results/merged_gwas/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/merged_with_prune.tsv.gz"
    output:
        "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/cfdr.tsv.gz"
    log:
        v = "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/v_and_sub_dat.RData"
    threads: 16
    resources:
        runtime = 300,
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        prin_col = 'P',
        aux_cols = lambda wildcards: [f'P.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('_and_'))],
        prune_cols = lambda wildcards: [f'prune_in.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('_and_'))],
        v_cols = lambda wildcards: [f'v.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('_and_'))],
        p_threshold = lambda wildcards: pow(10, -int(wildcards.cfdr_threshold)),
    group: "cfdr"
    #conda: env_path("cfdr.yaml")
    script:
        script_path("cfdr/iterative_cfdr.R")

checkpoint write_out_per_iteration_cfdr:
    input:
        "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/cfdr.tsv.gz"
    output:
        directory("results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/per_iteration")
    params:
        aux_traits = lambda w: w.aux_traits.split('_and_'),
        out_format = "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/per_iteration/it_%d.tsv.gz",
        bp_col = 'BP38',
        chr_col = 'CHR38',
        ref_col = 'REF',
        alt_col = 'ALT',
        prin_col = 'P',
        snp_col = 'SNPID',
        aux_cols = lambda wildcards: [f'P.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('_and_'))],
        prune_cols = lambda wildcards: [f'prune_in.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('_and_'))],
        v_cols = lambda wildcards: [f'v.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('_and_'))]
    group: "cfdr"
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script:
        script_path("cfdr/write_out_per_iteration_cfdr.R")

rule write_out_per_iteration_cfdr_stats:
    input:
        get_per_iteration_cfdr_files
    output:
        "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/per_iteration.done"
    localrule: True
    shell: "touch {output}"

use rule distance_clump_gwas as distance_clump_cfdr with:
    input:
        gwas = "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/cfdr.tsv.gz"
    output:
        "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.distance_clumped"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        p_col = lambda w: 'v.%d' % len(w.aux_traits.split('_and_')),
        index_threshold = lambda wildcards: 5e-8 if wildcards.clump_threshold == 'gws' else 1e-5,
        distance_window = 2e6

use rule distance_clump_cfdr as distance_clump_per_iteration_cfdr with:
    input:
        gwas = "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/per_iteration/it_{it}.tsv.gz"
    output:
        "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration/it_{it}_lead_snps.distance_clumped"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        p_col = lambda w: f'v.{w.it}',
        index_threshold = lambda wildcards: 5e-8 if wildcards.clump_threshold == 'gws' else 1e-5,
        distance_window = 2e6

use rule annotate_lead_snps as annotate_lead_snps_for_cfdr with:
    input:
        "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}"
    output:
        annotations = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}.annotations",
        rsIDs = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}.rsIDs"

use rule annotate_lead_snps_for_cfdr as annotate_lead_snps_for_per_iteration_cfdr with:
    input:
        "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration/it_{it}_lead_snps.{clumping_suffix}"
    output:
        annotations = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration/it_{it}_lead_snps.{clumping_suffix}.annotations",
        rsIDs = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration/it_{it}_lead_snps.{clumping_suffix}.rsIDs"

rule draw_manhattan_with_lead_snp_annotation_for_cfdr:
    input:
        gwas = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/cfdr.tsv.gz",
        rsIDs = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}.rsIDs"
    output:
        "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/annotated_manhattan.{clumping_suffix}.png"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = lambda w: 'v.%d' % len(w.aux_traits.split('_and_')),
        snp_col = 'SNPID',
        title = ''
    threads: 8
    resources:
        runtime = 20
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script:
        script_path("gwas/plot_gwas_manhattan.R")

use rule draw_manhattan_with_lead_snp_annotation_for_cfdr as draw_manhattan_with_lead_snp_annotation_for_per_iteration_cfdr with:
    input:
        gwas = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/per_iteration/it_{it}.tsv.gz",
        rsIDs = "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration/it_{it}_lead_snps.{clumping_suffix}.rsIDs"
    output:
        "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration/it_{it}_annotated_manhattan.{clumping_suffix}.png"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = lambda w: f'v.{w.it}',
        snp_col = 'SNPID',
        title = ''

rule draw_per_iteration_annotated_manhattans:
    input:
        get_per_iteration_cfdr_manhattans
    output:
        "results/cfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/per_iteration_annotated_manhattans.done"
    shell: "touch {output}"

checkpoint touch_cfdr_lead_snp_files:
    input:
        "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.distance_clumped"
    output:
        directory("results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps")
    params:
        out_dir = "results/cfdr/{prin_trait}_on_{aux_traits}/left/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t', header = 0)

        shell(f"mkdir {params.out_dir}")

        for x in daf['SNPID']:
            shell(f"touch {params.out_dir}/{x.replace(':', '_')}.done")

use rule subset_snps_for_lead_snp_window as subset_snps_for_cfdr_lead_snp_window with:
    output:
        temp(multiext("results/processed_gwas/{trait}/{snp_set}/{clump_threshold}/lead_snps/{snp_id}/neighbourhood", ".pgen", ".pvar.zst", ".psam")),
        range_file = temp("results/processed_gwas/{trait}/{snp_set}/{clump_threshold}/lead_snps/{snp_id}/range.txt")
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/qc/all/merged",
        out_stem = "results/processed_gwas/{trait}/{snp_set}/{clump_threshold}/lead_snps/{snp_id}/neighbourhood",
        window_width = 1e6

use rule distance_clump_gwas as distance_clump_fcfdr with:
    input:
        gwas = "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/binary_cfdr_output.tsv.gz"
    output:
        "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.distance_clumped"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        p_col = 'v.1',
        index_threshold = lambda wildcards: 5e-8 if wildcards.clump_threshold == 'gws' else 1e-5,
        distance_window = 2e6

use rule annotate_lead_snps as annotate_lead_snps_for_fcfdr with:
    input:
        "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}"
    output:
        annotations = "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}.annotations",
        rsIDs = "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}.rsIDs"

use rule draw_manhattan_with_lead_snp_annotation_for_cfdr as draw_manhattan_with_lead_snp_annotation_for_fcfdr with:
    input:
        gwas = "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/binary_cfdr_output.tsv.gz",
        rsIDs = "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/lead_snps.{clumping_suffix}.rsIDs"
    output:
        "results/fcfdr/{prin_trait}_on_{aux_traits}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{cfdr_threshold}/{clump_threshold}/annotated_manhattan.{clumping_suffix}.png"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = 'v.1',
        snp_col = 'SNPID',
        title = ''
