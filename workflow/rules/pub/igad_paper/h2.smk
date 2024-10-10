from functools import reduce

rule compile_per_annotations:
    input:
        snps = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/matching_ids.txt",
        imd = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_snps_1",
        iga = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/iga_snps_1",
        id = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/id_snps_1",
        bmi = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/bmi_snps_1",
        edu = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/edu_snps_1"
    output:
        "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/human_default/imd_iga_id_bmi_edu_annotations.tsv.gz"
    threads: 8
    localrule: True
    #conda: env_path('pid_cfdr_pipeline.yaml')
    script: script_path('pub/igad_paper/misc/compile_per_snp_annotations.R')

rule generate_gwas_and_cfdr_annotations:
    input:
        lead_snps = "results/pub/igad_paper/tables/gwas_and_cfdr_lead_snps.tsv",
        snps = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        in_gwas = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/gwas_snps_1",
        out_gwas = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/gwas_snps_2",
        in_cfdr = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/cfdr_snps_1",
        out_cfdr = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/cfdr_snps_2"
    params:
        window = 1e5
    threads: 8
    localrule: True
    #conda: env_path('pid_cfdr_pipeline.yaml')
    script: script_path('pub/igad_paper/misc/generate_gwas_and_cfdr_annotations.R')

rule generate_imd_and_non_imd_annotations:
    input:
        imd = "results/pub/igad_paper/resources/imd_ebi_association_coordinates.tsv",
        snps = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        in_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_snps_1",
        out_snps= "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_snps_2"
    threads: 8
    localrule: True
    #conda: env_path('pid_cfdr_pipeline.yaml')
    script: script_path('pub/igad_paper/misc/generate_imd_annotations.R')

use rule generate_imd_and_non_imd_annotations as generate_id_and_non_id_annotations with:
    input:
        imd = "results/pub/igad_paper/resources/infectious_disease_ebi_association_coordinates.tsv",
        snps = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        in_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/id_snps_1",
        out_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/id_snps_2"

use rule generate_imd_and_non_imd_annotations as generate_iga_and_non_iga_annotations with:
    input:
        imd = "results/pub/igad_paper/resources/iga_association_coordinates.tsv",
        snps = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        in_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/iga_snps_1",
        out_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/iga_snps_2"

use rule generate_imd_and_non_imd_annotations as generate_bmi_and_non_bmi_annotations with:
    input:
        imd = "results/pub/igad_paper/resources/bmi_ebi_association_coordinates.tsv",
        snps = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        in_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/bmi_snps_1",
        out_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/bmi_snps_2"

use rule generate_imd_and_non_imd_annotations as generate_edu_and_non_edu_annotations with:
    input:
        imd = "results/pub/igad_paper/resources/educational_attainment_ebi_association_coordinates.tsv",
        snps = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        in_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/edu_snps_1",
        out_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/edu_snps_2"

rule generate_id_excluding_imd_iga_annotations:
    input:
        "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/human_default/imd_iga_id_bmi_edu_annotations.tsv.gz",
    output:
        in_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/id-not-imd-iga_snps_1",
        out_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/id-not-imd-iga_snps_2"
    threads: 8
    localrule: True
    #conda: env_path('pid_cfdr_pipeline.yaml')
    script: script_path('pub/igad_paper/misc/generate_id_excluding_imd_annotations.R')

rule calculate_human_default_taggings_with_two_set_partition_annotations:
    input:
        multiext("results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"),
        in_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}_snps_1",
        out_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}_snps_2",
    output :
        tagging_file = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/merged.tagging"),
    log:
        log_file = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/merged.tagging.log"
    params:
        in_stem = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/merged",
        out_stem = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/merged",
        partition_stem = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}_snps_"
    threads: 8
    resources:
        runtime = 30
    group: "sumher"
    shell:
        # NB: weightings now only used with BLD-LDAK or BLD-LDAK+Alpha models
        # NB2: thinning only required for LDAK-Thin model
        """
        ldak --calc-tagging {params.out_stem} --bfile {params.in_stem} --partition-number 2 --partition-prefix {params.partition_stem} --ignore-weights YES --window-kb 1000 --power -0.25 --max-threads {threads} > {log.log_file}
        """

use rule estimate_h2_with_human_default as estimate_h2_with_gwas_and_cfdr_gene_annotations_and_human_default with:
    input:
        gwas = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/merged.assoc",
        tagging_file = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/merged.tagging"
    output:
        multiext("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/sumher.log"
    params:
        out_stem = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/sumher",
        prevalence = lambda w: get_metadata_field('bronson-finngen-igad', 'prevalence'),
        case_prop = lambda w: get_metadata_field('bronson-finngen-igad', 'case_prop')

rule copy_imd_and_iga_annotations_to_ldak_friendly_name_format:
    input:
        imd_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_snps_1",
        iga_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/iga_snps_1"
    output:
        imd_snps = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga_snps_1"),
        iga_snps = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga_snps_2")
    localrule: True
    shell:
        """
        cp {input.imd_snps} {output.imd_snps}
        cp {input.iga_snps} {output.iga_snps}
        """

rule copy_all_category_annotations_to_ldak_friendly_name_format:
    input:
        imd_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_snps_1",
        iga_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/iga_snps_1",
        id_snps = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/id_snps_1"
    output:
        imd_snps = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats_snps_1"),
        iga_snps = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats_snps_2"),
        id_snps = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats_snps_3")
    localrule: True
    shell:
        """
        cp {input.imd_snps} {output.imd_snps}
        cp {input.iga_snps} {output.iga_snps}
        cp {input.id_snps} {output.id_snps}
        """

rule calculate_human_default_taggings_with_imd_and_iga_annotations:
    input:
        multiext("results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"),
        imd_and_iga_snps_1 = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga_snps_1",
        imd_and_iga_snps_2 = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga_snps_2"
    output :
        tagging_file = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga/human_default/merged.tagging"),
    log:
        log_file = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga/human_default/merged.tagging.log"
    params:
        in_stem = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/merged",
        out_stem = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga/human_default/merged",
        annotation_stem = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/imd_and_iga_snps_",
        annotation_number = 2
    threads: 8
    resources:
        runtime = 30
    group: "sumher"
    shell:
        # NB: weightings now only used with BLD-LDAK or BLD-LDAK+Alpha models
        # NB2: thinning only required for LDAK-Thin model
        """
        ldak --calc-tagging {params.out_stem} --bfile {params.in_stem} --annotation-number {params.annotation_number} --annotation-prefix {params.annotation_stem} --ignore-weights YES --window-kb 1000 --power -0.25 --max-threads {threads} > {log.log_file}
        """

use rule calculate_human_default_taggings_with_imd_and_iga_annotations as calculate_human_default_taggings_with_all_cats_annotations with:
    input:
        multiext("results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"),
        all_cats_snps_1 = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats_snps_1",
        all_cats_snps_2 = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats_snps_2",
        all_cats_snps_3 = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats_snps_3"
    output :
        tagging_file = temp("results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats/human_default/merged.tagging"),
    log:
        log_file = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats/human_default/merged.tagging.log"
    params:
        in_stem = "results/processed_gwas/bronson-finngen-igad/{variant_set}/{variant_type}/merged",
        out_stem = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats/human_default/merged",
        annotation_stem = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/all_cats_snps_",
        annotation_number = 3

rule compile_h2_estimates_for_annotation_category:
    input:
        hers = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/sumher.hers",
        enrich = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/sumher.enrich",
        share = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/sumher.share",
        tagging = "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/merged.tagging"
    output:
        "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/{dataset}/human_default/compiled_estimates.tsv"
    localrule: True
    run:
        hers = pd.read_csv(input.hers, sep = ' ')
        hers = hers.assign(partition = ['in', 'out', 'all'])
        hers = hers[['partition', 'Heritability', 'Her_SD']]
        hers = hers.rename({'Heritability': 'h2.obs', 'Her_SD': 'h2.obs.se'}, axis = 1)

        enrich = pd.read_csv(input.enrich, sep = ' ')
        enrich = enrich.assign(partition = ['in', 'out'])
        enrich = enrich.rename({'Enrichment': 'enrichment', 'SD.1': 'enrichment.se'}, axis = 1)
        enrich = enrich[[ 'partition', 'enrichment', 'enrichment.se']]
        enrich = pd.concat([enrich, pd.DataFrame({'partition': ['all'], 'enrichment.se': [np.NaN], 'enrichment': [np.NaN]})])

        share = pd.read_csv(input.share, sep = ' ')
        share = share.rename({'Share': 'share'}, axis = 1)
        share = share.assign(partition = ['in', 'out'])
        share = share[['partition', 'share']]
        share = pd.concat([share, pd.DataFrame({'partition': ['all'], 'share': [1]})])

        with open(input.tagging, 'r') as tagging_file:
            for line in tagging_file:
                last_line = line

        matches = re.search(r'(\d+)\s+(\d+)', last_line)
        snp_counts = pd.DataFrame({'partition': ['in', 'out', 'all'], 'SNPs': [int(matches.group(1)), int(matches.group(2)), int(matches.group(1))+int(matches.group(2))]})

        daf = reduce(lambda left, right: pd.merge(left, right, on = 'partition'), [snp_counts, hers, enrich, share])

        daf = daf.assign(dataset = wildcards['dataset'])

        daf.to_csv(output[0], index = False, sep = '\t')

rule compile_h2_estimates_across_annotation_categories:
    input:
        expand("results/pub/igad_paper/misc/ldak/{{variant_set}}/{{variant_type}}/{dataset}/human_default/compiled_estimates.tsv", dataset = ['imd', 'gwas', 'cfdr', 'id', 'iga', 'bmi', 'edu'])
    output:
        "results/pub/igad_paper/misc/ldak/{variant_set}/{variant_type}/human_default/compiled_estimates.tsv"
    localrule: True
    run:
        dafs = [pd.read_csv(x, sep = '\t') for x in input]
        pd.concat(dafs)[['dataset', 'partition', 'SNPs', 'share', 'h2.obs', 'h2.obs.se', 'enrichment', 'enrichment.se']].to_csv(output[0], sep = '\t', index = False)
