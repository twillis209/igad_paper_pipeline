rule subset_summary_statistics_about_variant_for_igad_paper:
    input:
        "results/pub/igad_paper/misc/merged_meta_cfdr_and_aux_data.tsv.gz"
    output:
        sum_stats = temp("results/pub/igad_paper/misc/{variant_id}/{window_size}/sum_stats_with_maf.tsv.gz"),
        ids = temp("results/pub/igad_paper/misc/{variant_id}/{window_size}/ids.txt")
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        maf_col = 'ALT_FREQS',
        beta_cols = ['BETA.igad', 'BETA.ra', 'BETA.asthma', 'BETA.iga', 'BETA.igan'],
        se_cols = ['SE.igad', 'SE.ra', 'SE.asthma', 'SE.iga', 'SE.igan'],
        p_cols = ['P.igad', 'P.cfdr', 'P.ra', 'P.asthma', 'P.iga', 'P.igan'],
        variant_id = lambda w: w.variant_id.replace('_', ':') if '_' in w.variant_id else w.variant_id,
        window = lambda w: int(w.window_size.replace('kb', '')) * 1000
    threads: 8
    resources:
    script: script_path("pub/igad_paper/misc/subset_summary_statistics_about_variant_for_igad_paper.R")

def map_igad_paper_abbrv_to_metadata_abbrv(paper_abbrv):
    return {'igad': 'bronson-finngen-igad',
                  'igan': 'igan',
                  'asthma': 'asthma-ex',
                  'ra': 'ra-ishigaki',
                  'iga': 'liu-decode-lyons-dennis-iga'}[paper_abbrv]

rule run_coloc_for_two_traits_for_igad_paper:
    input:
        "results/pub/igad_paper/misc/{variant_id}/{window_size}/sum_stats_with_maf.tsv.gz"
    output:
        "results/pub/igad_paper/misc/{variant_id}/{window_size}/{first_trait,igad|ra|iga|igan|asthma}_{second_trait,igad|ra|iga|igan|asthma}_coloc.rds"
    params:
        first_controls = lambda w: int(get_metadata_field(map_igad_paper_abbrv_to_metadata_abbrv(w.first_trait), 'N0')),
        first_cases = lambda w: int(get_metadata_field(map_igad_paper_abbrv_to_metadata_abbrv(w.first_trait), 'N1')),
        first_beta = lambda w: f'BETA.{w.first_trait}',
        first_se = lambda w: f'SE.{w.first_trait}',
        second_controls = lambda w: int(get_metadata_field(map_igad_paper_abbrv_to_metadata_abbrv(w.second_trait), 'N0')),
        second_cases = lambda w: int(get_metadata_field(map_igad_paper_abbrv_to_metadata_abbrv(w.second_trait), 'N1')),
        second_beta = lambda w: f'BETA.{w.second_trait}',
        second_se = lambda w: f'SE.{w.second_trait}',
        maf_col = 'ALT_FREQS',
    localrule: True
    conda: env_path("coloc.yaml")
    script: script_path("pub/igad_paper/misc/run_coloc_for_igad_and_aux_trait.R")

rule run_coloc_at_runx3_locus:
    input:
        "results/pub/igad_paper/misc/1_24972350_C_T/100kb/igad_iga_coloc.rds",
        "results/pub/igad_paper/misc/1_24972350_C_T/100kb/igad_ra_coloc.rds",
        "results/pub/igad_paper/misc/1_24972350_C_T/100kb/igad_asthma_coloc.rds"

rule run_coloc_for_trait_and_iga_at_all_iga_gws_hits:
    input:
        [f"results/pub/igad_paper/misc/{y.replace(':', '_')}/100kb/{{trait}}_iga_coloc.rds" for y in [config.get('igad_paper').get('iga_loci').get(x).get('index_snp') for x in config.get('igad_paper').get('iga_loci')]]
    output:
        "results/pub/igad_paper/misc/coloc/{trait,igad|igan}_and_iga_at_gws_iga_hits.tsv"
    params:
        genes = config.get('igad_paper').get('iga_loci')
    localrule: True
    script: script_path("pub/igad_paper/misc/compile_coloc_results_at_index_snps.R")

use rule run_coloc_for_trait_and_iga_at_all_iga_gws_hits as run_coloc_for_igad_and_aux_traits_at_all_cfdr_hits with:
    input:
        [f"results/pub/igad_paper/misc/{y.replace(':', '_')}/100kb/igad_{{trait}}_coloc.rds" for y in [config.get('igad_paper').get('cfdr_loci').get(x).get('index_snp') for x in config.get('igad_paper').get('cfdr_loci')]]
    output:
        temp("results/pub/igad_paper/misc/coloc/igad_and_{trait,iga|asthma|ra}_at_cfdr_hits.tsv")
    params:
        genes = config.get('igad_paper').get('cfdr_loci')
    localrule: True

rule compile_igad_and_aux_trait_coloc_results:
    input:
        iga = "results/pub/igad_paper/misc/coloc/igad_and_iga_at_cfdr_hits.tsv",
        ra = "results/pub/igad_paper/misc/coloc/igad_and_ra_at_cfdr_hits.tsv",
        asthma = "results/pub/igad_paper/misc/coloc/igad_and_asthma_at_cfdr_hits.tsv"
    output:
        "results/pub/igad_paper/misc/coloc/igad_and_aux_traits_at_cfdr_hits.tsv"
    localrule: True
    run:
        iga = pd.read_csv(input.iga, sep = '\t')
        iga['aux_trait'] = 'iga'

        ra = pd.read_csv(input.ra, sep = '\t')
        ra['aux_trait'] = 'ra'

        asthma = pd.read_csv(input.asthma, sep = '\t')
        asthma['aux_trait'] = 'asthma'

        merged = pd.concat([iga, ra, asthma])

        merged.to_csv(output[0], sep = '\t', index = False)
