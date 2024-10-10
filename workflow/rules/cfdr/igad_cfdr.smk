rule drop_spurious_igad_cfdr_hits:
    input:
        "results/cfdr/bronson-finngen-igad-gc_on_asthma-ex_and_ra-ishigaki_and_liu-decode-lyons-dennis-iga/left/sans_mhc/snps_only/1000kb_1_0_2/3/cfdr.tsv.gz"
    output:
        "results/cfdr/bronson-finngen-igad-gc_on_asthma-ex_and_ra-ishigaki_and_liu-decode-lyons-dennis-iga/left/sans_mhc/snps_only/1000kb_1_0_2/3/postscreen/cfdr.tsv.gz"
    params:
        snp_col = 'SNPID',
        snps_to_drop = ['3:119733858:A:G']
    threads: 16
    resources:
        runtime = 20
    group: "gwas"
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/drop_spurious_hits.R")
