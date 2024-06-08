rule run_mr_on_iga_and_igad:
    input:
        "results/pub/igad_paper/misc/stats_at_gws_iga_snps.tsv"
    output:
        tsv = "results/pub/igad_paper/misc/mr/igad_and_iga_mr.tsv",
        png = "results/pub/igad_paper/misc/mr/igad_and_iga.png",
        png_ivw = "results/pub/igad_paper/misc/mr/igad_and_iga_ivw.png"
    params:
        x_beta = 'BETA.meta',
        x_se = 'SE.meta',
        y_beta = 'BETA.igad',
        y_se = 'SE.igad',
        snp_col = 'rsID',
        snps_to_exclude = ['rs16830188']
    localrule: True
    container: "docker://twillis209/r-mendelianrandomization"
    script: script_path("pub/igad_paper/misc/run_mr.R")

use rule run_mr_on_iga_and_igad as run_mr_on_iga_and_igan with:
    output:
        tsv = "results/pub/igad_paper/misc/mr/igan_and_iga_mr.tsv",
        png = "results/pub/igad_paper/misc/mr/igan_and_iga.png",
        png_ivw = "results/pub/igad_paper/misc/mr/igan_and_iga_ivw.png"
    params:
        x_beta = 'BETA.meta',
        x_se = 'SE.meta',
        y_beta = 'BETA.igan',
        y_se = 'SE.igan',
        snp_col = 'rsID',
        snps_to_exclude = ['rs3181356']
    localrule: True
