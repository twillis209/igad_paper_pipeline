rule draw_igad_locus_plot_for_igad_paper:
    input:
        "results/pub/igad_paper/misc/merged_meta_cfdr_and_aux_data.tsv.gz"
    output:
        "results/pub/igad_paper/figures/locus_plots/igad/{locus}/index_snp.png"
    params:
        index_snp_seqname = lambda w: int(config.get('igad_paper').get('igad_loci').get(w.locus).get('index_snp').split(':')[0]),
        index_snp_pos = lambda w: int(config.get('igad_paper').get('igad_loci').get(w.locus).get('index_snp').split(':')[1]),
        flank = lambda w: int(config.get('igad_paper').get('igad_loci').get(w.locus).get('flank').replace('kb', '')) * 1000
    threads: 2
    resources:
        runtime = 8
    #container: "docker://twillis209/r-locuszoomr:latest"
    script: script_path("pub/igad_paper/figures/draw_locus_plot_for_igad_paper.R")

use rule draw_igad_locus_plot_for_igad_paper as draw_iga_locus_plot_for_igad_paper with:
    output:
        "results/pub/igad_paper/figures/locus_plots/iga/{locus}/index_snp.png"
    params:
        index_snp_seqname = lambda w: int(config.get('igad_paper').get('iga_loci').get(w.locus).get('index_snp').split(':')[0]),
        index_snp_pos = lambda w: int(config.get('igad_paper').get('iga_loci').get(w.locus).get('index_snp').split(':')[1]),
        flank = lambda w: int(config.get('igad_paper').get('iga_loci').get(w.locus).get('flank').replace('kb', '')) * 1000

use rule draw_igad_locus_plot_for_igad_paper as draw_cfdr_locus_plot_for_igad_paper with:
    output:
        "results/pub/igad_paper/figures/locus_plots/cfdr/{locus}/index_snp.png"
    params:
        index_snp_seqname = lambda w: int(config.get('igad_paper').get('cfdr_loci').get(w.locus).get('index_snp').split(':')[0]),
        index_snp_pos = lambda w: int(config.get('igad_paper').get('cfdr_loci').get(w.locus).get('index_snp').split(':')[1]),
        flank = lambda w: int(config.get('igad_paper').get('cfdr_loci').get(w.locus).get('flank').replace('kb', '')) * 1000

rule draw_all_locus_plots_for_igad_paper:
    input:
        [f"results/pub/igad_paper/figures/locus_plots/igad/{x}/index_snp.png" for x in config.get('igad_paper').get('igad_loci')],
        [f"results/pub/igad_paper/figures/locus_plots/iga/{x}/index_snp.png" for x in config.get('igad_paper').get('iga_loci')],
        [f"results/pub/igad_paper/figures/locus_plots/cfdr/{x}/index_snp.png" for x in config.get('igad_paper').get('cfdr_loci')]

use rule draw_iga_locus_plot_for_igad_paper as draw_iga_locus_plot_for_igad_paper_with_yfloor with:
    output:
        "results/pub/igad_paper/figures/locus_plots/iga/{locus}/{floor}/index_snp.png"
    params:
        index_snp_seqname = lambda w: int(config.get('igad_paper').get('iga_loci').get(w.locus).get('index_snp').split(':')[0]),
        index_snp_pos = lambda w: int(config.get('igad_paper').get('iga_loci').get(w.locus).get('index_snp').split(':')[1]),
        flank = lambda w: int(config.get('igad_paper').get('iga_loci').get(w.locus).get('flank').replace('kb', '')) * 1000,
        yfloor = lambda w: pow(10, -int(w.floor))
