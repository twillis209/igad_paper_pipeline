rule plot_one_sample_qqplot_for_pad_meta_on_pid_and_non_pid:
    input:
        "results/pid/{window}kb/pad_with_pid_and_non_pid.tsv.gz"
    output:
        one_sample = "results/pid/{window}kb/qqplots/{snp_set}/pad_with_pid_and_non_pid_one_sample.png"
    threads: 8
    resources:
        runtime = 10
    script: script_path("gene_centric/pid/plot_one_sample_qqplot.R")

use rule plot_one_sample_qqplot_for_pad_meta_on_pid_and_non_pid as plot_one_sample_qqplot_for_igad_meta_on_pid_and_non_pid with:
    input:
        "results/pid/{window}kb/bronson-finngen-igad_with_pid_and_non_pid.tsv.gz"
    output:
        one_sample = "results/pid/{window}kb/qqplots/{snp_set}/bronson-finngen-igad_with_pid_and_non_pid_one_sample.png"

rule bootstrap_pi0_estimate:
    input:
        "results/pid/{window}kb/{dataset}_with_pid_and_non_pid.tsv.gz"
    output:
        "results/pid/{window}kb/{dataset}_with_pid_and_non_pid_pi0_boot.rds"
    threads: 8
    resources:
        runtime = 20
    script: script_path("gene_centric/pid/bootstrap_pi0_estimate.R")

rule run_auc_permutation_test_for_snp_samples:
    input:
        "results/pid/{window}kb/{dataset}_with_pid_and_non_pid.tsv.gz"
    output:
        "results/pid/{window}kb/{dataset}_with_pid_and_non_pid_auc.rds"
    threads: 8
    resources:
        runtime = 20
    script: script_path("gene_centric/pid/auc_permutation_test.R")

rule run_auc_permutation_for_gene_samples:
    input:
        all_genes = "results/pid/{window}kb_all_genes.rds",
        gwas = "results/igad_meta/meta.tsv.gz"
    output:
        "results/pid/{window}kb_{permutations}_{iei_to_all_ratio}/{seed}_auc_permutations.tsv.gz"
    params:
        no_of_permutations = lambda w: int(w.permutations),
        no_of_genes = lambda w: int(w.iei_to_all_ratio) * 448,
        seed = lambda w: int(w.seed)
    threads: 10
    resources:
        runtime = 20
    script: script_path("gene_centric/pid/auc_permutation_for_gene_samples.R")

rule run_gif_permutations_for_gene_samples:
    input:
        all_genes = "results/pid/{window}kb_all_genes.rds",
        gwas = "results/igad_meta/meta.tsv.gz"
    output:
        "results/pid/{window}kb_{permutations}_{iei_to_all_ratio}/{seed}_gif_permutations.tsv.gz"
    params:
        no_of_permutations = lambda w: int(w.permutations),
        no_of_genes = lambda w: int(w.iei_to_all_ratio) * 448,
        seed = lambda w: int(w.seed)
    threads: 10
    resources:
        runtime = 20
    script: script_path("gene_centric/pid/gif_permutation_for_gene_samples.R")
