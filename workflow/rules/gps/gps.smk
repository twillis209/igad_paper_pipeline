import os
import pandas as pd

rule compute_gps_for_trait_pair:
    input:
        "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/pruned.tsv"
    output:
        "results/gps/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/gps_value.tsv"
    localrule: True
    container: "dockerhub://twillis209/gps-cpp:latest"
    shell:
        "computeGpsCLI -i {input} -a P.A -b P.B -c {wildcards.trait_A} -d {wildcards.trait_B} -n {threads} -f pp -o {output}"

rule permute_trait_pair:
    input:
        "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/pruned.tsv"
    output:
        "results/gps/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{draws}_draws/permutations.tsv"
    threads: 20
    resources:
        runtime = get_permute_time,
    group: "permute"
    container: "dockerhub://twillis209/gps-cpp:latest"
    shell:"""
        permuteTraitsCLI -i {input} -o {output} -a P.A -b P.B -c {threads} -n {wildcards.draws}

        if grep -q inf {output}; then
            exit 1
        fi
    """

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
        gps_file = "results/gps/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/gps_value.tsv",
        perm_file = "results/gps/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{draws}_draws/permutations.tsv"
    output:
        "results/gps/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{draws}_draws/pvalue.tsv"
    localrule: True
    script: script_path("gps/fit_gev_and_compute_gps_pvalue.R")

rule run_gps_on_trait_and_imds:
    input:
        [f"results/gps/{{trait}}_and_{x}/inner/{{variant_set}}/{{variant_type}}/1000kb_1_0_8/3000_draws/pvalue.tsv" for x in imd_traits]
    output:
        "results/gps/combined/{trait}/inner/{variant_set}/{variant_type}/1000kb_1_0_8/3000_draws/pvalues.tsv"
    localrule: True
    run:
        d = []

        for i, x in enumerate(input):
            trait_A, trait_B = x.split('/')[2].split('_and_')
            snp_set = x.split('/')[4]

            with open(x, 'r') as infile:
                line = infile.readline()
                line = infile.readline()
                gps = float(line.split('\t')[0])
                pvalue = float(line.split('\t')[-1])

            d.append(
                {
                    'trait_A' : trait_A,
                    'trait_B' : trait_B,
                    'snp_set' : snp_set,
                    'gps': gps,
                    'pvalue': pvalue
                }
            )

        pd.DataFrame(d).to_csv(output[0], sep = '\t', index = False)
