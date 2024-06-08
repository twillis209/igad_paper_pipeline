import os
import pandas as pd

rule compute_gps_for_trait_pair:
    input:
        "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/pruned.tsv"
    output:
        "results/gps/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/gps_value.tsv"
    localrule: True
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c {wildcards.trait_A} -d {wildcards.trait_B} -n {threads} -f pp -o {output}"

rule permute_trait_pair:
    input:
        "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/pruned.tsv"
    output:
        "results/gps/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{draws}_draws/permutations.tsv"
    threads: 20
    resources:
        runtime = get_permute_time,
    group: "permute"
    shell:"""
        workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input} -o {output} -a P.A -b P.B -c {threads} -n {wildcards.draws}

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
    conda: env_path("pid_cfdr_pipeline.yaml")
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

rule compute_qvalues_for_gps_on_trait_and_imds:
    input:
        "results/gps/combined/{trait}/inner/sans_mhc/all/1000kb_1_0_8/3000_draws/pvalues.tsv"
    output:
        "results/gps/combined/{trait}/inner/sans_mhc/all/1000kb_1_0_8/3000_draws/qvalues.tsv"
    localrule: True
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gps/compute_qvalues_for_gps_pvalues.R")


rule run_pid_against_updated_iga_with_gps:
    input:
        "results/gps/bronson-finngen-igad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalue.tsv",
        "results/gps/10kG-finngen-li-bronson-ukb-pad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalue.tsv",
        "results/gps/10kG-finngen-li-ukb-cvid_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalue.tsv"
