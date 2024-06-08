import os
import pandas as pd
from itertools import chain
from scipy.stats import chi2

rule enumerate_imd_data_sets:
    input:
        [f"results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/realigned_alleles/{x}_table.tsv.gz" for x in imd_traits]

rule run_gps_on_cvid_meta_analysis_and_imds:
    input:
        "results/gps/combined/10kG-finngen-li-ukb-cvid/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalues.tsv"

rule run_gps_on_pad_meta_analysis_and_imds:
    input:
        "results/gps/combined/10kG-finngen-li-bronson-ukb-pad/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalues.tsv",
        "results/ldak/ldak-thin/combined/sans_mhc/snps_only/10kG-finngen-li-bronson-ukb-pad_and_imds.tsv"

rule run_gps_on_igad_meta_analysis_and_imds:
    input:
        "results/gps/combined/bronson-finngen-igad/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalues.tsv"

rule rerun_gps_on_pid_for_new_iga:
    input:
        "results/gps/10kG-finngen-li-ukb-cvid_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalue.tsv",
        "results/gps/bronson-finngen-igad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalue.tsv",
        "results/gps/10kG-finngen-li-bronson-ukb-pad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalue.tsv"

