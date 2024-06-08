rule get_ebi_associations_for_imd:
    output:
        "results/pub/igad_paper/resources/{trait,!bmi}_ebi_associations.tsv"
    params:
        efo_trait_code = lambda w: config.get('igad_paper').get('imd_trait_efo_codes').get(w.trait)
    resources:
        runtime = 10
    group: "igad_paper"
    conda: env_path("gwasrapidd.yaml")
    script: script_path('pub/igad_paper/misc/fetch_ebi_associations_for_imd.R')

#use rule get_ebi_associations_for_imd as get_ebi_assocations_for_bmi with:
#    output:
#        "results/pub/igad_paper/resources/bmi_ebi_associations.tsv"
#    params:
#        efo_trait_code = lambda w: config.get('igad_paper').get('bmi_trait_efo_code')
#    resources:
#        runtime = 20

rule create_ebi_associations_intervals_for_imd_panel:
    input:
        [f"results/pub/igad_paper/resources/{trait}_ebi_associations.tsv" for trait in ["asthma-ex", "hyperpara", "hypothy", "derm-ecz", "pbc", "addi", "psc", "ra", "jia", "sle", "crohns", "t1d", "ms", "uc-delange", "igan"]]
    output:
        "results/pub/igad_paper/resources/imd_ebi_association_coordinates.tsv"
    params:
        window = 1e5
    threads: 8
    localrule: True
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/misc/compile_imd_ebi_associations.R')

rule create_ebi_associations_intervals_for_infectious_disease:
    input:
        "results/pub/igad_paper/resources/infectious_disease_ebi_associations.tsv"
    output:
        "results/pub/igad_paper/resources/infectious_disease_ebi_association_coordinates.tsv"
    params:
        window = 1e5
    localrule: True
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/misc/create_ebi_associations_intervals_for_infectious_disease.R')

rule create_ebi_associations_intervals_for_educational_attainment:
    input:
        "results/pub/igad_paper/resources/educational_attainment_ebi_associations.tsv"
    output:
        "results/pub/igad_paper/resources/educational_attainment_ebi_association_coordinates.tsv"
    params:
        window = 1e5
    localrule: True
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/misc/create_ebi_associations_intervals_for_edu_attainment.R')

rule create_ebi_associations_intervals_for_bmi:
    input:
        "results/pub/igad_paper/resources/bmi_ebi_associations.tsv"
    output:
        "results/pub/igad_paper/resources/bmi_ebi_association_coordinates.tsv"
    params:
        window = 1e5
    localrule: True
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/misc/create_ebi_associations_intervals_for_bmi.R')

rule create_associations_intervals_for_iga:
    input:
        "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps_and_existing_associations.tsv"
    output:
        "results/pub/igad_paper/resources/iga_association_coordinates.tsv"
    params:
        window = 1e5
    localrule: True
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/misc/create_association_intervals_for_iga.R')
