rule generate_all_cvid_meta_outputs:
    input:
        "results/cvid_meta/with_ukb/prescreen/gws/forest_plots/forest_plots.done",
        "results/cvid_meta/with_ukb/prescreen/gws/annotated_manhattan.distance_clumped.png",
        "results/cvid_meta/without_ukb/prescreen/gws/forest_plots/forest_plots.done",
        "results/cvid_meta/without_ukb/prescreen/gws/annotated_manhattan.distance_clumped.png"
    output:
        "results/cvid_meta/all_meta_outputs.done"
    localrule: True
    shell: "touch {output}"

rule generate_all_pad_meta_outputs:
    input:
        "results/pad_meta/with_ukb/prescreen/gws/annotated_manhattan.distance_clumped.png",
        "results/pad_meta/without_ukb/prescreen/gws/annotated_manhattan.distance_clumped.png",
    output:
        "results/pad_meta/all_meta_outputs.done"
    localrule: True
    shell: "touch {output}"

rule generate_all_igad_meta_outputs:
    input:
        "results/igad_meta/gws/annotated_manhattan.distance_clumped.png",
        "results/igad_meta/gws/forest_plots/forest_plots.done"
    output:
        "results/igad_meta/all_meta_outputs.done"
    localrule: True
    shell: "touch {output}"

rule generate_all_pid_meta_outputs:
    input:
        "results/cvid_meta/all_meta_outputs.done",
        "results/pad_meta/all_meta_outputs.done",
        "results/igad_meta/all_meta_outputs.done"
