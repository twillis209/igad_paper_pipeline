rule download_tangye_table:
    output:
        "resources/pid_genes/tangye_table.xlsx"
    localrule: True
    shell:
        """
        wget -O {output} ncbi.nlm.nih.gov/pmc/articles/PMC9244088/bin/10875_2022_1289_MOESM2_ESM.xlsx
        """
