# README

This repository contains a `snakemake` pipeline to generate the results presented in our paper 'Leveraging pleiotropy identifies common variant associations with selective IgA deficiency'. It has been tested with `v8.11.3` of `snakemake`. 

We've used `conda` and `docker` to try to make things as reproducible as possible. `snakemake` should pull the necessary dependencies and images managed by these two tools. It's worth noting that building the larger `conda` environments, such as the `pid_cfdr_pipeline` one, can take over an hour.

When cloning this repo, use the `--recurse-submodules` flag to ensure the `gps_cpp` submodule is also cloned: this submodule contains the code for the GPS test software, which can be built by following the instructions [here](github.com/twillis209/gps_cpp).

As noted in the manuscript, we used `LDAK` (`v5.2`), which can be downloaded from Doug Speed's website [here](https://dougspeed.com/downloads2). The code in the pipeline expects `LDAK` to be on the `PATH`.

The `smk` files in the directory `workflow/rules/pub/igad_paper` contain the `snakemake` rules which generate outputs presented in the paper, such as figures (`figures.smk`, tables (`tables.smk`), and data sets.

As this pipeline has been spun out of a larger one containing work relating to several projects, it is possible that I have introduced regressions when downsizing repo content. Please open an issue if you encounter an error and I'll do my best to help you.
