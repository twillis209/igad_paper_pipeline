# README

## Setup

This repository contains a `snakemake` pipeline to generate the results presented in our paper 'Leveraging pleiotropy identifies common variant associations with selective IgA deficiency'. It has been tested with `v8.23` of `snakemake`.

We've used `conda` and `docker` to try to make things as reproducible as possible. `snakemake` should pull the necessary dependencies and images managed by these two tools. We can't manage installation of `snakemake` itself within the pipeline. I recommend `conda` to manage `snakemake` alongside the pipeline.

As noted in the manuscript, we used `LDAK` (`v5.2`), which can be downloaded from Doug Speed's website [here](https://dougspeed.com/downloads). The code in the pipeline expects `LDAK` to be on the `PATH`.

## Reproduction

The `smk` files in the directory `workflow/rules/pub/igad_paper` contain the `snakemake` rules which generate outputs presented in the paper, such as figures (`figures.smk`), tables (`tables.smk`), and data sets.

As this pipeline has been spun out of a larger one containing work relating to several projects, it is possible that I have introduced regressions when downsizing repo content. Please open an issue if you encounter an error and I'll do my best to help you.
