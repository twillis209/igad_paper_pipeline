# `snakemake` pipeline for the paper 'Leveraging pleiotropy identifies common variant associations with selective IgA deficiency'

This repository contains a `snakemake` pipeline to generate the results presented in the paper ['Leveraging pleiotropy identifies common variant associations with selective IgA deficiency'](https://github.com/twillis209/igad_paper_pipeline). It has been tested with `v8.23` of `snakemake`.

I've used `conda` and `docker` to try to make things as reproducible as possible. `snakemake` should pull the necessary dependencies and images managed by these two tools, but note that one can't manage installation of `snakemake` itself within the pipeline: I recommend `conda` to manage `snakemake` alongside the pipeline.

Note that the `master` branch contains the `docker`-enabled version; the `conda` branch contains a version in which the workflow is managed predominantly by `conda`. Each version uses `conda` and `docker` to some extent, e.g. [the `locuszoomr` package](https://github.com/myles-lewis/locuszoomr) is handled by `docker` in both using an image I maintain (owing to difficulties with getting `conda` to manage the R package's dependencies), whilst `ldsc` and its Python 2 environment is handled by `conda` (I did not want to add a `conda` environment to the `docker` image).

I think `docker` is preferable for managing workflows like this, but I maintain the `conda` version (and in fact used an earlier version of the same to produce the results in the paper) as my current HPC environment does not play nicely with the way `snakemake` uses `apptainer` ([see this GitHub issue if you'd like to know more](https://github.com/snakemake/snakemake/issues/2959)).

If you are using the `conda` version, I recommend running `snakemake -n all_igad_paper_artifacts` first to preview the workflow, then `snakemake --conda-create-envs-only all_igad_paper_artifacts` so that the prerequisite `conda` environments are built first. This can be a quite a lengthy process and `snakemake` doesn't log much output whilst it is happening. Note also that `ldak` is not provided by `conda`, you will have to install it and place it on the path.

## Reproduction

The `smk` files in the directory `workflow/rules/pub/igad_paper` contain the `snakemake` rules which generate outputs presented in the paper, such as figures (`figures.smk`), tables (`tables.smk`), and data sets. Targeting `snakemake` at the rule `all_igad_paper_artifacts` in `workflow/rules/pub/igad_paper/export_data.smk` should have it generate all the results in the paper, although *I warn you* to use `snakemake -n` or `snakemake --dry-run` to get a preview of just how much processing this will take; it won't be cheap in terms of CPU hours!

Some potential issues with this pipeline:
* I used `plink v2.00a6LM` in the paper and whilst this version is used in the `docker` image, the version made available through Bioconda (and used in the `conda` branch) is (as of 22/10/24) `v2.00a5.12`
  * The `docker` image contains LDAK `v6` wheres LDAK `v5` was used to produce the results in the paper

To my knowledge, neither of these versioning discrepancies should affect the results of this workflow in any significant way.

As this pipeline has been spun out of a larger one containing work relating to several projects, it is possible that I have introduced regressions when downsizing repo content. Please open an issue if you encounter an error and I'll do my best to help you.
