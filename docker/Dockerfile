FROM rocker/r-ver:latest

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    cmake \
    git \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error \
    -n 4 \
    data.table \
    R.utils \
    ggplot2 \
    remotes \
    gridExtra \
    devtools \
    coloc \
    evd \
    fitdistrplus \
    ggcorrplot \
    ggnewscale \
    ggrepel \
    ggtext \
    gridExtra \
    httr \
    kableExtra \
    LDlinkR \
    tidyverse \
    parallel \
    patchwork \
    pheatmap \
    MendelianRandomization \
    pracma \
    stringr \
    BiocManager

# Using -d flag with this caused a subscript out of bounds error
RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r --error \
    -n 4 \
    ensembldb \
    GenomicFeatures \
    GenomicRanges \
    EnsDb.Hsapiens.v75 \
    EnsDb.Hsapiens.v86 \
    org.Hs.eg.db \
    TxDb.Hsapiens.UCSC.hg38.knownGene

RUN R -e "remotes::install_github('myles-lewis/locuszoomr')"
RUN R -e "remotes::install_github('twillis209/cfdr')"

RUN wget -O /usr/local/bin/plink2.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20241007.zip
RUN unzip /usr/local/bin/plink2.zip -d /usr/local/bin

RUN wget -O /usr/local/bin/ldak https://github.com/dougspeed/LDAK/raw/main/ldak6.linux
RUN chmod u+x /usr/local/bin/ldak

ENTRYPOINT ["/bin/bash"]