FROM python:3.11-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    pkg-config \
    dpkg-dev \
    build-essential \
    libssl-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libdeflate-dev \
    gcc \
    g++ \
    autoconf \
    make \
    unzip \
    bzip2 \
    procps \
    curl \
    wget \
    git \
    jq \
    time

# htslib + samtools provide bgzip/tabix (index the regions bed) and samtools (decode CRAM -> BAM before genotyping).
ENV SAMTOOLS_VERSION="1.20"
RUN wget -q https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd htslib-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    && make clean

RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    && make clean

RUN python3 -m pip install --no-cache-dir --upgrade pip

# ATaRVa the tandem repeat genotyper. Installed from the bw2 fork at a pinned commit, which adds the
# reference_filename fix so CRAM reads decode reliably (SowpatiLab/ATaRVa#13); this lets the pipeline genotype CRAM
# input directly instead of converting it to BAM first. The prebuilt libssw.so ships in the repo, so no compile step
# is needed. Pinned to a commit so genotyping stays reproducible across image rebuilds.
RUN python3 -m pip install --no-cache-dir \
    'git+https://github.com/bw2/ATaRVa.git@3669f10aaab91b222388fcee6d38f235890174e8'

# str-analysis provides convert_atarva_vcf_to_expansion_hunter_json + combine_str_json_to_tsv. The atarva_pipeline
# runtime-installs a pinned commit (STR_ANALYSIS_COMMIT) on top of this so the converter can be updated without
# rebuilding the image; this baked copy is the fallback and pre-warms the dependencies.
RUN python3 -m pip install --no-cache-dir \
    'git+https://github.com/broadinstitute/str-analysis@3d5e3dc37161d41a1bc6f92afdcf3fc99e81b30a'

ENTRYPOINT ["/bin/bash"]
