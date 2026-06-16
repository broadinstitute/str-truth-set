FROM python:3.9-slim-bullseye

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

# pin cmake < 4: vamos vendors parasail-2.6.2, whose CMakeLists declares cmake_minimum_required < 3.5, which
# CMake >= 4 refuses to configure ("Compatibility with CMake < 3.5 has been removed")
RUN python3 -m pip install --no-cache-dir --upgrade pip "cmake<4"

RUN apt update && apt-get install --no-install-recommends -y \
	libboost-all-dev

RUN git clone https://github.com/bw2/vamos.git \
	&& cd vamos \
	&& git checkout 5c4a069847fbe6ca6067de1bff6098c9403db9ec \
	&& cd src/ \
	&& make \
	&& cp vamos /usr/bin/

# download the original, q10, q20, and q30 motif sets
RUN cd / \
	&& curl -fsSL "https://zenodo.org/record/8357361/files/original_motifs.set148.bed.gz?download=1" > original_motifs.set148.bed.gz \
    && curl -fsSL "https://zenodo.org/record/8357361/files/q-0.1_motifs.set148.bed.gz?download=1"  > q-0.1_motifs.set148.bed.gz \
	&& curl -fsSL "https://zenodo.org/record/8357361/files/q-0.2_motifs.set148.bed.gz?download=1"  > q-0.2_motifs.set148.bed.gz \
	&& curl -fsSL "https://zenodo.org/record/8357361/files/q-0.3_motifs.set148.bed.gz?download=1"  > q-0.3_motifs.set148.bed.gz


ENTRYPOINT ["/bin/bash"]