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

RUN python3 -m pip install --no-cache-dir --upgrade cmake pip

RUN apt update && apt-get install --no-install-recommends -y \
	libboost-all-dev

RUN git clone https://github.com/ChaissonLab/vamos.git \
	&& cd vamos*/src/ \
	&& git checkout 1.3.6 \
	&& make \
	&& cp vamos /usr/bin/

# download the original, q10, q20, and q30 motif sets
RUN cd / \
	&& curl "https://zenodo.org/record/8357361/files/original_motifs.set148.bed.gz?download=1" > original_motifs.set148.bed.gz \
    && curl "https://zenodo.org/record/8357361/files/q-0.1_motifs.set148.bed.gz?download=1"  > q-0.1_motifs.set148.bed.gz \
	&& curl "https://zenodo.org/record/8357361/files/q-0.2_motifs.set148.bed.gz?download=1"  > q-0.2_motifs.set148.bed.gz \
	&& curl "https://zenodo.org/record/8357361/files/q-0.3_motifs.set148.bed.gz?download=1"  > q-0.3_motifs.set148.bed.gz


ENTRYPOINT ["/bin/bash"]