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
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
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

ENV STRAGLR_VERSION=v1.5.0
RUN python3 -m pip install --no-cache-dir git+https://github.com/bcgsc/straglr.git@${STRAGLR_VERSION}#egg=straglr

#https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
ENV BLAST_VERSION="2.16.0+"
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${BLAST_VERSION}-x64-linux.tar.gz \
	&& tar xzf ncbi-blast-${BLAST_VERSION}-x64-linux.tar.gz \
	&& rm ncbi-blast-${BLAST_VERSION}-x64-linux.tar.gz \
	&& mv ncbi-blast-2.16.0+/bin/blastn /usr/bin/blastn \
    && rm  -rf ncbi-blast-2.16.0+ \
    && chmod 777 /usr/bin/blastn

RUN wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 \
	&& chmod 777 trf409.linux64 \
	&& mv trf409.linux64 /usr/bin/trf

ENTRYPOINT ["/bin/bash"]