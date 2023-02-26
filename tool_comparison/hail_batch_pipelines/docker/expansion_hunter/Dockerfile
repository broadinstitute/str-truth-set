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

ENV SAMTOOLS_VERSION="1.17"
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

RUN python3 -m pip install --upgrade cmake pip

# install the optimized version of ExpansionHunter which runs 2 to 3x faster and produces the same output
RUN git clone https://github.com/bw2/ExpansionHunter.git \
    && cd ExpansionHunter \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && make clean

RUN cp /ExpansionHunter/build/install/bin/ExpansionHunter /bin/

# Also, install the original, unmodified version of ExpansionHunter
RUN git clone https://github.com/Illumina/ExpansionHunter.git IlluminaExpansionHunter \
    && cd IlluminaExpansionHunter \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && make clean

RUN cp /IlluminaExpansionHunter/build/install/bin/ExpansionHunter /bin/IlluminaExpansionHunter

# install REViewer
RUN git clone https://github.com/bw2/REViewer \
    && cd REViewer  \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && make clean

RUN python3 -m pip install str_analysis==0.9.15

RUN apt-get install -y jq

ENTRYPOINT ["/bin/bash"]