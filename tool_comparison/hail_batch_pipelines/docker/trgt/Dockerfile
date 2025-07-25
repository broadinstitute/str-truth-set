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

ENV GCLOUD_SDK_VERISON="481.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV TRGT_VERSION="3.0.0"
RUN wget -q https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
    && tar xzf trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
    && chmod 777 trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt \
    && mv trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt /usr/bin/trgt


ENV PATH=/google-cloud-sdk/bin:$PATH

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN python3 -m pip install --no-cache-dir git+https://github.com/broadinstitute/str-analysis

ENTRYPOINT ["/bin/bash"]