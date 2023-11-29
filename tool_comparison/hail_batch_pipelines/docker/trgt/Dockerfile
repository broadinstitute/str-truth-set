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

RUN python3 -m pip install --upgrade cmake pip

# install the optimized version of ExpansionHunter which runs 2 to 3x faster and produces the same output
# wget -q https://github.com/PacificBiosciences/trgt/releases/download/v0.5.0/trvz-v0.5.0-linux_x86_64.gz \

RUN wget -q https://github.com/PacificBiosciences/trgt/releases/download/v0.5.0/trgt-v0.5.0-linux_x86_64.gz \
    && wget -q https://github.com/PacificBiosciences/trgt/releases/download/v0.5.0/trvz-v0.5.0-linux_x86_64.gz \
	&& gunzip trgt-v0.5.0-linux_x86_64.gz trvz-v0.5.0-linux_x86_64.gz \
    && chmod 777 trgt-v0.5.0-linux_x86_64 trvz-v0.5.0-linux_x86_64 \
    && mv trgt-v0.5.0-linux_x86_64 /usr/bin/trgt \
    && mv trvz-v0.5.0-linux_x86_64 /usr/bin/trvz

# install gcloud
ENV GCLOUD_SDK_VERISON="419.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

RUN python3 -m pip install str_analysis==1.2.5

ENTRYPOINT ["/bin/bash"]