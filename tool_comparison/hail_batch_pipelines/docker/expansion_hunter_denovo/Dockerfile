FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    curl \
    wget \
	pkg-config \
	procps \
	gcc \
	man-db \
	g++ \
	make \
	autoconf \
	unzip \
	bzip2 \
	zlib1g-dev \
	dpkg-dev \
	build-essential \
	libncurses5-dev \
	libcurl4-openssl-dev \
	libbz2-dev \
	liblzma-dev

# install ExpansionHunterDenovo binary from https://github.com/Illumina/ExpansionHunterDenovo
ENV EH_DENOVO_VERSION="v0.9.0"
RUN wget -q https://github.com/Illumina/ExpansionHunterDenovo/releases/download/${EH_DENOVO_VERSION}/ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && tar xzf ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && mv ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64/bin/ExpansionHunterDenovo /usr/local/bin/ExpansionHunterDenovo \
    && rm -r ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64 ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz

ENV SAMTOOLS_VERSION="1.9"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN python3 -m pip install str_analysis

ENTRYPOINT ["/bin/bash"]