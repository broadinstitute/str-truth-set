FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    curl \
    wget \
	pkg-config \
	procps \
    git \
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

ENV SAMTOOLS_VERSION="1.9"
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

RUN apt update && apt-get install --no-install-recommends -y libboost-all-dev

RUN wget https://github.com/Illumina/ExpansionHunter/archive/refs/heads/release-2.5.5.zip \
	&& unzip release-2.5.5.zip \
	&& rm release-2.5.5.zip \
    && cd ExpansionHunter-release-2.5.5 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && cp /ExpansionHunter-release-2.5.5/build/ExpansionHunter /bin/ \
    && make clean

RUN python3 -m pip install str_analysis==0.9.15

ENTRYPOINT ["/bin/bash"]