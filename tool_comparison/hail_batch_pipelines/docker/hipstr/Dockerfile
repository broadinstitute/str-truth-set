FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    curl \
    wget \
    git \
    time \
    autoconf \
    automake \
    libtool \
    libbz2-dev \
	zlib1g-dev \
	liblzma-dev \
	libssl-dev \
	libcurl4-openssl-dev \
	build-essential

RUN git clone https://github.com/HipSTR-Tool/HipSTR \
    && cd HipSTR \
    && make \
    && cp /HipSTR/HipSTR /bin/HipSTR \
    && make clean

ENV SAMTOOLS_VERSION="1.9"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN python3 -m pip install str_analysis==1.2.5

ENTRYPOINT ["/bin/bash"]
