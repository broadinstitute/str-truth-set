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

RUN python3 -m pip install --upgrade cmake pip

RUN git clone https://github.com/bw2/ExpansionHunter.git \
    && cd ExpansionHunter \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && make clean

RUN python3 -m pip install str_analysis

RUN cp /ExpansionHunter/build/install/bin/ExpansionHunter /bin/

ENTRYPOINT ["/bin/bash"]