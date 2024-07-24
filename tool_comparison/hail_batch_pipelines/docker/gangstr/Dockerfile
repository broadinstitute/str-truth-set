FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    curl \
    wget \
    git \
    autoconf \
    automake \
    libtool \
    libbz2-dev \
    time \
	zlib1g-dev \
	liblzma-dev \
	libssl-dev \
	libcurl4-openssl-dev \
	build-essential


ENV GANGSTR_VERSION="2.5.0-Source"
ENV HTSLIB_CFLAGS="-I./htslib-1.9"
ENV HTSLIB_LIBS=" -L$./htslib-1.9 -lhts"
ENV NLOPT_CFLAGS="-I./nlopt-2.4.2"
ENV NLOPT_LIBS=" -L$./nlopt-2.4.2 -lnlopt"
ENV CPPFLAGS="-std=c++98"

ENV LD_LIBRARY_PATH=/usr/local/lib/

RUN wget -q https://github.com/gymreklab/GangSTR/releases/download/v2.5/GangSTR-${GANGSTR_VERSION}.tar.gz \
    && tar xzf GangSTR-${GANGSTR_VERSION}.tar.gz \
    && rm GangSTR-${GANGSTR_VERSION}.tar.gz

WORKDIR /GangSTR-${GANGSTR_VERSION}

RUN python3 -m pip install --no-cache-dir --upgrade cmake pip

RUN mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && cmake --install . \
    && cp GangSTR /usr/local/bin/ \
    && cp thirdparty/htslib/lib/libhts.so.3 /usr/local/lib/ \
    && rm -rf /GangSTR-${GANGSTR_VERSION}

ENV SAMTOOLS_VERSION="1.20"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN python3 -m pip install --no-cache-dir git+https://github.com/broadinstitute/str-analysis

ENTRYPOINT ["/bin/bash"]
