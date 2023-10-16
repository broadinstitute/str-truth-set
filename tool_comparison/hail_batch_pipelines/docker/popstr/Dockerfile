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
	build-essential \
	unzip


ENV SAMTOOLS_VERSION="1.9"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

#RUN python3 -m pip install str_analysis==0.9.17

RUN apt-get install libboost-all-dev --no-install-recommends -y \
	&& wget https://github.com/DecodeGenetics/popSTR/archive/2.0.zip \
	&& unzip 2.0.zip \
	&& rm 2.0.zip \
	&& cd popSTR-2.0/ \
	&& chmod 777 ./install.sh ./setUp.sh \
	&& ./install.sh


ENTRYPOINT ["/bin/bash"]
