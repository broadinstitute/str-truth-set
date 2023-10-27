FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    calc \
    parallel \
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
	libncurses5-dev \
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

RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
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
	&& rm chr*markerInfo.gz \
	&& chmod 777 ./install.sh ./setUp.sh \
	&& ./install.sh

RUN apt update && apt-get install --no-install-recommends -y gdb less


#COPY popSTRMakefile /popSTR-2.0/Makefile
COPY msGenotyper.cpp /popSTR-2.0/msGenotyper.cpp
COPY msGenotyperDefault.cpp /popSTR-2.0/msGenotyperDefault.cpp
COPY run.sh /popSTR-2.0/run.sh
COPY runPerChrom.sh /popSTR-2.0/runPerChrom.sh

RUN cd /popSTR-2.0/ \
 	&& ./install.sh \
	&& tar xzf kernel.tar.gz


ENTRYPOINT ["/bin/bash"]
