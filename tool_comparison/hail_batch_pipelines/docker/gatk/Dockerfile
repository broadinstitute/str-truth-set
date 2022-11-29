FROM broadinstitute/gatk:4.3.0.0

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
RUN apt-get update && apt-get install --no-install-recommends -y \
    wget \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    build-essential
RUN apt-get update && apt-get install --no-install-recommends -y python3.7

ENV SAMTOOLS_VERSION="1.9"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN python3.7 -m pip install str_analysis

ENTRYPOINT ["/bin/bash"]