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

ENV SAMTOOLS_VERSION="1.20"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# install gcloud
ENV GCLOUD_SDK_VERISON="481.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN python3 -m pip install --no-cache-dir git+https://github.com/broadinstitute/str-analysis

ENTRYPOINT ["/bin/bash"]
