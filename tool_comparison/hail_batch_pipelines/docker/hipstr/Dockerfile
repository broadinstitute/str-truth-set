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


RUN python3 -m pip install str_analysis==0.9.17

ENTRYPOINT ["/bin/bash"]
