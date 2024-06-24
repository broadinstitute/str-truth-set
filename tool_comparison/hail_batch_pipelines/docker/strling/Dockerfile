FROM python:3.9-slim-bullseye

RUN apt-get update && apt-get install -y curl time

ENV STRLING_VERSION="v0.5.2"
RUN curl -L https://github.com/quinlan-lab/STRling/releases/download/${STRLING_VERSION}/strling -o /usr/local/bin/strling \
    && chmod 777 /usr/local/bin/strling

# install bedtools
ENV BED_TOOLS_VERSION="2.31.0"
RUN curl https://github.com/arq5x/bedtools2/releases/download/v${BED_TOOLS_VERSION}/bedtools.static -o bedtools.static \
	&& chmod 777 bedtools.static \
	&& mv bedtools.static /usr/bin/bedtools

RUN python3 -m pip install --no-cache-dir --upgrade pip \
 	&& python3 -m pip install --no-cache-dir str_analysis

ENTRYPOINT ["/bin/bash"]