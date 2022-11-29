FROM python:3.9-slim-bullseye

RUN apt-get update && apt-get install -y curl
ENV STRLING_VERSION="v0.5.1"
RUN curl -L https://github.com/quinlan-lab/STRling/releases/download/${STRLING_VERSION}/strling -o /usr/local/bin/strling \
    && chmod 777 /usr/local/bin/strling

RUN python3 -m pip install str_analysis

ENTRYPOINT ["/bin/bash"]