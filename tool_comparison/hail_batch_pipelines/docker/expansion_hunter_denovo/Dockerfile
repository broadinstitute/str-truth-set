FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y wget time

# install ExpansionHunterDenovo binary from https://github.com/Illumina/ExpansionHunterDenovo
ENV EH_DENOVO_VERSION="v0.9.0"
RUN wget -q https://github.com/Illumina/ExpansionHunterDenovo/releases/download/${EH_DENOVO_VERSION}/ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && tar xzf ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && mv ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64/bin/ExpansionHunterDenovo /usr/local/bin/ExpansionHunterDenovo \
    && rm -r ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64 ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz

ENTRYPOINT ["/bin/bash"]