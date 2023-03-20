FROM google/cloud-sdk:slim

# The base image current includes python3.9

RUN python3 -m pip install --upgrade \
	pandas==1.3.5 \
	numpy==1.21.5 \
	matplotlib==3.5.1 \
	seaborn==0.11.2 \
	str-analysis==0.9.19 \
	tqdm

COPY *.py .

ENTRYPOINT ["/bin/bash"]