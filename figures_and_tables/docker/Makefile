VERSION = v5
TAG = weisburd/truth-set-figures

all: build push

build:
	cd .. && docker build -t ${TAG}:${VERSION} -t ${TAG}:latest -f docker/Dockerfile .

push:
	docker push ${TAG}:${VERSION} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
