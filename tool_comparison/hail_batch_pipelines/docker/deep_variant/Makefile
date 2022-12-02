VERSION=1.4.0
TAG = weisburd/deep-variant

all: build push

build:
	docker build -t ${TAG}:${VERSION} -t ${TAG}:latest .

push:
	docker push ${TAG}:${VERSION}
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
