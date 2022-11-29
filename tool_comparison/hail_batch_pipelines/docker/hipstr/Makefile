VERSION = 0.6.2
TAG = weisburd/hipstr

all: build push

build:
	docker build -t ${TAG}:${VERSION} -t ${TAG}:latest .

push:
	docker push ${TAG}
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
