VERSION = 1.4.1
TAG = weisburd/straglr

all: build push

build:
	docker build -t ${TAG}:${VERSION} -t ${TAG}:latest .

push:
	docker push ${TAG}:${VERSION}
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/straglr@sha256:[^\"]*/"straglr@`cat sha256.txt`"/"  ../../*.py && rm ../../*.bak