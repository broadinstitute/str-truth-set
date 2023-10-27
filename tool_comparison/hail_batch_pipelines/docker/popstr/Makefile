VERSION = 2.0
TAG = weisburd/popstr

all: build push

build:
	docker build -t ${TAG}:${VERSION} -t ${TAG}:latest .

push:
	docker push ${TAG}
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/popstr@sha256:[^\"]*/"popstr@`cat sha256.txt`"/"  ../../*.py && rm ../../*.bak