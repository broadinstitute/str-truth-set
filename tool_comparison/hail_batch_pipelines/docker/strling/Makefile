VERSION = 0.5.2
TAG = weisburd/strling

all: build push

build:
	docker build -t ${TAG}:${VERSION} -t ${TAG}:latest .

push:
	docker push ${TAG}
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/strling@sha256:[^\"]*/"strling@`cat sha256.txt`"/"  ../../*.py && rm ../../*.bak