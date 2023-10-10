TAG = weisburd/advntr

all: build push

build:
	docker build -t ${TAG}:latest .

push:
	docker push ${TAG}:latest
	docker push ${TAG} | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
