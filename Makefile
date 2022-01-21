
VERSION := $(shell cat VERSION)
DOCKER_USERNAME := maplesond
DOCKER_IMAGENAME := quatradis
DOCKER_PATH := ${DOCKER_USERNAME}/${DOCKER_IMAGENAME}

dev:
	python3 setup.py develop

install:
	python3 setup.py install

test:
	cd tests && pytest

clean:
	rm -r build dist quatradis.egg-info

docker-build:
	docker build -t ${DOCKER_PATH}:${VERSION} .
	docker tag ${DOCKER_PATH}:${VERSION} ${DOCKER_PATH}:dev

docker-push: docker-build
	echo "${DOCKER_PASSWORD}" | docker login -u "${DOCKER_USERNAME}" --password-stdin
	docker tag ${DOCKER_PATH}:${VERSION} ${DOCKER_PATH}:latest
	docker push ${DOCKER_PATH}:latest
	docker push ${DOCKER_PATH}:${VERSION}
    
