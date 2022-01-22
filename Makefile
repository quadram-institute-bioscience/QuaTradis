
VERSION := $(shell cat VERSION)
DOCKER_USERNAME := quadraminstitute
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


update_develop_branch:
	git fetch origin
	git checkout develop
	git pull origin develop

bump_major_version: update_develop_branch
	python3 bump_version.py --mode=major

bump_minor_version: update_develop_branch
	python3 bump_version.py --mode=minor

bump_version: update_develop_branch
	python3 bump_version.py --mode=patch

release:
	git commit VERSION -m "chore(package): Bump version up to $(shell cat VERSION)"
	git push origin develop
	git checkout master
	git pull origin master
	git merge --no-ff --no-edit develop
	git tag "$(shell cat VERSION)"
	git push origin master
	git push origin master --tags


release_major: bump_major_version release
	echo "Released new major version"


release_minor: bump_minor_version release
	echo "Released new minor version"


release_patch: bump_version release
	echo "Released new patch"


