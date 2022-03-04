
VERSION := $(shell cat VERSION)
DOCKER_USERNAME := sbastkowski
DOCKER_IMAGENAME := quatradis
DOCKER_PATH := ${DOCKER_USERNAME}/${DOCKER_IMAGENAME}

dev:
	python3 setup.py develop

install:
	python3 setup.py install

unit_test:
	cd tests && pytest --cov --cov-report=xml --doctest-modules

script_test:
	./scripts/tradis --help 2>&1 | grep tradis > /dev/null
	./scripts/tradis tags --help | grep tags > /dev/null
	./scripts/tradis tags add --help | grep tags > /dev/null
	./scripts/tradis tags check --help | grep tags > /dev/null
	./scripts/tradis tags filter --help | grep tags > /dev/null
	./scripts/tradis tags remove --help | grep tags > /dev/null
	./scripts/tradis plot --help | grep plot > /dev/null
	./scripts/tradis plot create --help | grep plot > /dev/null
	./scripts/tradis plot combine --help | grep plot > /dev/null
	./scripts/tradis plot analyse --help | grep plot > /dev/null
	./scripts/tradis utils --help | grep utils > /dev/null
	./scripts/tradis utils index --help | grep utils > /dev/null
	./scripts/tradis pipeline --help | grep pipeline > /dev/null
	./scripts/tradis pipeline single --help | grep pipeline > /dev/null
	./scripts/tradis pipeline multiple --help | grep pipeline > /dev/null
	#./scripts/tradis_comparison.R --help
	#./scripts/tradis_essentiality.R --help



test: unit_test script_test


sonar: test
	sonar-scanner

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


update_master_branch:
	git fetch origin
	git checkout master
	git pull origin master

bump_major_version: update_master_branch
	python3 bump_version.py --mode=major

bump_minor_version: update_master_branch
	python3 bump_version.py --mode=minor

bump_version: update_master_branch
	python3 bump_version.py --mode=patch

release:
	git commit VERSION -m "chore(package): Bump version up to $(shell cat VERSION)"
	git tag "$(shell cat VERSION)"
	git push origin master
	git push origin master --tags


release_major: bump_major_version release
	echo "Released new major version"

release_minor: bump_minor_version release
	echo "Released new minor version"

release_patch: bump_version release
	echo "Released new patch"


