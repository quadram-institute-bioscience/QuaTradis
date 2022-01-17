

dev:
	python3 setup.py develop

install:
	python3 setup.py install

test:
	cd tests && pytest

clean:
	rm -r build dist quatradis.egg-info

docker-build:
	docker build -t quadraminstitute/quatradis:dev .

#docker-push:
#    docker push
