DISTFILES := $(shell find dist -name 'donuts*' | grep -v asc)
FULLTEST_CACHE_DIR := /tmp/donutscache

all: help

help:
	@echo "Tasks:"
	@echo "- test"
	@echo "- fulltest"
	@echo "- coverage"
	@echo "- package"
	@echo "- build-docs"
	@echo "- build"
	@echo "- tox"
	@echo "- sign"
	@echo "- upload"

test:
	python setup.py test

fulltest:
	mkdir -p $(FULLTEST_CACHE_DIR)
	python setup.py test -a "--cache-dir=$(FULLTEST_CACHE_DIR)" --remote-data

coverage:
	python setup.py test -a "--cache-dir=$(FULLTEST_CACHE_DIR) --cov donuts --cov-report html --cov-report term" --remote-data

build-docs:
	python setup.py build_docs

tox:
	ctox

package:
	@-rm -r dist 2>/dev/null
	@mkdir -p dist
	python setup.py sdist bdist_wheel

upload: package
	twine upload dist/*

sign:
	for file in $(DISTFILES); do gpg --detach-sign -a $$file; done

build:
	python setup.py build
