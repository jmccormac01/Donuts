DISTFILES := $(shell find dist -name 'donuts*' | grep -v asc)

all: help

help:
	@echo "Tasks:"
	@echo "- test"
	@echo "- fulltest"
	@echo "- package"
	@echo "- build-docs"
	@echo "- build"
	@echo "- tox"
	@echo "- sign"

test:
	python setup.py test

fulltest:
	python setup.py test -a '--cache-dir=/tmp/donutscache' --remote-data

build-docs:
	python setup.py build_docs

tox:
	ctox

package:
	@rm -r dist 2>/dev/null
	@mkdir -p dist
	python setup.py sdist bdist_wheel

sign:
	for file in $(DISTFILES); do gpg --detach-sign -a $$file; done

build:
	python setup.py build
