PYTHON ?= python3

.PHONY: install-dev test example lint

install-dev:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	PYTHONPATH=src $(PYTHON) -m pytest

example:
	PYTHONPATH=src $(PYTHON) examples/quickstart.py

lint:
	PYTHONPYCACHEPREFIX=/tmp/pycache $(PYTHON) -m compileall src tests examples scripts
