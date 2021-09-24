#!/bin/bash
#
# Run tests and code quality checks.

set -e

cd $(dirname $0)

source venv/bin/activate
mypy --strict src
pytest .
pylint src/ tests/
