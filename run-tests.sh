#!/bin/bash
#
# Run tests and code quality checks.

set -e

cd $(dirname $0)

source venv/bin/activate
# mypy src/  # Call me when mypy works for real numbers
pytest
pylint src/ tests/
