#!/bin/bash
#
# Run tests and code quality checks.

pytest .
pylint src/ tests/
