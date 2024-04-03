#!/usr/bin/env bash

# Ensure we have the latest golden input and output data from the S3 Bucket
aws s3 sync s3://aria-tools/tests/regression/ . --no-sign-request

# Run the extract test cases
./run_extract_test.py

# Run the tssetup test cases
./run_tssetup_test.py

# Run the test suite
python -m pytest
