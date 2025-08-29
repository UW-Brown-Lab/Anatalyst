#!/bin/bash

# copy pipeline code to the VM working dir
cp -r /app/* .

# the config should either reference data on /staging
# or data should be copied over first
# cp -r /staging/... .

# run the pipeline on the copied config.yml
python ./scripts/run_pipeline.py --config ./config.yml
