#!/usr/bin/env bash

set -eux

cat data/illumina-drones-aa.fa data/flongle-drones-aa.fa \
    | clustalo \
    -i - \
    --threads 8 \
    --dealign --full \
    --out output/all_drones.faa \
    --distmat-out output/all_drones.dist \
    &> output/clustalo.log
