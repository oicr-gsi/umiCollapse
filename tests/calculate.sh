#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find output files, return their md5sums to std out
find . -name "*umiCountsPerPosition.tsv" -xtype f -exec sh -c "cat {} | md5sum " \;
find . -name "*results.json" -xtype f -exec sh -c "cat {} | md5sum |sort" \;
find . -name "*.bam" -xtype f -exec sh -c "cat {} | md5sum " \;
ls | sed 's/.*\.//' | sort | uniq -c
