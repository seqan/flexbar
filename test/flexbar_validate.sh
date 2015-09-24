#!/bin/sh -e

echo ""

echo "Testing fasta:"
./flexbar_test_fasta.sh

echo "Testing fastq:"
./flexbar_test_fastq.sh

echo "Testing decompression:"
./flexbar_test_zip.sh

