#!/bin/sh -e

flexbar --reads test.fastq --target result_right --format i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.fastq result_right.fastq`

if ! $a ; then
echo "Error testing right mode fastq"
echo $a
exit -1
else
echo "Test 1 OK"
fi


flexbar --reads test.fastq --target result_left --format i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left.fastq result_left.fastq`

if ! $a ; then
echo "Error testing left mode fastq"
echo $a
exit -1
else
echo "Test 2 OK"
fi


flexbar --reads test.fastq --target result_any --format i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any.fastq result_any.fastq`

if ! $a ; then
echo "Error testing any mode fastq"
echo $a
exit -1
else
echo "Test 3 OK"
fi


flexbar --reads test.fastq --target result_left_tail --format i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail.fastq result_left_tail.fastq`

if ! $a ; then
echo "Error testing left_tail mode fastq"
echo $a
exit -1
else
echo "Test 4 OK"
fi


flexbar --reads test.fastq --target result_right_tail --format i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail.fastq result_right_tail.fastq`

if ! $a ; then
echo "Error testing right_tail mode fastq"
echo $a
exit -1
else
echo "Test 5 OK"
fi

echo ""

