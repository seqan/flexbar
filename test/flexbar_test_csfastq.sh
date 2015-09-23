#!/bin/sh -e

flexbar --reads test.csfastq --target result_right --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.csfastq result_right.csfastq`

if ! $a ; then
echo "Error testing right mode csfastq"
echo $a
exit 1
else
echo "Test 1 OK"
fi


flexbar --reads test.csfastq --target result_left --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left.csfastq result_left.csfastq`

if ! $a ; then
echo "Error testing left mode csfastq"
echo $a
exit 1
else
echo "Test 2 OK"
fi


flexbar --reads test.csfastq --target result_any --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any.csfastq result_any.csfastq`

if ! $a ; then
echo "Error testing any mode csfastq"
echo $a
exit 1
else
echo "Test 3 OK"
fi


flexbar --reads test.csfastq --target result_left_tail --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail.csfastq result_left_tail.csfastq`

if ! $a ; then
echo "Error testing left_tail mode csfastq"
echo $a
exit 1
else
echo "Test 4 OK"
fi


flexbar --reads test.csfastq --target result_right_tail --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail.csfastq result_right_tail.csfastq`

if ! $a ; then
echo "Error testing right_tail mode csfastq"
echo $a
exit 1
else
echo "Test 5 OK"
fi

echo ""

