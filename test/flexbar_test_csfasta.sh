#!/bin/sh -e

flexbar --reads test.csfasta --target result_right --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.csfasta result_right.csfasta`

if ! $a ; then
echo "Error testing right mode csfasta"
echo $a
exit 1
else
echo "Test 1 OK"
fi


flexbar --reads test.csfasta --target result_left --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left.csfasta result_left.csfasta`

if ! $a ; then
echo "Error testing left mode csfasta"
echo $a
exit 1
else
echo "Test 2 OK"
fi


flexbar --reads test.csfasta --target result_any --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any.csfasta result_any.csfasta`

if ! $a ; then
echo "Error testing any mode csfasta"
echo $a
exit 1
else
echo "Test 3 OK"
fi


flexbar --reads test.csfasta --target result_left_tail --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail.csfasta result_left_tail.csfasta`

if ! $a ; then
echo "Error testing left_tail mode csfasta"
echo $a
exit 1
else
echo "Test 4 OK"
fi


flexbar --reads test.csfasta --target result_right_tail --color-space --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail.csfasta result_right_tail.csfasta`

if ! $a ; then
echo "Error testing right_tail mode csfasta"
echo $a
exit 1
else
echo "Test 5 OK"
fi

echo ""

