#!/bin/sh -e

flexbar --reads reads.fasta --target result_right --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-error-rate 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.fasta result_right.fasta`

if ! $a ; then
echo "Error testing right mode fasta"
echo $a
exit 1
else
echo "Test 1 OK"
fi


flexbar --reads reads.fasta --target result_left --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-error-rate 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left.fasta result_left.fasta`

if ! $a ; then
echo "Error testing left mode fasta"
echo $a
exit 1
else
echo "Test 2 OK"
fi


flexbar --reads reads.fasta --target result_any --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-error-rate 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any.fasta result_any.fasta`

if ! $a ; then
echo "Error testing any mode fasta"
echo $a
exit 1
else
echo "Test 3 OK"
fi


flexbar --reads reads.fasta --target result_left_tail --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-error-rate 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail.fasta result_left_tail.fasta`

if ! $a ; then
echo "Error testing left_tail mode fasta"
echo $a
exit 1
else
echo "Test 4 OK"
fi


flexbar --reads reads.fasta --target result_right_tail --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-error-rate 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail.fasta result_right_tail.fasta`

if ! $a ; then
echo "Error testing right_tail mode fasta"
echo $a
exit 1
else
echo "Test 5 OK"
fi

echo ""

