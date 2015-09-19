
flexbar --reads test.fastq --target result_right --format fastq-i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.fastq result_right.fastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, right"
echo $a
exit -1
else
echo "Test 1 OK"
fi


flexbar --reads test.fastq --target result_left --format fastq-i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left.fastq result_left.fastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, left"
echo $a
exit -1
else
echo "Test 2 OK"
fi


flexbar --reads test.fastq --target result_any --format fastq-i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any.fastq result_any.fastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode any, left"
echo $a
exit -1
else
echo "Test 3 OK"
fi


flexbar --reads test.fastq --target result_left_tail --format fastq-i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail.fastq result_left_tail.fastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, left_tail"
echo $a
exit -1
else
echo "Test 4 OK"
fi


flexbar --reads test.fastq --target result_right_tail --format fastq-i1.5 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail.fastq result_right_tail.fastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode fasta, right_tail"
echo $a
exit -1
else
echo "Test 5 OK"
fi

echo ""

