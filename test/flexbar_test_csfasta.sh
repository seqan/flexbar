
flexbar --reads test.csfasta --target result_right --format csfasta --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.csfasta result_right.csfasta`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfasta, right"
echo $a
exit -1
else
echo "Test 1 OK"
fi


flexbar --reads test.csfasta --target result_left --format csfasta --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left.csfasta result_left.csfasta`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfasta, left"
echo $a
exit -1
else
echo "Test 2 OK"
fi


flexbar --reads test.csfasta --target result_any --format csfasta --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any.csfasta result_any.csfasta`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfasta, any"
echo $a
exit -1
else
echo "Test 3 OK"
fi


flexbar --reads test.csfasta --target result_left_tail --format csfasta --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail.csfasta result_left_tail.csfasta`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfasta, left_tail"
echo $a
exit -1
else
echo "Test 4 OK"
fi


flexbar --reads test.csfasta --target result_right_tail --format csfasta --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail.csfasta result_right_tail.csfasta`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfasta, right_tail"
echo $a
exit -1
else
echo "Test 5 OK"
fi

echo ""

