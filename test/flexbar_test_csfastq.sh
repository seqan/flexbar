
flexbar --reads test.csfastq --target result_right --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.csfastq result_right.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfastq, right"
echo $a
exit -1
else
echo "Test 1 OK"
fi


flexbar --reads test.csfastq --target result_left --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_left.csfastq result_left.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfastq, left"
echo $a
exit -1
else
echo "Test 2 OK"
fi


flexbar --reads test.csfastq --target result_any --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end ANY > /dev/null

a=`diff correct_result_any.csfastq result_any.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfastq, any"
echo $a
exit -1
else
echo "Test 3 OK"
fi


flexbar --reads test.csfastq --target result_left_tail --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end LEFT_TAIL > /dev/null

a=`diff correct_result_left_tail.csfastq result_left_tail.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfastq, left_tail"
echo $a
exit -1
else
echo "Test 4 OK"
fi


flexbar --reads test.csfastq --target result_right_tail --format csfastq --adapter-min-overlap 4 --adapters adapters_cs.fasta --min-read-length 10 --adapter-threshold 1 --adapter-trim-end RIGHT_TAIL > /dev/null

a=`diff correct_result_right_tail.csfastq result_right_tail.csfastq`

l1=`expr length "$a"`

if [ $l1 != 0 ]; then
echo "error testing mode csfastq, right_tail"
echo $a
exit -1
else
echo "Test 5 OK"
fi

echo ""

