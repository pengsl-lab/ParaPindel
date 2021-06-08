#!/bin/sh

mkdir test_results

../pindel -f ../test/SmallTest/hs_ref_chr20.fa -p ../test/SmallTest/COLO-829_20-p_ok-1.txt -o test_results/colo_test.out -k -s -l -T 1

../pindel -f ../test/SmallTest/sim1chrVs2.fa -i ../test/SmallTest/sim1chrVs2.conf_for_test -o test_results/simulated_test.out -k -s -l -T 1

../pindel2vcf -R SIMCHROM -r ../test/SmallTest/sim1chrVs2.fa -P test_results/simulated_test.out -d 00000000

echo
echo "Start of diffs"
diff -q gold_standard test_results
echo "End of diffs"

