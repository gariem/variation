#!/bin/bash

BIN_DIR=/home/egarcia/workspace/github/variation/bin
TEST_DIR=/home/egarcia/workspace/github/variation/tests/test_merge

#echo "GASM:===>>"
#cat $TEST_DIR/gasm_del.vcf
#echo "PILEUP:===>>"
#cat $TEST_DIR/pileup_del.vcf
#
#echo "DIFF:===>"
#python $BIN_DIR/merge.py --gasm_file=$TEST_DIR/gasm_del.vcf \
#        --pileup_file=$TEST_DIR/pileup_del.vcf --out=./final_test.vcf
#
#cat ./final_test.vcf

#echo "GASM:===>>"
#cat $TEST_DIR/gasm_ins.vcf
#echo "PILEUP:===>>"
#cat $TEST_DIR/pileup_ins.vcf
#
#echo "DIFF:===>"
#python $BIN_DIR/merge.py --gasm_file=$TEST_DIR/gasm_ins.vcf \
#        --pileup_file=$TEST_DIR/pileup_ins.vcf

python $BIN_DIR/merge.py --gasm_file=/home/egarcia/data/igv_workfiles/vcf/DBA_2J-gasm.vcf \
        --pileup_file=/home/egarcia/data/igv_workfiles/vcf/DBA_2J.snps-indels.pileup-cram.vcf --out=./final_test.vcf --log-level=INFO