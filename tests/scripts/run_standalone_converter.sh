#!/usr/bin/env bash

set -eu

input_file_base=${1}

TEST_INPUT_DIR=${TEST_DIR}/inputFiles
TEST_OUTPUT_DIR=testOutputs
mkdir -p ${TEST_OUTPUT_DIR}

input_file=${TEST_INPUT_DIR}/${input_file_base}
output_file=${TEST_OUTPUT_DIR}/${input_file_base/.slcio/.edm4hep.root}
patch_file=${TEST_OUTPUT_DIR}/${input_file_base/.slcio/_colls.txt}

echo "Creating the patch file for the standalone converter"
check_missing_cols --minimal ${input_file} > ${patch_file}

echo "Running the standalone converter"
lcio2edm4hep ${input_file} ${output_file} ${patch_file}

echo "Comparing the converted and original contents"
compare-contents ${input_file} ${output_file}
