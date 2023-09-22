#!/usr/bin/env bash

set -eu

TEST_INPUT_DIR=${TEST_DIR}/inputFiles
mkdir -p ${TEST_INPUT_DIR}

TEST_FILES=("ild_higgs_rec.slcio" "ild_higgs_dst.slcio")

for input_file in ${TEST_FILES[@]}; do
    if [ ! -f ${TEST_INPUT_DIR}/${input_file} ]; then
        echo "${input_file} test input file not yet present. Fetching it"
        wget https://key4hep.web.cern.ch/testFiles/k4EDM4hep2LcioConv/ILD/${input_file} -P ${TEST_INPUT_DIR}
    fi
done
