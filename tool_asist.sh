#!/bin/bash

planemo tool_init --id 'asist' --name 'ASIST' --version 1.0.1 \
          --description 'Antimicrobial Susceptibility standards based phenotypes' \
          --requirement 'pandas' \
          --example_command "\"${__tool_directory__}/asist_dynamic.py\" 'asist_input.csv'  'asist_output.csv'" \
          --example_input test-data/asist_input.csv \
          --example_output test-data/asist_output.csv \
          --test_case \
	  --version_command 'python -c "import pandas; print(pandas.__version__)"' \
          --doi 'https://ab-openlab.csir.res.in/asist'

planemo lint asist.xml
