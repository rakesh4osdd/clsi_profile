#!/bin/bash

planemo tool_init --id 'clsi_profile' --name 'CLSI Profile' \
          --description 'MIC profile using CLSI MIC breakpoints' \
          --requirement 'python' \
          --example_command "clsi_profile_type2.py 'input.csv' 'clsi.csv' 'output.csv'" \
          --example_input test-data/input.csv \
	  --example_input test-data/clsi.csv \
          --example_output test-data/output.csv \
          --test_case \
#          --version_command 'python3 -V' \
#          --help_from_command '' \
#          --doi ''

planemo lint
