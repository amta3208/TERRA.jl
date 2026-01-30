#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
matlab -nodisplay -r "addpath(genpath('$script_dir')); savepath; quit;"
