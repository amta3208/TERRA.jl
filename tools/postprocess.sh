#!/bin/bash

function postprocess_terra {
    local calling_folder
    calling_folder="$(pwd)"

    local script_dir matlab_tools_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    matlab_tools_dir="$script_dir/matlab"

    if [[ ! -d "$matlab_tools_dir" ]]; then
        echo "Expected MATLAB tooling directory not found: $matlab_tools_dir" >&2
        return 1
    fi

    local case_path case_name
    if [[ $# -eq 0 ]]; then
        case_path="$(pwd)"
        case_name="$(basename "$case_path")"
    elif [[ $# -eq 1 ]]; then
        if [[ ! -d "$1" ]]; then
            echo "Case directory not found: $1" >&2
            return 1
        fi
        case_path="$(cd "$1" && pwd)"
        case_name="$(basename "$case_path")"
    else
        echo "Unexpected additional argument. Only one argument is supported: case_directory." >&2
        return 1
    fi

    echo "Postprocessing output for: $case_name"
    echo "Case path is $case_path"

    # Escape single quotes for MATLAB single-quoted strings.
    local matlab_case_path matlab_case_name matlab_matlab_tools_dir
    matlab_case_path=${case_path//\'/\'\'}
    matlab_case_name=${case_name//\'/\'\'}
    matlab_matlab_tools_dir=${matlab_tools_dir//\'/\'\'}

    MATLABPATH="$matlab_tools_dir" matlab -nodisplay -r "addpath(genpath('$matlab_matlab_tools_dir')); postProcessTERRA('$matlab_case_path', '$matlab_case_name', true, true); quit;"

    echo "Finished postprocessing output for: $case_name"
    cd "$calling_folder"
}

function main {
    postprocess_terra
}

# Execute the function postprocess_terra if the bash script is called directly.
# Otherwise, don't run the function if this script is just being sourced.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
