#!/bin/bash

function postprocess_terra {
    local calling_folder=$(pwd)
    if [ $# == 0 ]; then # assuming this is called from inside the case folder 
        folder_name="${PWD##*/}"
        cd ..
    elif [ $# == 1 ]; then
        folder_name="$1"
    else 
        echo "Unexpected additional argument. Only one argument is supported: case_folder_name."
        exit 1
    fi
    case_name="$folder_name"
    echo "Postprocessing output for: $case_name"
    echo "Folder name is $folder_name"
    matlab -nodisplay -r "postProcessTERRA('$folder_name','$case_name', true, true); quit;"
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