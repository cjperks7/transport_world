#!/bin/bash

input_filename=$1
output_filename=$2

/usr/local/cmod/codes/efit/bin/fast_efitdd < "$input_filename" > "$output_filename"
status=$?

if [ $status -ne 0 ]; then
    echo "Error: Segfault occurred"
    # Handle the segfault (e.g., log the error, clean up, retry, etc.)
    exit 1
fi

exit 0
