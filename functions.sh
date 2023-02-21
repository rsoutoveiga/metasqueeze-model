#!/bin/bash

echo "the available functions are:"
echo "- parallel"
echo "- run_debug"
echo "- run_release"


parallel() {
    Rscript tools/parallel/meta-parallel.R $1 $2 $3 $4
}

run_debug() {
    echo "run_debug()"

    echo "three arguments (no argument runs default sim):"
    echo "argument 1: folder sim file name"
    echo "argument 2: sim file name"
    echo "argument 3: output file name (optional)"

    ./build/debug/metapop $1 $2 $3
}

run_release() {
    echo "run_release()"

    echo "three arguments (no argument runs default sim):"
    echo "argument 1: folder sim file name"
    echo "argument 2: sim file name"
    echo "argument 3: output file name (optional)"

    ./build/release/metapop $1 $2 $3
}

remove_debug() {
    rm -rf build/debug/*
    rm -rf build/debug/.* || true
}

remove_release() {
    rm -rf build/release/*
    rm -rf build/release/.* || true
}

remove_build() {
    rm -rf build/debug/*                                                        
    rm -rf build/debug/.* || true                                               
    rm -rf build/release/*                                                      
    rm -rf build/release/.* || true
} 


"$@"

#if declare -f "$1" > /dev/null
#then
#    # call arguments verbatim
#    "$@"
#else
#    # show a helpful error
#    echo "'$1' is not a known function name" >&2
#    exit 1
#fi

