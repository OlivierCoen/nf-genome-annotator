#!/bin/bash

database=""
nbthreads=""
args=""
prefix=""

usage() {
    echo "Usage: $0 --database <db_name> --threads <nb threads> --prefix <prefix>"
    exit 1
}


while [[ "$#" -gt 0 ]]; do
    case $1 in
        --database) database="$2"; shift ;;
        --threads) nbthreads="$2"; shift ;;
        --args) args="$2"; shift ;;
        --prefix) prefix="$2"; shift ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$database" || -z "$nbthreads" || -z "$prefix" ]]; then
    echo "Error: --database and --nbthreads and --prefix are required."
    usage
fi

exit_code=0

RepeatModeler \
    -database $database \
    -threads $nbthreads \
    $args \
    > repeatmodeler.log 2>&1 || exit_code=$?

no_family_msg="No families identified"

if grep -q "No families identified.  Perhaps the database is too small" repeatmodeler.log ; then
    echo $no_family_msg
    exit 0
elif grep -q "refined-cons.fa) does not exist" repeatmodeler.log ; then
    echo $no_family_msg
    exit 0
elif [[ $exit_code -ne 0 ]]; then
    echo "Unhandled error. RepeatModeler failed with exit code $exit_code"
    exit $exit_code
fi

exit 0
