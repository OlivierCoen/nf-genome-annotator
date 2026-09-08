#!/bin/bash

ARGS="$@"

exit_code=0

earlGrey \
    $ARGS \
    | tee earlgrey.log 2>&1 || exit_code=$?

no_family_msg="No families identified"

if grep -q "No families identified.  Perhaps the database is too small" earlgrey.log ; then
    echo $no_family_msg
    exit 0
elif grep -q "refined-cons.fa) does not exist" earlgrey.log ; then
    echo $no_family_msg
    exit 0
elif [[ $exit_code -ne 0 ]]; then
    echo "Unhandled error. earlGrey failed with exit code $exit_code"
    exit $exit_code
fi

exit 0
