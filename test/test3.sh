source "test/assert.sh"

mkdir -p test/tmp/
bin/scelestial -root 5 < test/data/sample.txt > test/tmp/sample-out.txt 2>/dev/null
assert_eq "$?" "1" "Error in test"

