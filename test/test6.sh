
source "test/assert.sh"

mkdir -p test/tmp/
bin/scelestial < test/data/sample1.txt 2>/dev/null
assert_eq "$?" "1" "Error in test"
