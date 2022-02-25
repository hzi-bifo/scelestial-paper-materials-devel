source "test/assert.sh"

mkdir -p test/tmp/
bin/scelestial -root 3 < test/data/sample2.txt > test/tmp/sample2-out.txt
diff test/data/sample2-out.txt test/tmp/sample2-out.txt
