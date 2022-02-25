
source "test/assert.sh"

mkdir -p test/tmp/
bin/scelestial -root 3 -no-internal-sample < test/data/sample.txt > test/tmp/sample-nis-r3.txt
diff test/tmp/sample-nis-r3.txt test/data/sample-nis-r3.txt
