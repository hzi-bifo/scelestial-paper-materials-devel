mkdir -p test/tmp/
bin/scelestial -root 4 < test/data/sample.txt > test/tmp/sample-reroot4.txt
diff test/tmp/sample-reroot4.txt test/data/sample-reroot4.txt
