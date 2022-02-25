mkdir -p test/tmp/
bin/scelestial < test/data/sample.txt > test/tmp/sample-out.txt
diff test/tmp/sample-out.txt test/data/sample-out.txt
