mkdir -p archive
rm -rf archive/$1.tgz
tar czf archive/$1.tgz src/scelestial.h src/scelestial.cc src/synthesis.cc src/synthesis.h test/data/ test/*.sh makefile
