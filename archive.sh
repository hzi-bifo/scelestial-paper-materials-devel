mkdir -p archive
rm -rf archive/$1.tgz
tar czf archive/$1.tgz src/scelestial.h src/scelestial.cc test/data/ test/*.sh makefile
