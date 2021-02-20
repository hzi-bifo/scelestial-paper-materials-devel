ALL: bin/scelestial bin/synthesis

bin/scelestial: src/scelestial.cc src/scelestial.h src/util.h
	g++ -O2 -Wconversion -Wno-sign-conversion -Wno-shorten-64-to-32 -std=c++0x src/scelestial.cc -o bin/scelestial -Wsign-compare

bin/synthesis: src/synthesis.cc src/synthesis.h
	g++ -std=c++11 src/synthesis.cc -o bin/synthesis -I/usr/local/Cellar/boost/1.69.0_2/include/ -I/home/hforoughmand/miniconda3/envs/gcc/include/ -lboost_program_options -L/home/hforoughmand/miniconda3/envs/gcc/lib/

bin/:
	mkdir bin

clean:
	rm -rf bin/
