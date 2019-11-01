ALL: bin/scelestial bin/synth

bin/scelestial: bin/ src/scelestial.cc src/scelestial.h src/util.h
	g++ -O2 -Wconversion -Wno-sign-conversion -Wno-shorten-64-to-32 -std=c++0x src/scelestial.cc -o bin/scelestial -Wsign-compare

bin/synth: src/synthesis.cc src/synthesis.h
	g++ -std=c++11 src/synthesis.cc -o bin/synth -I/usr/local/Cellar/boost/1.69.0_2/include/ -lboost_program_options

bin/:
	mkdir bin

clean:
	rm -rf bin/
