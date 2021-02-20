SCELESTIALSRC=../src
SCRIPT=bin/
BOOSTINCLUDE=/usr/local/Cellar/boost/1.69.0_2/include/

ALL:$(SCRIPT)/scelestial $(SCRIPT)/synthesis 

$(SCRIPT)/scelestial: $(SCELESTIALSRC)/scelestial.cc $(SCELESTIALSRC)/scelestial.h $(SCELESTIALSRC)/util.h
	g++ -O2 -Wconversion -Wno-sign-conversion -Wno-shorten-64-to-32 -std=c++0x $(SCELESTIALSRC)/scelestial.cc -o $(SCRIPT)/scelestial -Wsign-compare

$(SCRIPT)/synthesis: $(SCELESTIALSRC)/synthesis.cc $(SCELESTIALSRC)/synthesis.h
	g++ -std=c++11 $(SCELESTIALSRC)/synthesis.cc -o $(SCRIPT)/synthesis -I$(BOOSTINCLUDE) -lboost_program_options

clean:
	rm -rf result $(ALL)
