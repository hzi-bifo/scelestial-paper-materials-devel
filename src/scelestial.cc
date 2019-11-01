#define logger std::cerr
#include "scelestial.h"
using std::cerr;

void load(UniverseVertexSet& universeVertexSet, int argc, char* argv[]) {
	istream& f = cin;
	universeVertexSet.cells.clear();
	for (string l; getline(f, l); ) {
		//logger << "Line loaded " << l << endl;
		istringstream ss(l);
		string s;

		for (int i=0; ss >> s; i++) {
			if (i > 0) {
				if ((int)universeVertexSet.cells.size() <= i-1) {
					universeVertexSet.cells.push_back(Cell());
				}
				universeVertexSet.cells[i-1].append(allelCode(s));
			}
		}
	}

	for (int i=0; i<argc; i++) {
		if (string(argv[i]) == "-include-root") {
			if (i+1 < argc) {
				string allel = argv[i+1];
				Cell c;
				for (int j = 0; j < universeVertexSet.cells[0].size(); j++) {
					c.append(allelCode(allel));
				}
				universeVertexSet.add(c);
			} else {
				throw "No argument for root";
			}
		}
		if (string(argv[i]) == "-max-k") {
			if (i+1 < argc) {
				int kk = stoi(argv[i+1]);
				kRestrictionSteinerTreeMax = kk;
			} else {
				throw "No argument for k";
			}
		} 
		if (string(argv[i]) == "-min-k") {
			if (i+1 < argc) {
				int kk = stoi(argv[i+1]);
				kRestrictionSteinerTreeMin = kk;
			} else {
				throw "No argument for k";
			}
		} 
		if (string(argv[i]) == "-v") {
			logLevel = debug_option::DEBUG_ENABLE;
		}
		if (string(argv[i]) == "-vv") {
			logLevel = debug_option::DEBUG_ENABLE_V;
		}
	}

	assert(universeVertexSet.size() < MAX_SEQUENCE);
	assert(kRestrictionSteinerTreeMax < MAXTREELEAFS);
	// assert(kRestrictionSteinerTreeMin <= kRestrictionSteinerTreeMax);
}


int main(int argc, char* argv[]) {
	try {
		init();
		UniverseVertexSet universeVertexSet;
		load(universeVertexSet, argc, argv);
		if (logLevel > 0)
			logger << "Loaded" << endl;

		vector<int> cells;
		for (int i=0; i<universeVertexSet.size(); i++)
			cells.push_back(i);
		tuple<vector<EdgeWeight>, double> t = optimizeTree(universeVertexSet, cells, kRestrictionSteinerTreeMin, kRestrictionSteinerTreeMax );

		if (logLevel > 0)
			logger << "Tree optimized" << " cost=" << get<1>(t) << endl;

		map<int,Cell> imputation = calculateImputation(universeVertexSet, get<0>(t), cells);

		//following will change the graph!
		printResultAsGraph(std::cout, universeVertexSet, get<0>(t), get<1>(t), cells, imputation);

		if (logLevel > 0)
			logger << "Done" << endl;
	} catch(const ExitException& e)  {  
		cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

	return 0;
}
