#define logger std::cerr
#include "scelestial.h"
using std::cerr;

#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            throw ExitException(message); \
        } \
    } while (false)
            //std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
            //          << " line " << __LINE__ << ": " << message << std::endl; \


int newRoot = -1;
bool noInternalSample = false;

void load(UniverseVertexSet& universeVertexSet, int argc, char* argv[]) {
	bool addRoot = false;
	Cell addRootCell;

	for (int i=1; i<argc; i++) {
		//std::cerr << "arg: " << argv[i] << endl;
		if (string(argv[i]) == "-help") {
			std::cout << "usage: scelestial [options] <input >output" << std::endl
				  << "  -min-k arg             Sets min-k to arg. Default=3. " << std::endl
				  << "                         Scelestial considers all k-subsets of samples for min-k <= k <= max-k." << std::endl
				  << "  -max-k arg             Sets max-k to arg. Default=4. " << std::endl
				  << "                         Scelestial considers all k-subsets of samples for min-k <= k <= max-k." << std::endl
				  << "  -include-root arg      Adds a sample with all sites equal to arg as root of the tree." << std::endl
				  << "  -root root-index       Zero-based index of the new root." << std::endl
				  << "                         After inferring the tree, tree edges are directed toward new root." << std::endl
				  << "                         In the output line \"u v w\" u is the parent and v is the child node." << std::endl

				  << "  -no-internal-sample    Move all samples to leaf nodes." << std::endl
				  << "                         If -root is also present, a neighbor of root-index is chosen as the root." << std::endl;
			exit(0);
		} else if (string(argv[i]) == "-include-root") {
			if (i+1 < argc) {
				string allel = argv[i+1];
				Cell c;
				for (int j = 0; j < universeVertexSet.cells[0].size(); j++) {
					c.append(allelCode(allel));
				}
				addRoot = true;
				addRootCell = c;
				i++;
			} else {
				throw ExitException("No argument for root");
			}
		} else if (string(argv[i]) == "-root") {
			if (i+1 < argc) {
				newRoot = stoi(argv[i+1]);
				i++;
			} else {
				throw ExitException("Invalid root index argument");
			}
		} else if (string(argv[i]) == "-max-k") {
			if (i+1 < argc) {
				int kk = stoi(argv[i+1]);
				kRestrictionSteinerTreeMax = kk;
				i++;
			} else {
				throw ExitException("No argument for k");
			}
		} else if (string(argv[i]) == "-min-k") {
			if (i+1 < argc) {
				int kk = stoi(argv[i+1]);
				kRestrictionSteinerTreeMin = kk;
				i++;
			} else {
				throw ExitException("No argument for k");
			}
		} else if (string(argv[i]) == "-no-internal-sample") {
			noInternalSample = true;
		} else if (string(argv[i]) == "-v") {
			logLevel = debug_option::DEBUG_ENABLE;
		} else if (string(argv[i]) == "-vv") {
			logLevel = debug_option::DEBUG_ENABLE_V;
		} else {
			throw ExitException("Invalid argument " + string(argv[i]));
		}
	}

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

	if (addRoot) {
		universeVertexSet.add(addRootCell);
	}

	ASSERT(universeVertexSet.size() >= kRestrictionSteinerTreeMax, "Sample size should be not smaller than max-k");

	assert(universeVertexSet.size() < MAX_SEQUENCE);
	assert(kRestrictionSteinerTreeMax < MAXTREELEAFS);
	//assert(kRestrictionSteinerTreeMin <= kRestrictionSteinerTreeMax, "Min-k should be not more than Max-k");
}


int main(int argc, char* argv[]) {
	try {
		init();
		UniverseVertexSet universeVertexSet;
		load(universeVertexSet, argc, argv);
		if (logLevel > 0)
			logger << "Loaded" << endl;

		if (newRoot >= universeVertexSet.size() || newRoot < -1) {
			throw ExitException("Selected sample for re-rooting is invalid: " + std::to_string(newRoot));
		}

		vector<int> cells;
		for (int i=0; i<universeVertexSet.size(); i++)
			cells.push_back(i);
		tuple<vector<EdgeWeight>, double> t = optimizeTree(universeVertexSet, cells, kRestrictionSteinerTreeMin, kRestrictionSteinerTreeMax );

		if (logLevel > 0)
			logger << "Tree optimized" << " cost=" << get<1>(t) << endl;

		map<int,Cell> imputation = calculateImputation(universeVertexSet, get<0>(t), cells);

		//following will change the graph!
		printResultAsGraph(std::cout, universeVertexSet, get<0>(t), get<1>(t), cells, imputation, newRoot, noInternalSample);

		if (logLevel > 0)
			logger << "Done" << endl;
	} catch(const ExitException& e)  {  
		cerr << "Error: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
