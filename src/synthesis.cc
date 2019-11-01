#include <boost/program_options.hpp>
#include "synthesis.h"

using namespace std;
namespace po = boost::program_options;
using namespace synth;

ostream& operator<<(ostream& os, vector<bool> v) {
	for (auto b: v)
		os << b;
	return os;
}

Output removeDuplicates(const Output& o) {
	vector<bool> dup(o.sampleCount, false);
	for (int i=0; i<o.sampleCount; i++)
		for (int j=i+1; j<o.sampleCount; j++) {
			bool eq = true;
			for (unsigned k=0; eq && k<o.output.size(); k++) {
				eq &= (o.output[k][i] == o.output[k][j]);
			}
			if (eq)
				dup[j] = true;
		}
	set<int> samples;
	vector<int> newSampleIndex(o.sampleCount, -1);
	int sc = 0;
	for (int i=0; i<o.sampleCount; i++)
		if (!dup[i]) {
			newSampleIndex[i] = sc;
			sc++;
			samples.insert(i);
		}
	Output r(sc);
	for (unsigned i=0; i<o.cloneSamples.size(); i++) {
		r.cloneSamples.push_back(vector<int>());
		for (unsigned j=0; j<o.cloneSamples[i].size(); j++)
			if (samples.find(o.cloneSamples[i][j]) != samples.end())
				r.cloneSamples[i].push_back(newSampleIndex[o.cloneSamples[i][j]]);
	}
	for (unsigned i=0; i<o.sampleCount; i++)
		if (samples.find(i) != samples.end())
			r.sampleClone.push_back(o.sampleClone[i]);
	for (unsigned k=0; k<o.output.size(); k++) {
		r.output.push_back(vector<int>());
		for (int i=0; i<o.sampleCount; i++) {
			if (samples.find(i) != samples.end()) 
				r.output[k].push_back(o.output[k][i]);
		}
	}
	for (int i=0; i<n; i++) {
		r.compressParent(i);
	}
	return r;
}

int main(int argc, char* argv[]) {

	int seed;
	int seedError;
	int sampleCount;
	string cloneFileName, seqFileName, treeFileName, trueSeqFileName;
	bool removeDuplicateSequences = false;


	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("sample", po::value<int>(&sampleCount)->required(), "number of sampled cells")
		("step", po::value<int>(&step)->required(), "number of simulation steps")
		("locus", po::value<int>(&locusCount)->required(), "number of locuses")
		// ("locus-potential", po::value<int>(&potentialLocusCount)->default_value(-1), "number of potential locuses in genome")
		("aic", po::value<double>(&incAdvProb)->required(), "advantage increase rate")
		("adc", po::value<double>(&decAdvProb)->required(), "advantage decrease rate")
		("akc", po::value<double>(&keepAdvProb)->required(), "advantage keep rate")
		("seed", po::value<int>(&seed)->default_value(7), "random seed number")
		("seed-error", po::value<int>(&seedError)->default_value(7), "random seed for sequencing error")
		("step-mutation-rate", po::value<double>(&stepMutationRate)->default_value(100), "average number of mutations in each step")
		("ais", po::value<double>(&advIncStep)->default_value(1), "advantage increase step")
		("ads", po::value<double>(&advDecStep)->default_value(1), "advantage decrease step")
		("fclone", po::value<string>(&cloneFileName)->required(), "clones file name")
		("fseq", po::value<string>(&seqFileName)->required(), "seq file name")
		("fseqtrue", po::value<string>(&trueSeqFileName)->default_value(""), "seq (without error) filename")
		("ftree", po::value<string>(&treeFileName)->required(), "tree file name")
		("mvrate", po::value<double>(&missingValueRate)->default_value(0.2), "missing value rate")
		("zorate", po::value<double>(&zeroToOneRate)->default_value(0.1), "zero to one rate")
		("ozrate", po::value<double>(&oneToZeroRate)->default_value(0.2), "one to zero rate")
		("remove-dup", po::value<bool>(&removeDuplicateSequences)->default_value(false), "remove duplicate sequences")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}
	po::notify(vm);


	seedError = seed;
	srand(seed);
	generator.seed(seedError);

	double probSum = incAdvProb + decAdvProb + keepAdvProb;
	incAdvProb /= probSum;
	decAdvProb /= probSum;
	keepAdvProb /= probSum;

	simulate();
	Output o = sample(sampleCount);
	if (removeDuplicateSequences) {
		o = removeDuplicates(o);
	}
	printSample(o, cloneFileName, seqFileName, treeFileName, trueSeqFileName);

	return 0;
}
