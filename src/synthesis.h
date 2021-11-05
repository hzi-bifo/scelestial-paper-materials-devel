#include <iostream>
#include <cstdlib>
#include <fstream>
#include <random>

using namespace std;

namespace synth {


std::default_random_engine generator;
std::uniform_int_distribution<int> uint_dist(0, RAND_MAX);

double doubleRand(double mx) {
   return double(uint_dist(generator)) / (double(RAND_MAX) + 1.0) * mx;
}

int intRand(int mx) {
	return (uint_dist(generator) * 1.0) / (RAND_MAX) * mx;
}


const int MAXN = 10000;
const int MAXLEN = 2000;
const int MAX_SAMPLE = 1000;

size_t n;
vector<double> adv, advCumSum;
// parent of a node, distance of node to its parent (not use for it now)
vector<int> parent, parentDistance;
vector<vector<int>> sequence;
// number of children
vector<size_t> childCount;

size_t locusCount, step;
double incAdvProb, decAdvProb, keepAdvProb;
double advIncStep, advDecStep;
double missingValueRate, zeroToOneRate, oneToZeroRate;

double stepMutationRate;

// when growing the clone tree, only allow nodes with less than treeGrowFilterOutChilderSize to have new child
int treeGrowFilterOutChilderSize = 0;

enum SamplingMethod {
    ADVANTAGE, UNIFORM
};
std::istream& operator>>(std::istream& in, SamplingMethod& samplingMethod) {
	string token;
	in >> token;
	if (token == "advantange")
		samplingMethod = SamplingMethod::ADVANTAGE;
	else if (token == "uniform")
		samplingMethod = SamplingMethod::UNIFORM;
	else
		throw runtime_error("invalid sampling method " + token);
	return in;
}

SamplingMethod samplingMethod;


int sampleWithAdvantage() {
	if (treeGrowFilterOutChilderSize == 0) {
		double prolR = doubleRand(advCumSum[n-1]);
		vector<double>::iterator prolIt = lower_bound(advCumSum.begin(), advCumSum.end(), prolR);
		int prolIdx = prolIt - advCumSum.begin();
		return prolIdx;
	} else {
		double advSum = 0;
		for (size_t i=0; i<parent.size(); i++) {
			if ((int)childCount[i] < treeGrowFilterOutChilderSize)
				advSum += adv[i];
		}
		double prolR = doubleRand(advSum), advSum2 = 0;
		for (size_t i=0; i<parent.size(); i++) {
			if ((int)childCount[i] < treeGrowFilterOutChilderSize) {
				advSum2 += adv[i];
				if (advSum2 > prolR)
					return i;
			}
		}
		return -1;
	}
}

void simulate() {
	assert(step < MAXN);
	//everything is not on vectors, so no more need to check the sizes here
	//assert(locusCount < MAXLEN);
	//init:
	sequence.clear();
	sequence.push_back(vector<int>(locusCount, 0));
	adv.clear();
	adv.push_back(1);
	advCumSum.clear();
	advCumSum.push_back(adv[0]);
	parent.clear();
	parent.push_back(-1);
	parentDistance.clear();
	parentDistance.push_back(0);
	childCount.clear();
	childCount.push_back(0);
	n = 1;

	std::poisson_distribution<int> distribution(stepMutationRate);

	
	for (size_t s = 0; s < step; s++) {
		size_t prolIdx = sampleWithAdvantage();

		assert(prolIdx < n);
		assert(0 <= prolIdx);

		sequence.push_back(sequence[prolIdx]);
		size_t cnt = distribution(generator);
		for (size_t j=0; j<cnt; j++) {
			int l = intRand(locusCount);
			sequence[n][l] = 1-sequence[n][l];
		}

		double advR = doubleRand(1);
		if (advR < incAdvProb) 
			adv.push_back(adv[prolIdx] * (1 + advIncStep));
		else if (advR < decAdvProb)
			adv.push_back(adv[prolIdx] * (1 - advDecStep));
		else
			adv.push_back(adv[prolIdx]);
		advCumSum.push_back(advCumSum[n-1] + adv[n]);

		parent.push_back(prolIdx);
		parentDistance.push_back(1);
		childCount.push_back(0);
		childCount[prolIdx] ++;
		n++;
	}

}


struct Output {
	size_t sampleCount;
	vector<vector<int>> cloneSamples;
	vector<int> sampleClone;
	vector<vector<int>> output;
	vector<int> parentCompressed,
		parentCompressedDistance;


	Output(size_t _sampleCount) : sampleCount(_sampleCount) {
		cloneSamples = vector<vector<int>>(n, vector<int>());
		parentCompressedDistance = vector<int>(n, -1);
		parentCompressed = vector<int>(n,-1);
	}

	void compressParent(int v) {
		if (parentCompressedDistance[v] != -1)
			return;
		if (parent[v] == -1) {
			parentCompressed[v] = -1;
			parentCompressedDistance[v] = 0;
			return;
		}
		compressParent(parent[v]);
		if (parent[parent[v]] != -1 && cloneSamples[parent[v]].size() == 0) {
			parentCompressed[v] = parentCompressed[parent[v]];
			parentCompressedDistance[v] = parentCompressedDistance[parent[v]] + 1;
		} else {
			parentCompressed[v] = parent[v];
			parentCompressedDistance[v] = 1;
		}
	}

};


Output sample(size_t sampleCount) {
	assert(sampleCount < MAX_SAMPLE);

	Output o(sampleCount);

	if (samplingMethod == SamplingMethod::ADVANTAGE) {
		for (size_t i=0; i<sampleCount; i++) {
			int prolIdx = sampleWithAdvantage();

			o.cloneSamples[prolIdx].push_back(i);
			o.sampleClone.push_back(prolIdx);
		}
	} else if (samplingMethod == SamplingMethod::UNIFORM) {
		if (sampleCount != childCount.size()) {
			cerr << "Warning: Sample count and child count should be equal for this sampling method " << sampleCount << " " << childCount.size() << endl;
		}
		for (size_t i=0; i<sampleCount; i++) {
			int prolIdx = i % childCount.size();
			o.cloneSamples[prolIdx].push_back(i);
			o.sampleClone.push_back(prolIdx);
		}
	} else {
		throw runtime_error("Invalid samplingMethod");
	}


	std::uniform_real_distribution<double> distribution(0.0,1.0);

	o.output = vector<vector<int>>(locusCount, vector<int>(sampleCount, -1));

	for (size_t i=0; i<sampleCount; i++) {
		for (size_t j=0; j<locusCount; j++) {
			int v = 3;
			if (distribution(generator) < missingValueRate) {
				v = 3;
			} else {
				int lv = sequence[o.sampleClone[i]][j];
				if (lv == 0) {
					if (distribution(generator) < zeroToOneRate)
						lv = 1;
				} else {
					if (distribution(generator) < oneToZeroRate)
						lv = 0;
				}
				v = lv;
			}
			o.output[j][i] = v;
		}
	}

	for (size_t i=0; i<n; i++) {
		o.compressParent(i);
	}
	return o;
}

void printSample(const Output& o, string cloneFileName, string seqFileName, string treeFileName, string trueSeqFileName) {

	ofstream fclone(cloneFileName);
	for (size_t i=0; i<n; i++) {
		fclone << i << " ";
		for (auto j: o.cloneSamples[i]) {
			fclone << j+1 << " ";
		}
		fclone << endl;
	}

	ofstream fseq(seqFileName);
	for (size_t j=0; j<locusCount; j++) {
		for (size_t i=0; i<o.sampleCount; i++) {
			fseq << o.output[j][i] << " ";
		}
		fseq << endl;
	}

	ofstream ftree(treeFileName);
	for (size_t i=0; i<n; i++) {
		if (parent[i] == -1 || o.cloneSamples[i].size() > 0)
			ftree << i << " ";
	}

	ftree << endl;
	for (size_t i=0; i<n; i++) {
		if (parent[i] != -1 && o.cloneSamples[i].size() > 0)
			ftree << i << "->" << o.parentCompressed[i] << " " << o.parentCompressedDistance[i] << endl;
	}

	if (trueSeqFileName != "") {
		ofstream fseq(trueSeqFileName);
		for (size_t j=0; j<locusCount; j++) {
			for (size_t i=0; i<o.sampleCount; i++) {
				fseq << sequence[o.sampleClone[i]][j] << " ";
			}
			fseq << endl;
		}
	}

}


}


