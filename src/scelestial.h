#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <set>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <chrono>

using std::string;
using std::ostream;
using std::vector;
using std::set;
using std::map;
using std::tuple;
using std::ostringstream;
using std::to_string;
using std::istream;
using std::ostream;
using std::cin;
using std::get;
using std::istringstream;
using std::stoi;
using std::stod;
using std::min;

class ExitException: public std::exception {
	public:
	ExitException(string _reason = "") : reason(_reason) {}
	virtual const char* what() const throw() {
		return reason.c_str();
	}
	private:
	string reason;
};

template<typename T, typename C>
ostream& operator<<(ostream& os, const set<T, C>& c) {
	for (auto const &i : c)
		os << i << " ";
	return os;
}

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& c) {
	for (auto const &i : c)
		os << i << " ";
	return os;
}

template<typename T, typename K>
ostream& operator<<(ostream& os, const map<T, K>& c) {
	for (auto &i : c)
		os << i.first << ":" << i.second << " ";
	return os;
}

template<typename T>
string array2string(const T* a, const T* b) {
	ostringstream os;
	for (const T* i=a; i<b; i++)
		os << *i << " ";
	return os.str();
}

enum debug_option
{
    DEBUG_DISABLE,
    DEBUG_ENABLE,
    DEBUG_ENABLE_V
};
#ifdef DEBUG_LEVEL_1
	debug_option logLevel = DEBUG_ENABLE;
	#define IFDEBUG(A) {if(logLevel > 0) { logger << A << endl; }}
#elif  DEBUG_LEVEL_2
	debug_option logLevel = DEBUG_ENABLE_V;
	#define IFDEBUG(A) {if(logLevel > 0) { logger << A << endl; }}
#else
	debug_option logLevel = DEBUG_DISABLE;
	#define IFDEBUG(A) {}
#endif
using std::endl;

const int MAXTREELEAFS = 10;
const int MAXNODE = 3 * MAXTREELEAFS;
const double EPSILON = 1e-7;

char nucleotideAcids[] = {'A', 'T', 'C', 'G'};
const int nucleotideAcidCount = 4;
const int allelCodingSize = 11;
string allelCoding[][2] = {
	{"A","A/A"},		
						{"T","T/T"},
						{"C","C/C"},
						{"G","G/G"},
						{"K","A/C"},
						{"L","A/G"},
						{"M","C/T"},
						{"N","C/G"},
						{"O","T/G"},
						{"P","T/A"},
						{"X","./."}};
map<string, char> allelCodingMap;
map<char,int> allelCoding2Int;
char int2AllelCoding[allelCodingSize];
int xCode;
void init() {
	for (int i=0; i<allelCodingSize; i++) {
		allelCodingMap[allelCoding[i][1]] = allelCoding[i][0][0];
		string r = string("") + allelCoding[i][1][2] + allelCoding[i][1][1] + allelCoding[i][1][0];
		allelCodingMap[r] = allelCoding[i][0][0];
	}

	for (int i=0; i<allelCodingSize; i++) {
		allelCoding2Int[allelCoding[i][0][0]] = i;
		int2AllelCoding[i] = allelCoding[i][0][0];
	}

	xCode = allelCoding2Int['X'];
}

char allelCode(string allel) {
	if (allelCodingMap.find(allel) == allelCodingMap.end()) {
		if (logLevel > 0)
			logger << ("Invalid Allel " + allel) << endl;
		throw ExitException("Invalid Allel " + allel);
	}
	//logger << "C " << allel << "->" << allelCodingMap[allel];
	return allelCodingMap[allel];
}

long long power(int a, int b, int mod) {
	if (b == 0)
		return 1;
	long long p = power(a, b/2, mod);
	p = (p * p) % mod;
	if (b % 2 == 1)
		p = p * a;
	return p % mod;
}

const double IMPUTATION_COST = 0.50001;
//const double IMPUTATION_X_X_FACTOR = 2 * IMPUTATION_COST;
const double IMPUTATION_X_X_FACTOR = IMPUTATION_COST;
const int MAX_SEQUENCE = 1000;
int kRestrictionSteinerTreeMax = 4, kRestrictionSteinerTreeMin = 3;
// We know that MST/2 < Stiener-tree, thus by setting this to 0.5, this is only a valid bound,
// If we let this value to be more than 0.5, we may miss some subsets, but the algorithm would finish faster. 
double steinerTreeMSTLowerBoundRate = 0.5;

inline double charDistance(char a, char b) {
	if (a == 'X' && b == 'X')
		return IMPUTATION_X_X_FACTOR * IMPUTATION_COST;
	else
		return (a=='X' || b=='X') ? IMPUTATION_COST : a == b ? 0 : 1;
}

struct Cell {
	vector<char> s;

	Cell() {
	}

	int size() const {
		return (int) s.size();
	}

	void load(string str) {
		s.clear();
		for (auto& c: str)
			append(c);
	}

	void append(char a) {
		s.push_back(a);
	}

	void set(char a, int pos) {
		s[pos] = a;
	}

	char get(int pos) const {
		return s[pos];
	}

	double distance(const Cell& c) const {
		double cnt = 0;
		for (int i=0; i<(int)s.size(); i++) {
			cnt += charDistance(s[i], c.s[i]);
		}
		return cnt;
	}

	string toString() const {
		string r(s.begin(), s.end());
		return r;
	}

};


ostream& operator<<(ostream& os, const Cell& s) {
	return os << s.toString();
}



struct DisjointSetArray {
	mutable vector<int> parent;

	DisjointSetArray() : parent(MAX_SEQUENCE) {
		for(int i=0; i<MAX_SEQUENCE; i++)
			parent[i] = -1;
	}

	template<typename It>
	DisjointSetArray(It begin, It end) : parent(MAX_SEQUENCE)  {
		for (It i = begin; i != end; i++) {
			parent[*i] = *i;
		}
	}

	void add(int t) {
		parent[t] = t;
	}

	int parentStar(int t) const {
		if (parent[t] == t)
			return t;
		return parent[t] = parentStar(parent[t]);
	}

	void join(int a, int b) {
		parent[parentStar(a)] = parent[parentStar(b)];
	}

	bool isJoint(int a, int b) const {
		return parentStar(a) == parentStar(b);
	}
};

ostream& operator<<(ostream& os, const DisjointSetArray& ds) {
	os << "{";
	for (int i=0; i<MAX_SEQUENCE; i++) {
		if (ds.parent[i] != -1) {
			os << i << ":" << ds.parent[i] << " ";
		}
	}
	return os << "}";
}


struct GenerateAllTrees;

ostream& operator<<(ostream& os, const GenerateAllTrees& a);

struct GenerateAllTrees {
	//parameters
	int n;
	bool unRooted;

	//inner variables
	vector<int> nodeIndex;
	int nodeCount;
	vector<set<int>> availableLeafs;
	/** children of node v in tree is in childs[v]. If childs[v].size() == 1, v is a leaf and childs[v][0] is the index of the leaf, 
	 * on which this tree is built. 
	 */ 
	vector<vector<int>> childs;

	//output
	vector<string> trees;
	char treeRepresentation[100 * 8 * 3]; // MAXN * 3 * MAXLGN
	int treeRepresentationEnd;
	vector< vector<vector<int>> > outChilds;

	//debug
	int depth;

	GenerateAllTrees(int _n, bool _unRooted = true) : n(_n), unRooted(_unRooted) {
	}

	vector<vector<vector<int>>> run() {
		trees.clear();
		depth = 0;
		treeRepresentationEnd = 0;
		availableLeafs.push_back(set<int>());
		for (int i=0; i<n; i++) {
			availableLeafs.back().insert(i);
		}
		treeRepresentationAppend("(");
		childs.clear();
		childs.push_back(vector<int>());
		nodeIndex.push_back(0);
		nodeCount = 1;
		outChilds.clear();
		rec();
		return outChilds;
	}

	int treeRepresentationAppend(string s) {
		int l = treeRepresentationEnd;
		for (int i=0, j=treeRepresentationEnd; i< (int)s.length(); i++, j++) {
			treeRepresentation[j] = s[i];
		}
		treeRepresentationEnd += (int) s.length();
		return l;
	}

	void treeRepresentationMoveBack(int lt) {
		treeRepresentation[lt] = 0;
		treeRepresentationEnd = lt;
	}

	void forAllSubsets(const set<int>& refSet, set<int>::const_iterator i, set<int>& s1, set<int>& s2) {
		//logger << " 2^S: " << *i << " [" << s1 << "] [" << s2 << "]" << endl;
		if (i == refSet.end()) {
			if (s2.size() > 0 || refSet.size() == 1) {
				//Set is ready.
				availableLeafs.push_back(s2); // replaces the current available leafs for the current node with the remaining nodes, s2. 
				//Node that, current set, refSet, is already removed from the list.
				availableLeafs.push_back(s1); // leafs for the newly generated node.

				childs[nodeIndex.back()].push_back(nodeCount);
				childs.push_back(vector<int>());
				nodeIndex.push_back(nodeCount);
				nodeCount++;
				int lt = treeRepresentationAppend("(");
				rec();
				treeRepresentationMoveBack(lt);
				nodeCount--;
				nodeIndex.pop_back();
				childs.pop_back();
				childs[nodeIndex.back()].pop_back();

				availableLeafs.pop_back();
				availableLeafs.pop_back();
			}
			return;
		}
		set<int>::const_iterator j = i;
		j++;
		s1.insert(*i);
		forAllSubsets(refSet, j, s1, s2);
		s1.erase(s1.find(*i));
		if (i != refSet.begin()) {
			//always put first to the first set
			s2.insert(*i);
			forAllSubsets(refSet, j, s1, s2);
			s2.erase(s2.find(*i));
		}
	}

	void rec() {
		depth++;
		//logger << "Rec: \n" << *this << endl;

		if (availableLeafs.size() == 0) {
			//Tree is ready
			if (!unRooted || childs[0].size() >= 3) {
				trees.push_back(treeRepresentation);
				outChilds.push_back(childs);
			}
		} else if (availableLeafs.back().size() == 0) {
			//we return to a node, all its available childs are already handled.
			int p = nodeIndex.back();
			nodeIndex.pop_back();
			availableLeafs.pop_back();
			int lt = treeRepresentationAppend(")");
			rec();
			treeRepresentationMoveBack(lt);
			availableLeafs.push_back(set<int>());
			nodeIndex.push_back(p);
		} else if (availableLeafs.back().size() == 1 && childs[nodeIndex.back()].size() == 0) {
			//this is a new child.
			int n = *availableLeafs.back().begin();
			availableLeafs.pop_back();

			int p = nodeIndex.back();
			nodeIndex.pop_back();
			childs[p].push_back(n); // no need to be removed after rec, it will be removed from forAllSubsets imediately after
			int lt = treeRepresentationAppend(to_string(n) + ")");
			rec();
			treeRepresentationMoveBack(lt);
			nodeIndex.push_back(p);

			set<int> r;
			r.insert(n);
			availableLeafs.push_back(r);
		} else {
			set<int> r = availableLeafs.back();
			availableLeafs.pop_back();

			set<int> s1, s2;
			forAllSubsets(r, r.begin(), s1, s2);

			availableLeafs.push_back(r);
		}

		depth--;
	}
};

ostream& operator<<(ostream& os, const GenerateAllTrees& a) {
	char tab[100];
	for (int i=0; i<a.depth * 2; i++)
		tab[i] = ' ';
	tab[a.depth * 2] = 0;

	int j=0;
	for (auto i: a.availableLeafs) {
		os << tab << j << ": " << i << endl;
		j++;
	}

	os << tab << "  T:   " << a.treeRepresentation << endl;
	return os;
}

struct UniverseVertexSet {

	vector<Cell> cells;

	int length() const {
		if (cells.size() == 0)
			return 0;
		return cells[0].size();
	}

	int size() const {
		return (int) cells.size();
	}

	int add(const Cell& c) {
		cells.push_back(c);
		return (int)(cells.size()) - 1;
	}

	Cell& getVertex(int index) {
		return cells[index];
	}

	const Cell& getVertex(int index) const {
		return cells[index];
	}

	double distance(int v, int u) const {
		return cells[v].distance(cells[u]);
	}

	int precalculatedPairwiseDistanceSize;
	double precalculatedPairwiseDistance[MAX_SEQUENCE][MAX_SEQUENCE];

	void precalculatePairwiseDistances() {
		assert(cells.size() < MAX_SEQUENCE);
		precalculatedPairwiseDistanceSize = (int) cells.size();
		for (int i=0; i<(int)cells.size(); i++) {
			precalculatedPairwiseDistance[i][i] = 0;
			for (int j=i+1; j<(int)cells.size(); j++)
				precalculatedPairwiseDistance[i][j] = precalculatedPairwiseDistance[j][i] = distance(i, j);
		}
	}

	double inputSequencesDistance(int v, int u) const {
		if (v >= precalculatedPairwiseDistanceSize || u >= precalculatedPairwiseDistanceSize) {
			if (logLevel > 0)
				logger << "!! inputSeqDistance requested but not pre calculated " << v << " " << u << " " << precalculatedPairwiseDistanceSize << endl;
			return distance(v, u);
		}
		return precalculatedPairwiseDistance[v][u];
	}

};

struct EdgeWeight;
ostream& operator<<(ostream& os, const EdgeWeight& e);

struct EdgeWeight {
	int v, u;
	double w;

	EdgeWeight(int _v=0, int _u=0, double _w=0) : v(_v), u(_u), w(_w) {}

	bool operator<(const EdgeWeight& b) const {
		const EdgeWeight& a = *this;
		long long w1 = (long long)(a.w/EPSILON),
			w2 = (long long)(b.w/EPSILON);
		if (w1 != w2) {
			return w1 < w2;
		}
		if (a.v != b.v) {
			return a.v < b.v;
		}
		return a.u < b.u;
	}
};

ostream& operator<<(ostream& os, const EdgeWeight& e) {
	return os << "(" << e.v << " " << e.u << " " << e.w << ")";
}


/**
 */ 
tuple<vector<EdgeWeight>, double> imputeTree(
		const vector<int>& leaves, 
		const vector<vector<int>>& tree, 
		UniverseVertexSet& universeVertexSet, 
		bool onlyCost) {
	IFDEBUG("imputeTree called ...")
	int len = universeVertexSet.length();

	double finalCost = 0;
	vector<vector<int>> imputedTree(tree.size(), vector<int>(len, -1));

	IFDEBUG("imputeTree started ...")

	for (int l=0; l<len; l++) {
#ifdef SMALL_STACK_CONFIG
		vector<vector<double>> cost(MAXNODE, vector<double>(allelCodingSize, 0));
		vector<vector<vector<int>>> path(MAXNODE, vector<vector<int>>(allelCodingSize, vector<int>(MAXNODE, 0)));
#else
		double cost[MAXNODE][allelCodingSize];
		int path[MAXNODE][allelCodingSize][MAXNODE];
#endif
		for (int v=(int)(tree.size())-1; v>=0; v--) {
			for (int a = 0; a<allelCodingSize; a++) {
				if (tree[v].size() == 1) {
					// a leaf!
					int c = tree[v][0];
					const Cell& leaf = universeVertexSet.getVertex(leaves[c]);

					// HADI APPROX: 
					cost[v][a] = charDistance(leaf.s[l], int2AllelCoding[a]);
					// (leaf.sIndex[l] == a) ? 0 : ((leaf.s[l] == 'X') ? IMPUTATION_COST : 100000);
				} else {
					// HADI APPROX: 
					double cstSum = 0;
					for (auto c: tree[v]) {
						#ifndef DP_DOUBLE_TO_INT_APPROX
						// HADI APPROX: 
						double cst = cost[c][a];
						#else
						long long cst = cost[c][a];
						#endif
						path[v][a][c] = a;
						for (int aa=0; aa<allelCodingSize; aa++)
							if (aa != xCode && cst > cost[c][aa] + 1) {
								cst = cost[c][aa]+1;
								path[v][a][c] = aa;
							}
						cstSum += cst;
					}
					cost[v][a] = cstSum;
				}
			}
		}

		int minA = -1;
		for (int a=0; a<allelCodingSize; a++) {
			if (a != xCode && (minA == -1 || cost[0][a] < cost[0][minA]))
				minA = a;
		}

		// HADI APPROX: 
		finalCost += cost[0][minA];

		if (!onlyCost) {
			vector<int> minAs(tree.size(), -1);
			minAs[0] = minA;
			for (int v=0; v<(int)tree.size(); v++) {
				imputedTree[v][l] = minAs[v];
				if (tree[v].size() == 1) {
				} else {
					for (auto c: tree[v]) {
						minAs[c] = path[v][minAs[v]][c];
					}
				}
			}
		}

	}
	// return imputedTree;
	
	IFDEBUG("imputed tree cost calculated")


	if (!onlyCost) {
		vector<EdgeWeight> gEdges;
		vector<int> treeNode2UniverseIndex(tree.size(), -1);
		for (int v=(int)(imputedTree.size())-1; v>=0; v--) {
			Cell c;
			for (auto loc : imputedTree[v]) {
				c.append(int2AllelCoding[loc]);
			}
			int h = universeVertexSet.add(c);
			treeNode2UniverseIndex[v] = h;

			if (tree[v].size() == 1) {
				int leafIndex = leaves[tree[v][0]];
				gEdges.push_back(EdgeWeight(
						leafIndex, h, 
						universeVertexSet.distance(leafIndex, h)));
			} else {
				for (auto u : tree[v]) {
					int uIndex = treeNode2UniverseIndex[u];
					gEdges.push_back(EdgeWeight(
						h, uIndex,
						universeVertexSet.distance(h, uIndex)
					));
				}
			}
		}
		return make_tuple(gEdges, finalCost);
	}
	return make_tuple(vector<EdgeWeight>(), finalCost);
}

tuple<vector<EdgeWeight>, double> imputeTree(
		const vector<int>& leaves, 
		const vector<vector<int>>& tree, 
		UniverseVertexSet& universeVertexSet) {
	return imputeTree(leaves, tree, universeVertexSet, false);
}

double imputeTreeCost(
		const vector<int>& leaves, 
		const vector<vector<int>>& tree, 
		UniverseVertexSet& universeVertexSet) {
	IFDEBUG("imputeTreeCost called")
	//return get<1>(imputeTree(leaves, tree, universeVertexSet, true));
	tuple<vector<EdgeWeight>, double> r = imputeTree(leaves, tree, universeVertexSet, true);
	return std::get<1>(r);
}


template<typename T>
struct SubsetIterator {
	vector<int> l; 
	int n, k;

	const vector<T>& referenceSet;

	SubsetIterator(int _n, int _k, const vector<T>& _referenceSet) : n(_n), k(_k), referenceSet(_referenceSet) {
		init();
	}

	void init() {
		l.clear();
		for (int i=0; i<k; i++)
			l.push_back(i);
	}

	bool isValid() const {
		return l.size() > 0;
	}

	void next() {
		for (int nn=n; l.size() > 0 && l.back() == nn-1; ) {
			l.pop_back();
			nn--;
		}
		if (l.size() == 0)
			return;
		l.back()++;
		while ((int)l.size() < k)
			l.push_back(l.back()+1);
		return;
	}

	vector<T> get() {
		vector<T> li;
		for (auto i: l){
			li.push_back(referenceSet[i]);
		}
		return li;
	}

	static int size(int n, int k1, int k2) {
		int sum = 0;
		for (int k=k1; k<=k2; k++) {
			int soorat = 1, makhraj = 1;
			for (int i=0; i<k; i++) {
				soorat *= (n-i);
				makhraj *= (i+1);
			}
			sum += soorat / makhraj;
		}
		return sum;
	}
};

double treeCost(const vector<EdgeWeight>& edges, UniverseVertexSet& universeVertexSet) {
	double cst = 0;
	for (auto e: edges) {
		cst += e.w;
	}
	return cst;
}

template<typename T, typename C>
bool includesInSet(const set<T>& A, const set<T, C>& B) {
	for (auto i: A) {
		if (B.find(i) == B.end())
			return false;
	}
	return true;
}

template<typename T, typename C>
bool includesInSet(const vector<T>& A, const set<T, C>& B) {
	for (auto i: A) {
		if (B.find(i) == B.end())
			return false;
	}
	return true;
}

template<typename T, typename C>
void subtractEq(set<T, C>& A, const vector<T>& B) {
	for (auto i: B) {
		if (A.find(i) != A.end())
			A.erase(A.find(i));
	}
}


template<typename T>
void add_to_container(vector<T>& container, const T& element) {
	container.push_back(element);
}

template<typename T, typename C>
void add_to_container(set<T, C>& container, const T& element) {
	container.insert(element);
}

template<typename T>
void sort_container(vector<T>& container) {
	sort(container.begin(), container.end());
}


template<typename T, typename C>
void sort_container(set<T, C>& container) {
}


// V == vector<EdgeWeigt> or set<EdgeWeight, ???>
template<typename V>
double mstEq(V& M) {
	DisjointSetArray ds;
	set<int> vertices;
	for (auto e : M) {
		ds.add(e.v);
		ds.add(e.u);
		vertices.insert(e.v);
		vertices.insert(e.u);
	}
	int n = (int) vertices.size();
	sort_container(M);

	double cost = 0;
	V R;
	int addedEdge = 0;
	for (auto e: M) {
		if (addedEdge >= n - 1)
			break;
		if (!ds.isJoint(e.v, e.u)) {
			ds.join(e.v, e.u);
			add_to_container(R, e);
			cost += e.w;
		}
	}
	M = R;
	return cost;
}

tuple<vector<EdgeWeight>, double> bermenApplyCandidateTrees(UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, const vector<EdgeWeight>& T, 
			const vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>>& candidateStack,
			bool applyOnUniverse);

double bermenCandidateTreeCost(UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, const vector<EdgeWeight>& T, 
			const vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>>& candidateStack) {
	auto t = bermenApplyCandidateTrees(universeVertexSet,
		input, n, T, 
		candidateStack,
		false);
	return get<1>(t);
}

vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>> bermenGenerateCandidateTrees(
	UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, int minkk, int kk, vector<EdgeWeight>& T
) {

	vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>> candidateStack;

	if (logLevel > 0)
		logger << "Init Tree: " << T << " for k=" << kk << " cost=" << bermenCandidateTreeCost(
					universeVertexSet,
						input, n, T, 
						candidateStack
				) << endl; 


	int stepsPassed = 0, stepsLBDd = 0;
	int stepsTotal = SubsetIterator<int>::size(n, minkk, kk);
	auto startTime = std::chrono::steady_clock::now();

	for (int k = minkk; k<=kk; k++) {

		if (logLevel > 0) {
			logger << "Generating all trees " << k << " ..." << endl;
		}
		GenerateAllTrees alg(k);
		if (logLevel > 0) {
			logger << "Generating all trees " << k << " created" << endl;
		}
		vector<vector<vector<int>>> trees = alg.run();
		if (logLevel > 0) {
			logger << "Generating all trees executed" << endl;
		}

		//fill candidateStack
		SubsetIterator<int> choose(n, k, input);
		IFDEBUG("choose created")


		for (int tt=0; choose.isValid(); tt++, choose.next()) {


			if (logLevel > 0) {
				stepsPassed++;
				if (stepsPassed % 1000 == 0) {
				//if (stepsPassed % 1 == 0) {
					auto end = std::chrono::steady_clock::now();

					//Following is too time consuming, 
					// but is good to have a quote about the performance of the algorithm during steps
					double cost = bermenCandidateTreeCost(
						universeVertexSet,
							input, n, T, 
							candidateStack
					);

					logger << "step " << stepsPassed << "/" << stepsTotal << " lbd=" << stepsLBDd << " time=" << std::chrono::duration_cast<std::chrono::seconds>(end - startTime).count() << " k=" << k << " cost=" << cost << endl; 
					// return candidateStack; //BAD BAD
				}
			}

			vector<int> li = choose.get();

			sort(T.begin(), T.end());
			IFDEBUG("Edges sorted " << T.size())
			/**
			 * terminalTreeDiscardingEdges = bridges
			 */ 

			vector<EdgeWeight> terminalTreeRemainingEdges, terminalTreeDiscardingEdges;
			double bridgeCost = 0; // the cost of separating l from each other in T
			{
				DisjointSetArray ds;
				for (auto n : input) {
					ds.add(n);
				}
				for (int i=0; i+1<(int)li.size(); i++) {
					ds.join(li[i], li[i+1]);
				}
				IFDEBUG("ds joins done")

				for (auto &e: T) {
					if (!ds.isJoint(e.v, e.u)) {
						ds.join(e.v, e.u);
						terminalTreeRemainingEdges.push_back(e);
					} else {
						terminalTreeDiscardingEdges.push_back(e);
						bridgeCost += e.w;
					}
				}
				IFDEBUG("edges joined")
			}

			IFDEBUG("Choose done")

			double minC = 1000000000;
			IFDEBUG("calculating minC " << minC << "minT ...")
			vector<vector<int>> minT;
			IFDEBUG("calculating minC " << minC << "minT (2) ...")

			for (vector<vector<int>> &tree : trees) {
			//for (auto &tree : trees) {
				IFDEBUG("impute tree " << tree << " ... ")
				double c = imputeTreeCost(li, tree, universeVertexSet);
				IFDEBUG("impute tree " << tree << " done ")
				if (c < minC) {
					minC = c;
					minT = tree;
				}
			}
			IFDEBUG("minC " << minC << "minT found")

			{

				double gain = bridgeCost - minC;
				if (gain > EPSILON) {
					DisjointSetArray ds, ds2;
					for (auto n : input) {
						ds.add(n);
						ds2.add(n);
					}
					for (int i=0; i+1<(int)li.size(); i++) {
						ds.join(li[i], li[i+1]);
					}

					for (auto &e: T) {
						if (!ds.isJoint(e.v, e.u)) {
							ds.join(e.v, e.u);
							ds2.join(e.v, e.u);
						} else {
						}
					}

					//returns terminal of a component
					map<int, int> dsRoot2TerminalMap;
					for (auto t : li) {
						dsRoot2TerminalMap[ds2.parentStar(t)] = t;
					}

					vector<EdgeWeight> terminalNewEdges;
					for (auto e: terminalTreeDiscardingEdges) {
						if (ds2.isJoint(e.v, e.u)) {
							if (logLevel > 0)
								logger << "Discarding edges are connected " << e << endl;
							ostringstream os;
							os << "Discarding edges are connected " << e;
							throw ExitException(os.str());
						}

						int t1 = dsRoot2TerminalMap[ds2.parentStar(e.v)],
							t2 = dsRoot2TerminalMap[ds2.parentStar(e.u)];

						//We can keep the track of the tree here too!
						EdgeWeight ne = EdgeWeight(t1, t2, e.w - gain);
						terminalNewEdges.push_back(ne);
					}

					// A good thing to see on logger
					
					candidateStack.push_back(make_tuple(minT, terminalTreeDiscardingEdges, terminalNewEdges, li));
					

					T = terminalTreeRemainingEdges;
					T.insert(T.end(), terminalNewEdges.begin(), terminalNewEdges.end());

					if (logLevel > 0) {
						double cost = bermenCandidateTreeCost(
							universeVertexSet,
								input, n, T, 
								candidateStack
						);

						logger << " candid: " << li << " gain=" << gain;
						logger << "  s=" << stepsPassed << " cost=" << cost << endl; 
					}

					assert((int)T.size() == n-1);


				} 
			}
			IFDEBUG("gain calculated")

		}
			
	}


	if (logLevel > 0) {
		auto end = std::chrono::steady_clock::now();
		logger << "step " << stepsPassed << "/" << stepsTotal << " lbd=" << stepsLBDd << " time=" << std::chrono::duration_cast<std::chrono::seconds>(end - startTime).count() << endl; 
	}
	if (logLevel > 0) {
		logger << "Candidate trees calculated" << endl;
	}

	return candidateStack;
}

tuple<vector<EdgeWeight>, double> bermenApplyCandidateTrees(UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, const vector<EdgeWeight>& T, 
			const vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>>& candidateStack,
			bool applyOnUniverse) {
	
	vector<EdgeWeight> gEdges;
	double cost = 0;

	//check candidateStack
	//and apply candidateStack
	{
		DisjointSetArray ds;

		set<EdgeWeight> E = set<EdgeWeight>(T.begin(), T.end()),
			M = E;

		for (auto n : input) {
			ds.add(n);
		}


		for (int i=(int)(candidateStack.size())-1; i >=0; i--) {
			auto t = candidateStack[i];

			vector<vector<int>> tree = get<0>(t);
			vector<EdgeWeight> terminalTreeDiscardingEdges = get<1>(t), 
				terminalNewEdges = get<2>(t);
			vector<int> leavesIndices = get<3>(t);

			assert(includesInSet(M, E));
			
			// E += terminalTreeDiscardingEdges;
			E.insert(terminalTreeDiscardingEdges.begin(), terminalTreeDiscardingEdges.end());
			if (includesInSet(terminalNewEdges, M)) {
				//apply
				if (applyOnUniverse) {
					tuple<vector<EdgeWeight>, double> t = imputeTree(
						leavesIndices,
						tree,
						universeVertexSet,
						false
					);
					vector<EdgeWeight> ee = get<0>(t);
					gEdges.insert(gEdges.end(), ee.begin(), ee.end());
					
					cost += get<1>(t);
				} else {
					double t = imputeTreeCost(
						leavesIndices,
						tree,
						universeVertexSet
					);
					cost += t;
				}
				for (auto e: terminalNewEdges) {
					ds.join(e.v, e.u);
				}

			} else {
				subtractEq(E, terminalNewEdges);
				M = E; 
				mstEq(M);
			}

			assert((int)M.size() == n-1);
		}


		for (auto e: M) {
			if (!ds.isJoint(e.v, e.u)) {
				double ew = universeVertexSet.distance(e.v, e.u);
				gEdges.push_back(EdgeWeight(
					e.v, e.u,
					ew
				));
				cost += ew;


				ds.join(e.v, e.u);
			}
		}


	}

	return make_tuple(gEdges, cost);
}

/**
 * Implementation of: mproved Approximations for the Steiner Tree Problem
 * Authors: Piotr Berman and Viswanathan Ramaiyer
 */ 
tuple<vector<EdgeWeight>, double> optimizeTree(
		UniverseVertexSet& universeVertexSet,
		const vector<int>& input, int minkk, int kk) {
	/**
	 * T always contain a tree between terminal nodes.
	 */

	assert(kk < MAXTREELEAFS);

	double distanceMatrix[MAX_SEQUENCE][MAX_SEQUENCE];

	vector<EdgeWeight> T;
	for (auto v: input) {
		for (auto u: input) {
			double w = universeVertexSet.distance(v, u);
			T.push_back(EdgeWeight(v, u, w));
			distanceMatrix[v][u] = distanceMatrix[u][v] = w;
		}
	}
	mstEq(T);

	int n = (int)input.size();

	//Tree, distartedEdges, newEdges, treeLeaveIndices(tau)
	vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>> candidateStack = 
		bermenGenerateCandidateTrees(universeVertexSet,
			input, n, minkk, kk, T);

	tuple<vector<EdgeWeight>, double> t = bermenApplyCandidateTrees(
		universeVertexSet,
			input, n, T, 
			candidateStack,
			true
	);
	if (logLevel > 0) {
		logger << "Trees applied" << endl; 
	}


	return t;
}

// Removes degree two edges which does not provide any other information regarding input vertices
tuple<vector<EdgeWeight>, vector<int>> compressGraph(UniverseVertexSet& universeVertexSet, const vector<EdgeWeight>& edges, const vector<int>& inputCells) {
	DisjointSetArray ds;
	for (int i=0; i<universeVertexSet.size(); i++) {
		ds.add(i);
	}

	for (auto e: edges) {
		if (e.w == 0) {
			ds.join(e.v, e.u);
		}
	}

	map<int, int> minOfSet;
	for (int i=0; i<universeVertexSet.size(); i++) {
		int p = ds.parentStar(i);
		if (minOfSet.find(p) == minOfSet.end())
			minOfSet[p] = p;
		minOfSet[p] = min(minOfSet[p], i);
	}

	vector<EdgeWeight> E;
	for (auto e: edges) {
		if (e.w != 0) {
			E.push_back(EdgeWeight(
				minOfSet[ds.parentStar(e.v)],
				minOfSet[ds.parentStar(e.u)],
				e.w
			));
		}
	}

	set<int> V;
	for (auto m: minOfSet) {
	   V.insert(m.second); 
	}

	//We should add input vertices, even there are more than one with equal sequences.
	for (auto v: inputCells) {
		if (V.find(v) == V.end()) {
			V.insert(v);
			E.push_back(EdgeWeight(
				v,
				minOfSet[ds.parentStar(v)],
				0
			));
		}
	}

	return make_tuple(E, vector<int>(V.begin(), V.end()));
}

void printResultAsGraph(ostream& os, UniverseVertexSet& universeVertexSet, const vector<EdgeWeight>& edges, double cost, const vector<int>& cells, const map<int, Cell>& imputation) {

	tuple<vector<EdgeWeight>, vector<int>> tt = compressGraph(universeVertexSet, edges, cells);

	if (logLevel > 0)
		logger << "Tree compressed" << endl;

	vector<EdgeWeight> e = get<0>(tt);
	vector<int> v = get<1>(tt);

	set<int> inputCells;
	for (auto c: cells) {
		inputCells.insert(c);
	}
	
	os << v.size() << endl;
	for (int j=0; j<(int)v.size(); j++) {
		int i = v[j];
		os << i << " " << 
			(inputCells.find(i) == inputCells.end() ? 0 : 1) << " " << 
			universeVertexSet.getVertex(i).toString() << " " << 
			(imputation.find(i) != imputation.end() ? imputation.find(i)->second.toString() : "-") << endl;
	}
	os << e.size() << endl;
	for (auto ee: e) {
		os << ee.v << " " << ee.u << " " << ee.w << endl;
	}
}

void rewriteCellXValuesTo(Cell& to, const Cell& from) {
	for (int i=0; i<to.size(); i++) {
		if (to.s[i] == 'X' && from.s[i] != 'X')
			to.s[i] = from.s[i];
	}
}


map<int,Cell> calculateImputation(UniverseVertexSet& universeVertexSet, const vector<EdgeWeight>& edges, const vector<int>& inputCells) {
	set<int> input(inputCells.begin(), inputCells.end());
	map<int, Cell> imputation;
	for (auto v: inputCells) {
		imputation[v] = universeVertexSet.getVertex(v);
	}
	for (auto e: edges) {
		if (input.find(e.v) != input.end()) {
			rewriteCellXValuesTo(imputation[e.v], universeVertexSet.getVertex(e.u));
		} else if (input.find(e.u) != input.end()) {
			rewriteCellXValuesTo(imputation[e.u], universeVertexSet.getVertex(e.v));
		}
	}
	Cell aCell;
	for (int i=0; i<universeVertexSet.length(); i++) {
		aCell.append('A');
	}
	for (auto v: inputCells) {
		rewriteCellXValuesTo(imputation[v], aCell);
	}
	return imputation;
}
