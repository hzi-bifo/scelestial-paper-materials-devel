# Scelestial: Single Cell Lineage Tree Inference based on a Steiner Tree Approximation Algorithm

Scelestial is a fast application for phylogeny reconstruction on single-cell data. The algorithm is based on an approximation algorithm for the Steiner problem and could be considered as a generalization of neighbor-joining method in which more than two samples are chosen for merge.

Table of Contents
=================
   * [Installation](#installation)
   * [Running Scelestial](#running-scelestial)
      * [Sample input &amp; output](#sample-input--output)
   * [Running Scelestial Explained](#running-scelestial-explained)
      * [Input format](#input-format)
      * [Output format](#output-format)
      * [Moving all the samples to the leaf nodes](#moving-all-the-samples-to-the-leaf-nodes)
      * [Re-rooting the tree](#Re-rooting-the-tree)
   * [Running Scelestial Explained](#running-scelestial-explained)
   * [Running Pre-defined Tests](#running-pre-defined-tests)
   * [Running Scelestial on Simulated Data](#running-scelestial-on-simulated-data)
      * [Evaluating the results](#evaluating-the-results)
         * [Comparing sample distances](#comparing-sample-distances)
         * [Comparing partition similarity](#comparing-partition-similarity)
      * [Generating PDF](#generating-pdf)
   * [Reproducing results of paper](#reproducing-results-of-paper)
      * [Preparation](#preparation)
      * [Test scripts](#test-scripts)

## Installation:
Scelestial could be installed via conda under bioconda channel. For bioconda-based installation, first install [miniconda](https://docs.conda.io/en/latest/miniconda.html) (or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)) and then 
install via following command
```bash
conda install -c conda-forge libcxx 'boost-cpp=1.74.0' 'boost=1.74.0'
conda install -c bioconda scelestial
```
Then installation could be checked by executing following command
```
scelestial -help
```

Alternatively, scelestial could be installed from source. Download source codes via
```
git clone git@github.com:hzi-bifo/scelestial-paper-materials-devel.git scelestial
```
Change folder and make
```
cd scelestial
make bin/scelestial
```
Now you can run scelestial with relative or aboslute address
```
bin/scelestial -help
```

To build `synthesis` address of header files and library folder of boost program options should be fixed in the makefile. Then following command builds `synthesis`
```
make bin/synthesis
```
And following command executes it `synthesis` 
```
bin/synthesis --help
```


## Running Scelestial
In short, run scelestial with
```
bin/scelestial <[input-file] >[output-file]
```
The detail of input/output formats and how to
generate synthetic data follows.

### Sample input & output:
Input (`data/sample.txt`):
```
1 G/G ./. ./. C/C G/G
2 C/G C/G C/G C/G ./.
3 C/G C/G ./. ./. ./.
4 G/G ./. A/A G/G G/G
5 C/C A/A ./. ./. A/A
6 A/A T/T ./. ./. ./.
7 ./. T/T ./. ./. A/A
8 ./. A/A A/A C/C ./.
9 ./. C/C ./. ./. ./.
10 T/T T/T G/G ./. ./.
11 ./. G/G G/G ./. ./.
12 C/C ./. ./. ./. ./.
13 ./. ./. A/A C/C C/C
14 ./. ./. G/G C/C C/C
15 ./. T/T ./. T/T ./.
```
The input consists of 5 samples with 15 sites.

Output:
```
6
0 1 GNNGCAXXXTXCXXX GNNGCAAACTGCCCT
1 1 XNNXATTACTGXXXT GNNGATTACTGCCCT
2 1 XNXAXXXAXGGXAGX ANAAAAAAAGGAAGA
3 1 CNXGXXXCXXXXCCT CNNGAAACCTGCCCT
4 1 GXXGAXAXXXXXCCX GNNGAAAACTGCCCT
5 0 GNNGAAAACTGCCCT -
5
4 5 4.50009
3 5 6.00008
1 5 4.50005
0 5 4.50007
1 2 6.50013
```

The inferred tree is shown in the following figure
<img src="https://github.com/hzi-bifo/scelestial-paper-materials-devel/blob/master/output-impute.png?raw=true" width="50%">


The output represents a tree with 6 nodes. The first 5 nodes are the input nodes (specified by second column, 1: input samples, 0: new nodes). Imputed sequences are shown as the 4th column. For example sample 0 is imputed as `GNNGCAAACTGCCCT`. In this sequence we use folloing coding to represent duplets `"A":"A/A", "T":"T/T", "C":"C/C", "G":"G/G", "K":"A/C", "L":"A/G", "M":"C/T", "N":"C/G", "O":"T/G", "P":"T/A", "X":"./."`.





## Running Scelestial Explained:

Scelestial is easy to be executed. It accepts input from standard input and prints output to the standard output. There are a few optional arguments for customizing behavior of Scelestial.

```
usage: scelestial [options] <input >output
  -min-k arg             Sets min-k to arg. Default=3. 
                         Scelestial considers all k-subsets of samples for min-k <= k <= max-k.
  -max-k arg             Sets max-k to arg. Default=4. 
                         Scelestial considers all k-subsets of samples for min-k <= k <= max-k.
  -include-root arg      Adds a sample with all sites equal to arg as root of the tree.
  -root root-index       Zero-based index of the new root.
                         After inferring the tree, tree edges are directed toward new root.
                         In the output line "u v w" u is the parent and v is the chine node.
  -no-internal-sample    Move all samples to leaf nodes.
                         If -root is also present, a neighbor of root-index is chosen as the root.
```


### Input format:
Input consists of a matrix with cells as columns and loci as rows. Each row starts with the name of the locus. Elements of the matrix represent sequencing results for one locus and one cell in a 10-state representation. Thus matrix elements are one of the following cases
A/A, T/T, C/C, G/G, A/T, A/C, A/G, T/C, T/G, C/G
Also missing value is represented as `./.`.

### Output format:
The output of scelestial represents the inferred evolutionary tree as well as the location of samples within the tree. The file starts with a number, the number of nodes in the inferred tree. Then, for each cell there is a line consists of two numbers and two strings. The first number is the id of the node, the second number shows if the node represents a node from the input (1) or not (0). The next string shows the original sequence of the sample and the next one is the imputed sequence for the sample.

After the description of the nodes of the tree, there is a part to represent the tree itself. For each edge there is a line consists of two endpoints of the edge followed by the length of the edge.

### Moving all the samples to the leaf nodes
Scelestial may assign samples to internal nodes. If the infered phylogeny is required to contain sample nodes as leaf nodes, the option `-no-internal-sample` could be used. With this option all the samples which are assigned to internal nodes are moved to new leaf nodes connected to their corresponding internal nodes. 

### Re-rooting the tree
The option `-root root-index` specifies the root of the tree. Scelestial internally infers an unrooted tree. With this option Scelestial move the root-index sample to the root of the tree via changing direction of the edges. The root index which is given to this option as an argument is the zero-based index of the input sample (column of the input matrix).

Note that if option `-no-internal-sample` is used in combination with option `-root root-index`, the root of the final tree would be parent of the sample root-index.

## Running Pre-defined Tests
### Test 2: A simple and small test

On the root of the folder execute following command:
```
sh bin/test/test1.sh
```

### Test 2: Rerooting
The option `-root 4` (in the bash script file) asks Scelestial to re-root the resulting tree to the sample 4 (which is 5-th column of the matrix).
On the root of the folder execute following command:
```
sh bin/test/test2.sh
```

### Test 3: Rerooting and moving samples to leaf nodes
With default options, Scelestial may place samples in internal nodes. In this representation simplifies the tree and thus the resulting visualization would be easier to understand. Scelestial has an option `-no-internal-sample` for moving samples to leaf nodes. In this case for each sample which is assigned to internal nodes, a new leaf node which is connected to the internal node with branch length zero is created and the samples is moved to this new node. 

The combination of `-no-internal-sample` option and `-root` option may make a confusion. Option `-root [ROOT]` is supposed to put sample `[ROOT]` at the root of the phylogeny, while option `-no-internal-sample` asks Scelestial to put all the samples in leaf nodes. When these two options are present, Scelestial puts parent of sample `[ROOT]` as the root of the tree. This test shows the result of executing Scelestial with both these two options.

On the root of the folder execute following command:
```
sh bin/test/test5.sh
```

### Other Tests
Some other tests are created for testing Scelestial. Details could be found in following files
* `bin/test/test3.sh`
* `bin/test/test4.sh`


## Generating simulated data:
The simulation is done by simulating an evolutionary tree. The tree starts with a node. The tree generation is done through some steps (number of steps is defined by --step parameter). At each step, a node from the evolutionary tree is selected and a child is added to it. Each node has a relative advantage. Nodes are selected based on their relative advantages. Advantage of a child is defined based on its parent advantage, it increases/decreases/remains as its parent with probabilities defined by parameters --aic/--adc/--akc after normalization (division by sum of these values). The amount of increment/decrement is defined by --ais and --ads parameters.

The tumor simulator accepts a lot of parameters. List of parameters are as follows:
```
  --sample arg                        number of sampled cells
  --step arg                          number of simulation steps
  --locus arg                         number of loci
  --aic arg                           advantage increase rate
  --adc arg                           advantage decrease rate
  --akc arg                           advantage keep rate
  --seed arg (=7)                     random seed number
  --seed-error arg (=7)               random seed for sequencing error
  --step-mutation-rate arg (=100)     average number of mutations in each step
  --ais arg (=1)                      advantage increase step
  --ads arg (=1)                      advantage decrease step
  --fclone arg                        clones file name
  --fseq arg                          seq file name
  --fseqtrue arg                      seq (without error) filename
  --ftree arg                         tree file name
  --mvrate arg (=0.20000000000000001) missing value rate
  --zorate arg (=0.10000000000000001) zero to one rate
  --ozrate arg (=0.20000000000000001) one to zero rate
  --remove-dup arg (=0)               remove duplicate sequences
```
The output consists of four files:
1. --fseq: the sequencing result (containing errors and missing values)
2. --fseqtrue: The true sequences for the sampled cells
3. --ftree: The simulated evolutionary tree
4. --fclone: Shows the simulation node for each sample.

A sample simulation could be
```bash
# Running the simulator
bin/synthesis --sample 5 --step 10 --locus 15 --aic 1 --adc 1 --akc 10 --ais 1 --ads 1 --seed 7 --fclone data/synth01-clone.txt --fseq data/synth01-seq.txt --fseqtrue data/synth01-seq-true.txt --ftree data/synth01-tree.txt --mvrate 0.5 --zorate 0.1 --ozrate 0.2
# Converting simulated data to scelestial's format
python src/convert-input.py data/synth01-seq.txt data/synth01-scelestial.txt /dev/null
# Generating some names for cells
for ((i=1; i<=5; i++)); do echo "C$i"; done > data/synth01-cell-names.txt
```

### Running Scelestial on Simulated Data:
The easiest part is to run the scelestial as follows:
```bash
bin/scelestial <data/synth01-scelestial.txt >data/synth01-scelestial-tree-clone.txt
```

### Evaluating the results:
To evaluate the results, first, we separate the tree from clones in the output of scelestia.
```bash
python src/steiner-to-clone-tree.py data/synth01-scelestial.txt data/synth01-scelestial-tree-clone.txt data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt
```

#### Comparing sample distances
Then we calculate distances between pairs of samples in both true tree and scelestial's tree.
```bash
# Generating the distance matrix for the true tree
python src/tree-clone-2-distance-matrix.py data/synth01-tree.txt data/synth01-clone.txt > data/synth01-true-distance-matrix.txt
# Generating the distance matrix for the inferred tree
python src/tree-clone-2-distance-matrix.py data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt > data/synth01-scelestial-distance-matrix.txt
```
Then we generate the distance matrix between true and scelestial's distance matrices as follows
```bash
python src/calc-matrix-distance.py data/synth01-true-distance-matrix.txt TRUE  data/synth01-scelestial-distance-matrix.txt SCEL
# Result:
# TRUE SCEL
# TRUE 0.00 0.41
# SCEL 0.41 0.00
```
The result of the previous command is a matrix. The only important element in the matrix is the element in the row TRUE and column SCEL (or row SCEL and column TRUE). This element shows an overall difference between distances between pairs of samples in two trues. This is a value between 0 and 2 and higher values show more difference between the trees.

#### Comparing partition similarity
```bash
# Generating partitions for the true tree
python src/tree-part-sim.py data/synth01-tree.txt data/synth01-clone.txt > data/synth01-true-part.txt
# Generating partitions for the inferred tree
python src/tree-part-sim.py data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt > data/synth01-scelestial-part.txt
```

Then we can compare partitions produced by the true tree and the inferred tree as follows
```bash
# Calculating partition similarity between partitions of the two trees
python src/calc-part-distance.py --match --normalize data/synth01-true-part.txt TRUE data/synth01-scelestial-part.txt SCEL
# Result:
# TRUE SCEL
# TRUE 1.0 0.6
# SCEL 0.6 1.0
```
The element in row TRUE and column SCEL shows the similarity between TRUE and SCEL with respect to the partitions they create on samples by the trees. The value is a number between 0 and 1.

### Generating PDF
A tree could be represented as PDF with the following command
```bash
python src/clone-tree-to-mu-tree-imput.py data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt data/synth01-seq.txt /dev/null data/synth01-cell-names.txt data/synth01-out --compress
```




## Reproducing results of paper

All the codes and scripts generating results appearing in the paper are in the "testing" folder. To perform the tests, some installations should be done.

### Preparation

To reproduce results in the paper, following applications should be installed for R:
```
library(oncoNEM)
library(igraph)
```

Following command should switch to conda environment with python 2 (for biopython):
```
conda py2
```

Also install java, boost-devel and graphviz for python. Makefile variables in the testing folder should be filled correctly.


### Test scripts

Testing scripts perform tests and put outputs to folders defined at the beginning of bash script files with name DIR. Several scripts run bin/test-gen-eval-all.sh with various parameters. Detail of tests are as follows:



| Script name | Location in paper | Output file |
| :--- | :--- | :--- |
| run-cmp30-1.sh | Fig. 7 (a) | result/TEST-CMP30-1/run-cmp30-1-avg.txt, result/TEST-CMP30-1/run-cmp30-1.txt |
| run-cmp30-17.sh | Table 1 | result/TEST-CMP30-17/run-cmp30-17-avg.txt, result/TEST-CMP30-17/run-cmp30-17.txt |
| run-cmp30-3b.sh  | Fig. 7 (b) | result/TEST-CMP30-3bs/run-cmp30-3b.txt, result/TEST-CMP30-3bs/run-cmp30-3b-avg.txt
| run-cmp30-5e3.sh | Fig. 1 | result/TEST-CMP30-5E3/run-cmp30-5e3.txt
| run-cmp30-7.sh | Fig. 4 | result/TEST-CMP30-7/[12]/*.pdf
| run-cmp30-8.sh | Fig. 5 | result/TEST-CMP30-8/1/*.pdf
| run-sens30-all4.sh | Fig. 2, Fig. 3 | result/TEST-SENS30-ALL4/LOCUS2.txt, <br/> result/TEST-SENS30-ALL4/MV2.txt, <br/> result/TEST-SENS30-ALL4/ZO2.txt, <br/> result/TEST-SENS30-ALL4/OZ2.txt, <br/> result/TEST-SENS30-ALL4/SAMPLE2.txt |


# License
This project is under GPL 3.0 license.
