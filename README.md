# single-cell lineage tree reconstruction
To make
make

In short, to run scelestial
bin/scelestial <[input-file] >[output-file]
The detail of input/output formats and how to 
generate synthetic data follows.

input format:
Input consists of a matrix with cells as columns and loci as rows. Each row
starts with the name of the locus. Elements of matrix represent sequencing 
result for one locus and one cell in a 10-state representation. Thus matrix
elements are one of following cases
A/A, T/T, C/C, G/G, A/T, A/C, A/G, T/C, T/G, C/G
Also missing value is represented as ``./.''.

output format:
The output of scelestial represents the infered evolutionary tree as well as
location of samples within the tree. The file starts with a number, number of 
nodes in the infered tree. Then, for each cell there is a line consists of two numbers
and two strings. The first number is the id of the node, the second number shows
if the node represents a node from the input (1) or not (0). Next string
shows the original sequence of the sample and the next one is the imputed sequence for the sample.

After representation of the nodes of the tree, there is a part for representation of the tree itself.
For each edge there is a line consists of two endpoints of the edge followed by the length of the edge.


Generating simulated data:
The simulation is done by simulating an evolutionary tree. The tree starts with a
node. The tree generation is done through some steps (number of steps is defined
by --step parameter). At each step, a node from evolutionary tree is selected and a child is added to
it. Each node has a relative advantage. Nodes are selected based on their relative
advantages. Advantage of a child is defined based on its parent advantage, it 
increases/dicreases/remains as its parrent with probabilities defined by parameters
--aic/--adc/--akc after normalization (division by sum of these values). The amount of
increament/decrement is defined by --ais and --ads parameters.

The tumor simulator accepts a lot of parameters. List of parameters are as follows:
  --sample arg                        number of sampled cells
  --step arg                          number of simulation steps
  --locus arg                         number of locuses
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

The output constist of four files:
1) --fseq: the sequencing result (containing errors and missing values)
2) --fseqtrue: The true sequences for the sampled cells
3) --ftree: The simulated evolutionary tree
4) --fclone: Shows the simulation node for each sample.

A sample simulation could be
bin/synth --sample 5 --step 10 --locus 15 --aic 1 --adc 1 --akc 10 --ais 1 --ads 1 --seed 7 --fclone data/synth01-clone.txt --fseq data/synth01-seq.txt --fseqtrue data/synth01-seq-true.txt --ftree data/synth01-tree.txt --mvrate 0.5 --zorate 0.1 --ozrate 0.2 
python src/convert-input.py data/synth01-seq.txt data/synth01-scelestial.txt /dev/null
for ((i=1; i<=5; i++)); do echo "BC-$i"; done > data/synth01-cell-names.txt

Running Scelestial:
bin/scelestial <data/synth01-scelestial.txt >data/synth01-scelestial-tree-clone.txt

Evaluating the results:
Separating tree from clones in the output of scelestia.
python src/steiner-to-clone-tree.py data/synth01-scelestial.txt data/synth01-scelestial-tree-clone.txt data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt

The distance between two trees could be calculated as follows
Generating distance matrix for the true tree
python src/tree-clone-2-distance-matrix.py data/synth01-tree.txt data/synth01-clone.txt > data/synth01-true-distance-matrix.txt
Generating distance matrix for the infered tree
python src/tree-clone-2-distance-matrix.py data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt > data/synth01-scelestial-distance-matrix.txt
Calculating distances between matrices
python src/calc-matrix-distance.py data/synth01-true-distance-matrix.txt TRUE  data/synth01-scelestial-distance-matrix.txt SCEL

Generating partitions for the true tree
python src/tree-part-sim.py data/synth01-tree.txt data/synth01-clone.txt > data/synth01-true-part.txt
Generating partitions for the infered tree
python src/tree-part-sim.py data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt > data/synth01-scelestial-part.txt
Calculating partition similarity between partitions of the two trees
python src/calc-part-distance.py --match --normalize data/synth01-true-part.txt TRUE data/synth01-scelestial-part.txt SCEL

A tree could be represented as PDF with following command
python src/clone-tree-to-mu-tree-imput.py data/synth01-scelestial-tree.txt data/synth01-scelestial-clone.txt data/synth01-seq.txt /dev/null data/synth01-cell-names.txt data/synth01-out --compress
