#!/bin/bash
cd /home/hforoughmand/sc-phyl-infer-steiner/test/test-setup-1

source /home/hforoughmand/miniconda3/etc/profile.d/conda.sh
conda activate py3
clear

LIB=L21
NAME=ss
DIR=TEST-L21-short

#rm -rf $DIR
#mkdir $DIR

#GENERATE
echo "Generating inputs, or copying it ..."
DATA=../../data/
#cat $DATA/$LIB-filter.csv | python $DATA/cleaner/filter-matrix.py --lib $LIB --cell-type-file $DATA/Larvae_Seurat_batch_r_out_cells_3.csv --class-list 0-4 --cell-filter 0.5 --row-filter 0.8 --scite-out $DIR/input-scite.txt --impute-out $DIR/input-impute.txt --type-out $DIR/true-type.txt --class-resample-cnt 10

#GENEREATE - CONTINUED
python convert-input.py $DIR/input-scite.txt $DIR/input-impute-useless.txt $DIR/input-bp.txt
# python tree-clone-2-distance-matrix.py $DIR/true-tree.txt $DIR/true-clones.txt > $DIR/true-distance-matrix
# python visualize-tree.py TRUE $DIR/true-tree.txt $DIR/true-clones.txt $DIR/output-true

# Running OncoNEM
echo "Running oncoNEM ..."
time /usr/bin/Rscript run-onconeme.R $DIR/input-scite.txt `pwd` $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt
python tree-clone-2-distance-matrix.py $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt > $DIR/onconeme-distance-matrix
python visualize-tree.py Onco $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt $DIR/output-onconeme


# Running Bermen
echo "Running Bermen's alg ..."
time ../../phyl/impute -min-k 3 -max-k 3 < $DIR/input-impute.txt > $DIR/impute-tree-clone.txt 
python steiner-to-clone-tree.py $DIR/input-scite.txt $DIR/impute-tree-clone.txt $DIR/impute-tree.txt $DIR/impute-clones.txt
python tree-clone-2-distance-matrix.py $DIR/impute-tree.txt $DIR/impute-clones.txt > $DIR/impute-distance-matrix
python visualize-tree.py Impute $DIR/impute-tree.txt $DIR/impute-clones.txt $DIR/output-impute --type-file $DIR/true-type.txt
#python summ-classification.py 


# Running BitPhylogeny
echo "Running BitPhylogeny ..."
conda activate py2
time python ../../bitphylogeny/python/bin/BitPhylogeny analyse $DIR/input-bp.txt $DIR/temp/bit-output -n 200 -b 10 -t 5 -mode "mutation" -seed 1234
conda deactivate
python bitphylogeny-2-tree-clone.py $DIR/temp/bit-output/input-bp.txt/treescripts/nodes-`ls $DIR/temp/bit-output/input-bp.txt/treescripts/nodes-*.graphml | sed 's/^.*-\([0-9]*\)\.graphml$/\1/' | sort -nr | head -n 1`.graphml $DIR/bp-tree.txt $DIR/bp-clones.txt
python tree-clone-2-distance-matrix.py $DIR/bp-tree.txt $DIR/bp-clones.txt > $DIR/bp-distance-matrix
python visualize-tree.py BPhyl $DIR/bp-tree.txt $DIR/bp-clones.txt $DIR/output-bp --type-file $DIR/true-type.txt

#Running SCITE
echo "Running SCITE ..."
time ../../SCITE/scite -i $DIR/input-scite.txt -n `cat $DIR/input-scite.txt | wc -l` -m `head -1 $DIR/input-scite.txt | sed 's/[^,]//g' | wc -c` -r 1 -l 200000 -fd 0.2 -ad 0.21545 0.1 -o $DIR/temp/scite- -max_treelist_size 1 -a -seed 7
python scite-2-tree-clone.py $DIR/temp/scite-_ml0.gv $DIR/input-scite.txt $DIR/scite-tree.txt $DIR/scite-clones.txt
python tree-clone-2-distance-matrix.py $DIR/scite-tree.txt $DIR/scite-clones.txt > $DIR/scite-distance-matrix
python visualize-tree.py SCITE $DIR/scite-tree.txt $DIR/scite-clones.txt $DIR/output-scite


# COMPARISON
echo "Comparison ..."
python calc-matrix-distance.py $DIR/onconeme-distance-matrix Onco $DIR/impute-distance-matrix Impt $DIR/bp-distance-matrix BPhl $DIR/scite-distance-matrix SCTE

echo -n "OncoNEM:"
python tree-crossing-number.py $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt $DIR/true-type.txt
echo -n "SCITE:"
python tree-crossing-number.py $DIR/scite-tree.txt $DIR/scite-clones.txt $DIR/true-type.txt
echo -n "BP:"
python tree-crossing-number.py $DIR/bp-tree.txt $DIR/bp-clones.txt $DIR/true-type.txt
echo -n "impute:"
python tree-crossing-number.py $DIR/impute-tree.txt $DIR/impute-clones.txt $DIR/true-type.txt
