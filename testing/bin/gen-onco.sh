DIR=$1
SAMPLE_CNT=$2
CLONE_CNT=$3
LOCUS_CNT=$4
SEED=$5
UNOBSERVED=$6
FPR=$7
FNR=$8
MVR=$9

Rscript $SCRIPT/gen-onconem-new.R `pwd` $DIR/input-scite.txt $DIR/true-clones.txt $DIR/true-tree.txt $DIR/true-seq.txt $SAMPLE_CNT $CLONE_CNT $LOCUS_CNT $SEED $UNOBSERVED $FPR $FNR $MVR
python $SCRIPT/convert-input.py $DIR/input-scite.txt $DIR/input-impute.txt $DIR/input-bp.txt 
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/true-tree.txt $DIR/true-clones.txt > $DIR/true-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/true-tree.txt $DIR/true-clones.txt > $DIR/true-part
