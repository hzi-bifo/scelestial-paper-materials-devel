DIR=$1
SAMPLE_CNT=$2
LOCUS_CNT=$3
STEP_CNT=$4
MISSING_VALUE_RATE=$5
ZORATE=$6
OZRATE=$7
STEPMUTATIONRATE=$8
SEED=$9
REMOVE_DUP_OPT="${10} ${11}"
AIC="${12}"
ADC="${13}"
AKC="${14}"
AIS=0.1
ADS=0.1
OPTS="${15}"

echo " synthesis $SCRIPT/synthesis/synthesis --seed $SEED --seed-error $SEED --sample $SAMPLE_CNT --locus $LOCUS_CNT --aic $AIC --adc $ADC --akc $AKC --step $STEP_CNT --ais $AIS --ads $ADS --fclone $DIR/true-clones.txt --fseq $DIR/input-scite.txt --ftree $DIR/true-tree.txt --mvrate $MISSING_VALUE_RATE --zorate $ZORATE --ozrate $OZRATE --fseqtrue $DIR/true-seq.txt  --step-mutation-rate $STEPMUTATIONRATE $REMOVE_DUP_OPT $OPTS " 
$SCRIPT/synthesis/synthesis --seed $SEED --seed-error $SEED --sample $SAMPLE_CNT --locus $LOCUS_CNT --aic $AIC --adc $ADC --akc $AKC --step $STEP_CNT --ais $AIS --ads $ADS --fclone $DIR/true-clones.txt --fseq $DIR/input-scite.txt --ftree $DIR/true-tree.txt --mvrate $MISSING_VALUE_RATE --zorate $ZORATE --ozrate $OZRATE --fseqtrue $DIR/true-seq.txt  --step-mutation-rate $STEPMUTATIONRATE $REMOVE_DUP_OPT $OPTS
#Rscript gen-onconem.R `pwd` $DIR/input-scite.txt $DIR/true-clones.txt $DIR/true-tree.txt $DIR/true-seq.txt $SAMPLE $CLONE $LOCUS $((7+77777*CNT)) 
python $SCRIPT/convert-input.py $DIR/input-scite.txt $DIR/input-impute.txt $DIR/input-bp.txt 
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/true-tree.txt $DIR/true-clones.txt > $DIR/true-distance-matrix
#python tree-clone-uw-2-distance-matrix.py $DIR/true-tree.txt $DIR/true-clones.txt > $DIR/true-distance-matrix-uw
python $SCRIPT/tree-part-sim.py $DIR/true-tree.txt $DIR/true-clones.txt > $DIR/true-part
#python visualize-tree.py TRUE $DIR/true-tree.txt $DIR/true-clones.txt $DIR/output-true

#DEBUG:
echo "Gen Synthesis seed: $SEED $DIR $REMOVE_DUP_OPT"
head -n1 $DIR/input-scite.txt
