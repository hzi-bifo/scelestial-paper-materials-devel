# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

source init.sh

DIR=result/TEST-SENS30-ALL4P
OUTPUT=run-sens30-all4p
LEN=10
OPT=$1

rm -rf $DIR
mkdir -p $DIR

OP_SET=( SAMPLE LOCUS MV ZO OZ )
FILE_NAME_SET=( sample loc mv zo oz )

for REP in `seq 1 $LEN`; do 
#export -f run
#for idx in "${!OP_SET[@]}"; do
	#OP=${OP_SET[$idx]}
	#FILE_NAME=${FILE_NAME_SET[$idx]}

	#if [ "$OPT" != "" ]; then
	#	if [ "$OPT" != "$OP" ]; then
	#		continue
	#	fi
	#fi
	#echo "RUN-SENS-ALL $OP $FILE_NAME ($idx)"

	#SAMPLE_OPT="50 50 1"
	#if [ "$OP" = "SAMPLE" ]; then
		SAMPLE_OPT="5 200 20"
	#fi
	#LOC_OPT="200 200 1"
	#if [ "$OP" == "LOCUS" ]; then
		LOC_OPT="20 1000 100"
	#fi
	#MV_OPT="0.07 0.071 1"
	#if [ "$OP" == "MV" ]; then
		MV_OPT="0.00 0.501 0.1"
	#fi
	#ZO_OPT="0.015 0.0151 0.1"
	##ZO_OPT="0.000000 0.0000651 0.2"
	#if [ "$OP" == "ZO" ]; then
		ZO_OPT="0.000 0.1001 0.01"
	#fi
	#OZ_OPT="0.10 0.101 0.1"
	#OZ_OPT="0.00 0.031 0.1"
	#if [ "$OP" == "OZ" ]; then
		OZ_OPT="0.000 0.500 0.10"
	#fi
		
	#rm -rf $DIR/$OP
	#mkdir -p $DIR/$OP
	mkdir $DIR/$REP
	sem -j 10 $SCRIPT/test-gen-eval-all.sh $DIR/$REP -rep 1 -mut 20 20.1 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 20 20 5 -step 5 5 1 -gen synth -run I -seed $((77777*REP)) -timeout 10h 36000 &>$DIR/$REP/$OUTPUT.eo

done
sem --wait


for REP in `seq 1 $LEN`; do 
	cat $DIR/$REP/$OUTPUT.eo
done > $DIR/$OUTPUT.eo


cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" locus "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s I?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pI?" | python $SCRIPT/eq-to-csv.py --avg impt `python $SCRIPT/repeat.py 'I?' 0 $LEN` --var vimpt `python $SCRIPT/repeat.py 'I?' 0 $LEN`	--avg5 pimpt `python $SCRIPT/repeat.py 'pI?' 0 $LEN` --var5 vpimpt `python $SCRIPT/repeat.py 'pI?' 0 $LEN` > $DIR/$OUTPUT.txt
#cat $DIR/$OUTPUT.txt > $DIR/2.txt


exit

#Genereate $OUTPUT-1.eo (before everything is finished)
LEN=10

STEP_MIN=5
STEP_MAX=5
STEP_INC=1
DATA_FILES=NONE
CLONE_MIN=20
CLONE_MAX=20
CLONE_INC=1
SAMPLE_MIN=5
SAMPLE_MAX=200
SAMPLE_INC=20

LOCUS_MIN=20
LOCUS_MAX=1000
LOCUS_INC=100
MISSING_VALUE_RATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange(0.00, 0.501, 0.1)]));"`
ZORATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange(0.000, 0.1001, 0.01)]));"`
OZRATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange(0.000, 0.500, 0.10)]));"`
STEPMUTATIONRATE_RANGE="20"
MATRIX_MATRIX_DISTANCE_NORMALIZATION=
REPEAT_COUNT=1
for REP in `seq 1 $LEN`; do 
	ls $DIR/$REP | while read a; do
ALLCNT=0
for DATA_FILE in ${DATA_FILES//,/ } ; do
for ((STEP=STEP_MIN; STEP<=STEP_MAX; STEP+=STEP_INC)); do
#ONLY FOR ONCONEM-GEN
for ((CLONE=CLONE_MIN; CLONE<=CLONE_MAX; CLONE+=CLONE_INC)); do
for ((STEP=STEP_MIN; STEP<=STEP_MAX; STEP+=STEP_INC)); do
#for ((CLONE=CLONE_MIN; CLONE<=CLONE_MAX; CLONE+=CLONE_INC)); do
for ((SAMPLE=SAMPLE_MIN; SAMPLE<=SAMPLE_MAX; SAMPLE+=SAMPLE_INC)); do
for ((LOCUS=LOCUS_MIN; LOCUS<=LOCUS_MAX; LOCUS+=LOCUS_INC)); do

for MISSING_VALUE_RATE in $MISSING_VALUE_RATE_RANGE; do
for ZORATE in $ZORATE_RANGE; do
for OZRATE in $OZRATE_RANGE; do

#for ((STEPMUTATIONRATE=STEPMUTATIONRATE_MIN; STEPMUTATIONRATE<=STEPMUTATIONRATE_MAX; STEPMUTATIONRATE+=STEPMUTATIONRATE_INC)); do
for STEPMUTATIONRATE in $STEPMUTATIONRATE_RANGE; do
	ALLCNT=$((ALLCNT+1))
	echo "DIR=$DIR SAMPLE=$SAMPLE LOCUS=$LOCUS STEP=$STEP MVR=$MISSING_VALUE_RATE 0->1=$ZORATE 1->0:$OZRATE STEPMUT:$STEPMUTATIONRATE CLONE:$CLONE ALLCNT:$ALLCNT DIV_RATE=$DIV_RATE DEATH_RATE=$DEATH_RATE ITER=$REP"
	for ((CNT=0; CNT<REPEAT_COUNT; CNT++)); do
		IDIR=$DIR/$REP/$ALLCNT/$CNT
		CMP_STR="$IDIR/true-distance-matrix True$CNT "
		CMP_PART_STR="$IDIR/true-part pTrue$CNT "
		OUT_STR="True$CNT $IDIR/true-tree.txt $IDIR/true-clones.txt"
		CMP_STR="$CMP_STR $IDIR/impute-distance-matrix I$CNT "
		CMP_PART_STR="$CMP_PART_STR $IDIR/impute-part pI$CNT "
		OUT_STR="$OUT_STR Steiner$CNT $IDIR/impute-tree.txt $IDIR/impute-clones.txt"
		python $SCRIPT/calc-matrix-distance.py $MATRIX_MATRIX_DISTANCE_NORMALIZATION $CMP_STR
		python $SCRIPT/calc-part-distance.py --match --normalize $CMP_PART_STR
	done
done
done; done; done; done; done; done; done; done; done; 

done; done > $DIR/$OUTPUT-1.eo

LEN=1
cat $DIR/$OUTPUT-1.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" locus "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s I?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pI?" | grep 'pI0:' | grep 'pI0:' | python $SCRIPT/eq-to-csv.py --avg impt `python $SCRIPT/repeat.py 'I?' 0 $LEN` --var vimpt `python $SCRIPT/repeat.py 'I?' 0 $LEN`	--avg5 pimpt `python $SCRIPT/repeat.py 'pI?' 0 $LEN` --var5 vpimpt `python $SCRIPT/repeat.py 'pI?' 0 $LEN` > $DIR/$OUTPUT.txt
#Then create tables via Libre Office



exit

for ((CNT=1; CNT<=11; CNT++)); do 
for ((RCNT=0; RCNT<10; RCNT++)); do
	LDIR=$DIR/MV/$CNT/$RCNT/
	for ((i=0; i<50; i++)); do echo C-$i; done > $LDIR/cell-names.txt
	for ((i=0; i<50; i++)); do echo 1; done > $LDIR/cell-types.txt
	touch  $LDIR/mutation-info.txt
	for b in impute true ; do
		python $SCRIPT/visualize-tree.py $b $LDIR/$b-tree.txt $LDIR/$b-clones.txt $LDIR/output-$b --type-file $LDIR/cell-types.txt; 
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $LDIR/$b-tree.txt $LDIR/$b-clones.txt $LDIR/input-scite.txt $LDIR/mutation-info.txt $LDIR/cell-names.txt $LDIR/output-$b-mut2 --compress 
		cp $LDIR/output-$b-mut2.pdf $DIR/MV/output-$b-mut2-$CNT-$RCNT.pdf
	done
done
done
