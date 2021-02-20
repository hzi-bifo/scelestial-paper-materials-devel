#!/bin/bash
## Mimincs old run-sens-1-1.sh
#  run-sens-all2.sh -- run-sens-all2.sh & run-sens-all2.sh

########################################################
## FOR DEBUG

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

source init.sh

DIR=result/TEST-SENS30-ALL2
OUTPUT=run-sens30-all2
LEN=10

rm -rf $DIR
mkdir -p $DIR

OP_SET=( SAMPLE LOCUS MV ZO OZ ZOOZ )
FILE_NAME_SET=( sample loc mv zo oz zooz )
for idx in "${!OP_SET[@]}"; do
	OP=${OP_SET[$idx]}
	FILE_NAME=${FILE_NAME_SET[$idx]}
	if [ "$1" != "" ]; then
		if [ "$1" != "$OP" ]; then
			continue
		fi
	fi
	echo "RUN-SENS-ALL $OP $FILE_NAME"
	SAMPLE_OPT="50 50 1"
	if [ "$OP" = "SAMPLE" ]; then
		SAMPLE_OPT="5 200 20"
	fi
	LOC_OPT="200 200 1"
	if [ "$OP" == "LOCUS" ]; then
		LOC_OPT="20 1000 100"
	fi
	MV_OPT="0.5 0.51 1"
	if [ "$OP" == "MV" ]; then
		MV_OPT="0.0 0.61 0.1"
	fi
	ZO_OPT="0.000065 0.0000651 0.2"
	#ZO_OPT="0.000000 0.0000651 0.2"
	if [ "$OP" == "ZO" ]; then
		ZO_OPT="0.000000 1.0 0.1"
	fi
	OZ_OPT="0.33 0.331 0.1"
	#OZ_OPT="0.00 0.031 0.1"
	if [ "$OP" == "OZ" ]; then
		OZ_OPT="0.0 1.01 0.1"
	fi
		
	rm -rf $DIR/$OP
	mkdir -p $DIR/$OP
	$SCRIPT/test-gen-eval-all.sh $DIR/$OP -rep $LEN -mut 20 20.1 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 20 20 5 -step 5 5 1 -gen synth -run I &>$DIR/$OP/$OUTPUT.eo

	cat $DIR/$OP/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" locus "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s I?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pI?" | python $SCRIPT/eq-to-csv.py 		--avg impt `python $SCRIPT/repeat.py 'I?' 0 $LEN` 		--var vimpt `python $SCRIPT/repeat.py 'I?' 0 $LEN`	--avg5 pimpt `python $SCRIPT/repeat.py 'pI?' 0 $LEN` --var5 vpimpt `python $SCRIPT/repeat.py 'pI?' 0 $LEN` > $DIR/$OP/$OUTPUT.txt
	cat $DIR/$OP/$OUTPUT.txt > $DIR/"$OP"2.txt

done


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
