#!/bin/bash
## Mimincs old run-cmp6-all12-3-5

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

source init.sh

DIR=result/TEST-CMP30-8
DIRIN=result/TEST-CMP30-8/test/
OUTPUT=run-cmp30-8-sample-on-leaf
LEN=1
MATRIX_DISTANCE_NORMALIZATION=
TIMEOUT=

#rm -rf $DIR
#mkdir -p $DIR

#SAMPLE_OPT="100 100 20"
#LOC_OPT="200 200 1"
#MV_OPT="0.05 0.051 1"
#ZO_OPT="0.00005 0.000051 0.2"
#OZ_OPT="0.18 0.181 0.1"
#./test-gen-eval-all.sh $DIR/$OP -rep $LEN -mut 10 10.01 4 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 0 1 10 -step 5 5 1 -gen data -run2 onco,steiner,bitphyl,scite,phiscs,sasc,sciphi,sifit,siclonefit -data-file ../../data/Luang/Luang-1.txt,../../data/Luang/Luang-2.txt &>$DIR/$OP/$OUTPUT.eo
#test-data-eval-all.sh $DIR -run2 onco,steiner,bitphyl,scite,phiscs,sasc,sciphi,sifit,siclonefit -data-file ../../data/Luang/Luang-1.txt,../../data/Luang/Luang-2.txt &>$DIR/$OUTPUT.eo
OLD_DIR=$DIR
SDIR=$DIR
DATA_FILES=data/Li-input.txt:data/Li-cell-types.txt:data/Li-cell-names.txt:data/Li-mutation-info.txt:BN-3
#RUN2=,onco,steiner,bitphyl,scite,sasc,sciphi,sifit,siclonefit,
RUN2=,steiner,
#phiscs
CNT=0
for DATA_FILE_TYPE in ${DATA_FILES//,/ } ; do
	CNT=$((CNT+1))
	DIR=$SDIR/$CNT
	DIRTEST=$DIRIN/$CNT
	mkdir -p $DIRTEST
	#mkdir -p $DIR
	DATA_FILE=$(echo $DATA_FILE_TYPE | cut -f1 -d:)
	DATA_TYPE=$(echo $DATA_FILE_TYPE | cut -f2 -d:)
	DATA_NAME=$(echo $DATA_FILE_TYPE | cut -f3 -d:)
	DATA_MUT_INFO=$(echo $DATA_FILE_TYPE | cut -f4 -d:)
	DATE_ROOT_NAME=$(echo $DATA_FILE_TYPE | cut -f5 -d:)
	DATA_ROOT_INDEX=$(grep -n $DATE_ROOT_NAME $DATA_NAME | cut -d: -f1)
	DATA_ROOT_ZBINDEX=$(( DATA_ROOT_INDEX-1 ))
	#cat $DATA_FILE > $DIR/input-scite.txt
	#python $SCRIPT/convert-input.py $DIR/input-scite.txt $DIR/input-impute.txt $DIR/input-bp.txt 
	#cp $DATA_NAME $DIR/cell-names.txt
	#cp $DATA_MUT_INFO $DIR/mutation-info.txt
	echo "PARAMETERS: $DIR $DATA_FILE $DATA_TYPE"

	if [[ $RUN =~ .*I.* || $RUN2 =~ .*,steiner,.* ]]; then
		timestamp
		echo time $SCRIPT/scelestial -min-k 3 -max-k 3 -no-internal-sample -root $DATA_ROOT_ZBINDEX $DIR/input-impute.txt $DIRTEST/impute-tree-clone.txt 
		time $SCRIPT/scelestial -min-k 3 -max-k 3 -no-internal-sample -root $DATA_ROOT_ZBINDEX < $DIR/input-impute.txt > $DIRTEST/impute-tree-clone.txt 
		python $SCRIPT/steiner-to-clone-tree.py $DIR/input-scite.txt $DIRTEST/impute-tree-clone.txt $DIRTEST/impute-tree.txt $DIRTEST/impute-clones.txt 
		python $SCRIPT/steiner-to-seq.py $DIR/input-scite.txt $DIRTEST/impute-tree-clone.txt > $DIRTEST/impute-seq.txt
		#$TIMEOUT $SCRIPT/run-sciphyl.sh $DIR 3 $MATRIX_DISTANCE_NORMALIZATION 2>&1
		#CMP_STR="$CMP_STR $DIR/impute-distance-matrix I$CNT "
		#CMP_PART_STR="$CMP_PART_STR $DIR/impute-part pI$CNT "
		#OUT_STR="$OUT_STR Steiner$CNT $DIR/impute-tree.txt $DIR/impute-clones.txt"
		echo "RUN_END STEINER"
	fi

	#for b in 'impute' 'bp' 'onconeme' 'sasc' 'sciphi' 'scite' 'siclonefit' 'sifit' ; do
	for b in 'impute'; do
		echo $b:
		python $SCRIPT/visualize-tree.py $b $DIRTEST/$b-tree.txt $DIRTEST/$b-clones.txt $DIRTEST/output-$b --type-file $DATA_TYPE; 
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIRTEST/$b-tree.txt $DIRTEST/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info.txt $DIR/cell-names.txt $DIRTEST/output-$b-mut2 --compress 
	done

done > $SDIR/$OUTPUT.eo

DIR=$OLD_DIR


echo "Done" 
