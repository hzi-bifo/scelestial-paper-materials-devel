#!/bin/bash
## Mimincs old run-cmp6-all12-3-5


# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

source init.sh

DIR=result/TEST-CMP30-7P
DIRIN=result/TEST-CMP30-7P/test/
OUTPUT=run-cmp30-7-sample-on-leaf
LEN=1
MATRIX_DISTANCE_NORMALIZATION=
TIMEOUT=

#rm -rf $DIR
#mkdir -p $DIR $DIR/1 $DIR/2

OLD_DIR=$DIR
SDIR=$DIR
DATA_FILES=data/Luang/Luang-1.txt:data/Luang/Luang-1-types.txt:data/Luang/Luang-1-names.txt:data/Luang/Luang-1-mut-info.txt:N-37,data/Luang/Luang-2.txt:data/Luang/Luang-2-types.txt:data/Luang/Luang-2-names.txt:data/Luang/Luang-2-mut-info.txt:N-27
#RUN2=,onco,steiner,bitphyl,scite,sasc,sciphi,sifit,siclonefit,
RUN2=,steiner,
#phiscs,

SEM=""
SEM_WAIT=""


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
	#cp $DATA_TYPE $DIR/cell-type.txt
	echo "PARAMETERS: $DIR $DATA_FILE $DATA_TYPE"

	if [[ $RUN =~ .*I.* || $RUN2 =~ .*,steiner,.* ]]; then
		timestamp
		time $SCRIPT/scelestial -min-k 3 -max-k 3 -no-internal-sample -root $DATA_ROOT_ZBINDEX < $DIR/input-impute.txt > $DIRTEST/impute-tree-clone.txt 
		python $SCRIPT/steiner-to-clone-tree.py $DIR/input-scite.txt $DIRTEST/impute-tree-clone.txt $DIRTEST/impute-tree.txt $DIRTEST/impute-clones.txt 2>/dev/null
		python $SCRIPT/steiner-to-seq.py $DIR/input-scite.txt $DIRTEST/impute-tree-clone.txt > $DIRTEST/impute-seq.txt
		#$SEM $TIMEOUT $SCRIPT/run-sciphyl.sh $DIR 3 $MATRIX_DISTANCE_NORMALIZATION 2>&1
		#CMP_STR="$CMP_STR $DIR/impute-distance-matrix I$CNT "
		#CMP_PART_STR="$CMP_PART_STR $DIR/impute-part pI$CNT "
		#OUT_STR="$OUT_STR Steiner$CNT $DIR/impute-tree.txt $DIR/impute-clones.txt"
		echo "RUN_END STEINER"
	fi

	echo "waiting for semaphores ... "
	$SEM_WAIT
	echo "waiting for semaphores done"

	#for b in 'impute' 'bp' 'onconeme' 'sasc' 'sciphi' 'scite' 'siclonefit' 'sifit' ; do
	for b in 'impute' ; do
		if [ -f $DIRTEST/$b-tree.txt ]; then
			echo $a $b:
			echo -e ""
			echo python $SCRIPT/visualize-tree.py $b $DIRTEST/$b-tree.txt $DIRTEST/$b-clones.txt $DIRTEST/output-$b --type-file $DATA_TYPE --cell-name $DIR/cell-names.txt --color $DIR/colors.txt 
			python $SCRIPT/visualize-tree.py $b $DIRTEST/$b-tree.txt $DIRTEST/$b-clones.txt $DIRTEST/output-$b --type-file $DATA_TYPE --cell-name $DIR/cell-names.txt --color $DIR/colors.txt 
			python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIRTEST/$b-tree.txt $DIRTEST/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info.txt $DIR/cell-names.txt $DIRTEST/output-$b-mut2 --compress --cell-type $DIR/cell-type.txt

			python $SCRIPT/tree-to-fig-tree.py $DIRTEST/$b-tree.txt $DIRTEST/$b-clones.txt $DIR/input-scite.txt $DIRTEST/$b-tree
		else
			echo "Warnning: no $b files"
		fi
	done
done > $SDIR/$OUTPUT.eo

DIR=$OLD_DIR


echo "Done" 


CNT=0
for DATA_FILE_TYPE in `echo $DATA_FILES | sed 's/,/\n/g'`  ; do
	CNT=$((CNT+1))
	DIR=$SDIR/$CNT
	DATA_FILE=$(echo $DATA_FILE_TYPE | cut -f1 -d:)
	DATA_TYPE=$(echo $DATA_FILE_TYPE | cut -f2 -d:)
	DATA_NAME=$(echo $DATA_FILE_TYPE | cut -f3 -d:)
	DATA_MUT_INFO=$(echo $DATA_FILE_TYPE | cut -f4 -d:)
	for b in 'impute'; do
		python $SCRIPT/visualize-tree.py $b $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/output-$b-mut-cmmall --type-file $DATA_TYPE --cell-name $DIR/cell-names.txt --color $DIR/colors.txt --seq $DIR/input-scite.txt --mark-mutations common=1:all:ignore-leaf
		python $SCRIPT/visualize-tree.py $b $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/output-$b-mut-cmmpos --type-file $DATA_TYPE --cell-name $DIR/cell-names.txt --color $DIR/colors.txt --seq $DIR/input-scite.txt --mark-mutations common=1:positive:ignore-leaf
		python $SCRIPT/visualize-tree.py $b $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/output-$b-mut-cmm0.9all --type-file $DATA_TYPE --cell-name $DIR/cell-names.txt --color $DIR/colors.txt --seq $DIR/input-scite.txt --mark-mutations common=0.9:all:ignore-leaf
		#mutations common in a node and its direct children
		python $SCRIPT/visualize-tree.py $b $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/output-$b-mut-dcpos --type-file $DATA_TYPE --cell-name $DIR/cell-names.txt --color $DIR/colors.txt --seq $DIR/input-scite.txt --mark-mutations direct-children:positive:ignore-leaf 
	done
done

