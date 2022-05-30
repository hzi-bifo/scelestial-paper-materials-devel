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
OUTPUT=run-cmp30-8
LEN=1
MATRIX_DISTANCE_NORMALIZATION=
TIMEOUT=

#rm -rf $DIR
mkdir -p $DIR

#SAMPLE_OPT="100 100 20"
#LOC_OPT="200 200 1"
#MV_OPT="0.05 0.051 1"
#ZO_OPT="0.00005 0.000051 0.2"
#OZ_OPT="0.18 0.181 0.1"
#./test-gen-eval-all.sh $DIR/$OP -rep $LEN -mut 10 10.01 4 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 0 1 10 -step 5 5 1 -gen data -run2 onco,steiner,bitphyl,scite,phiscs,sasc,sciphi,sifit,siclonefit -data-file ../../data/Luang/Luang-1.txt,../../data/Luang/Luang-2.txt &>$DIR/$OP/$OUTPUT.eo
#test-data-eval-all.sh $DIR -run2 onco,steiner,bitphyl,scite,phiscs,sasc,sciphi,sifit,siclonefit -data-file ../../data/Luang/Luang-1.txt,../../data/Luang/Luang-2.txt &>$DIR/$OUTPUT.eo
OLD_DIR=$DIR
SDIR=$DIR
DATA_FILES=data/Li-input.txt:data/Li-cell-types.txt:data/Li-cell-names.txt:data/Li-mutation-info.txt
#RUN2=,onco,steiner,bitphyl,scite,sasc,sciphi,sifit,siclonefit,
RUN2=,steiner,
#phiscs
CNT=0
for DATA_FILE_TYPE in ${DATA_FILES//,/ } ; do
	CNT=$((CNT+1))
	DIR=$SDIR/$CNT
	mkdir -p $DIR
	DATA_FILE=$(echo $DATA_FILE_TYPE | cut -f1 -d:)
	DATA_TYPE=$(echo $DATA_FILE_TYPE | cut -f2 -d:)
	DATA_NAME=$(echo $DATA_FILE_TYPE | cut -f3 -d:)
	DATA_MUT_INFO=$(echo $DATA_FILE_TYPE | cut -f4 -d:)
	cat $DATA_FILE > $DIR/input-scite.txt
	python $SCRIPT/convert-input.py $DIR/input-scite.txt $DIR/input-impute.txt $DIR/input-bp.txt 
	cp $DATA_NAME $DIR/cell-names.txt
	cp $DATA_MUT_INFO $DIR/mutation-info.txt
	echo "PARAMETERS: $DIR $DATA_FILE $DATA_TYPE"

#Running beast
		if [[ $RUN =~ .*E.* || $RUN2 =~ .*,beast,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-beast.sh $DIR $SAMPLE $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/ONC1-O-distance-matrix P$CNT $DIR/ONC1-BD-distance-matrix E$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/ONC1-O-part pP$CNT $DIR/ONC1-BD-part pE$CNT "
			OUT_STR="$OUT_STR Beast$CNT $DIR/ONCO1-O-tree.txt $DIR/ONCO1-O-clones.txt"
			echo "RUN_END BEAST"
		fi

# Running OncoNEM
		if [[ $RUN =~ .*O.* || $RUN2 =~ .*,onco,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-onconem.sh $DIR $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/onconeme-distance-matrix O$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/onconeme-part pO$CNT "
			OUT_STR="$OUT_STR Onco$CNT $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt"
			echo "RUN_END ONCONEM"
		fi

# Running Bermen
		if [[ $RUN =~ .*I.* || $RUN2 =~ .*,steiner,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-sciphyl.sh $DIR 3 $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/impute-distance-matrix I$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/impute-part pI$CNT "
			OUT_STR="$OUT_STR Steiner$CNT $DIR/impute-tree.txt $DIR/impute-clones.txt"
			echo "RUN_END STEINER"
		fi

# Running BitPhylogeny
		if [[ $RUN =~ .*B.* || $RUN2 =~ .*,bitphyl,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-bp.sh $DIR $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/bp-distance-matrix B$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/bp-part pB$CNT "
			OUT_STR="$OUT_STR Bitphyl$CNT $DIR/bp-tree.txt $DIR/bp-clones.txt"
			echo "RUN_END BP"
		fi

#Running SCITE
		if [[ $RUN =~ .*S.* || $RUN2 =~ .*,scite,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-scite.sh $DIR $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/scite-distance-matrix S$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/scite-part pS$CNT "
			OUT_STR="$OUT_STR SCITE$CNT $DIR/scite-tree.txt $DIR/scite-clones.txt"
			echo "RUN_END SCITE"
		fi

#Running BSC [Beast Single-cell Plugin]
		if [[ $RUN =~ .*U.* || $RUN2 =~ .*,scplug,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-bsc.sh $DIR $BSCFOLDER $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/bsc-distance-matrix U$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/bsc-part pU$CNT "
			OUT_STR="$OUT_STR SCPlug$CNT $DIR/bsc-tree.txt $DIR/bsc-clones.txt"
			echo "RUN_END SCPLUG"
		fi

#Running BSC [Beast Single-cell Plugin]
		if [[ $RUN =~ .*L.* || $RUN2 =~ .*,scplugwl,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-bsc-loc.sh $DIR $BSCFOLDER $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/bscloc-distance-matrix L$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/bscloc-part pL$CNT "
			OUT_STR="$OUT_STR SCPlugLoc$CNT $DIR/bscloc-tree.txt $DIR/bscloc-clones.txt"
			echo "RUN_END SCPLUGWL"
		fi

#Running phiscs
		if [[ $RUN2 =~ .*,phiscs,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-phiscs.sh $DIR 600 $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/phiscs-distance-matrix phiscs$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/phiscs-part pphiscs$CNT "
			OUT_STR="$OUT_STR PhISCS$CNT $DIR/phiscs-tree.txt $DIR/phiscs-clones.txt"
			echo "RUN_END PHISCS"
		fi

#Running sasc
		if [[ $RUN2 =~ .*,sasc,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-sasc.sh $DIR $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/sasc-distance-matrix sasc$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/sasc-part psasc$CNT "
			OUT_STR="$OUT_STR SASC$CNT $DIR/sasc-tree.txt $DIR/sasc-clones.txt"
			echo "RUN_END SASC"
		fi

#Running sciphi
		if [[ $RUN2 =~ .*,sciphi,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-sciphi.sh $DIR $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/sciphi-distance-matrix sciphi$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/sciphi-part psciphi$CNT "
			OUT_STR="$OUT_STR SCIPhI$CNT $DIR/sciphi-tree.txt $DIR/sciphi-clones.txt"
			echo "RUN_END SCIPHI"
		fi

#Running sifit
		if [[ $RUN2 =~ .*,sifit,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-sifit.sh $DIR 10000 $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/sifit-distance-matrix sifit$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/sifit-part psifit$CNT "
			OUT_STR="$OUT_STR SiFit$CNT $DIR/sifit-tree.txt $DIR/sifit-clones.txt"
			echo "RUN_END SIFIT"
		fi

#Running siclonefit
		if [[ $RUN2 =~ .*,siclonefit,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-siclonefit.sh $DIR 1 $MATRIX_DISTANCE_NORMALIZATION 2>&1
			CMP_STR="$CMP_STR $DIR/siclonefit-distance-matrix siclonefit$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/siclonefit-part psiclonefit$CNT "
			OUT_STR="$OUT_STR SiCloneFit$CNT $DIR/siclonefit-tree.txt $DIR/siclonefit-clones.txt"
			echo "RUN_END SICLONEFIT"
		fi

	for b in 'impute' 'bp' 'onconeme' 'sasc' 'sciphi' 'scite' 'siclonefit' 'sifit' ; do
		echo $a $b:
		python $SCRIPT/visualize-tree.py $b $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/output-$b --type-file $DATA_TYPE; 
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info.txt $DIR/cell-names.txt $DIR/output-$b-mut2 --compress 

		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info.txt $DIR/cell-names.txt $DIR/output-$b-mut-cmmall --mark-mutations common=1:all:ignore-leaf
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info.txt $DIR/cell-names.txt $DIR/output-$b-mut-cmmpos --mark-mutations common=1:pos:ignore-leaf
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info.txt $DIR/cell-names.txt $DIR/output-$b-mut-cmm0.9all --mark-mutations common=0.9:all:ignore-leaf
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info.txt $DIR/cell-names.txt $DIR/output-$b-mut-dcpos --mark-mutations direct-children:positive:ignore-leaf
	done

done > $SDIR/$OUTPUT.eo

DIR=$OLD_DIR


echo "Done" 

awk -F '\t' 'BEGIN { OFS="\t" }{$12 = "-"; print}' $DIR/mutation-info.txt > $DIR/mutation-info-empty.txt

CNT=0
for DATA_FILE_TYPE in ${DATA_FILES//,/ } ; do
	CNT=$((CNT+1))
	DIR=$SDIR/$CNT
	DATA_FILE=$(echo $DATA_FILE_TYPE | cut -f1 -d:)
	DATA_TYPE=$(echo $DATA_FILE_TYPE | cut -f2 -d:)
	DATA_NAME=$(echo $DATA_FILE_TYPE | cut -f3 -d:)
	DATA_MUT_INFO=$(echo $DATA_FILE_TYPE | cut -f4 -d:)

	for b in 'impute'; do
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info-empty.txt $DIR/cell-names.txt $DIR/output-$b-mut-cmmall --mark-mutations common=1:all:ignore-leaf
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info-empty.txt $DIR/cell-names.txt $DIR/output-$b-mut-cmmpos --mark-mutations common=1:pos:ignore-leaf
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info-empty.txt $DIR/cell-names.txt $DIR/output-$b-mut-cmm0.9all --mark-mutations common=0.9:all:ignore-leaf
		python $SCRIPT/clone-tree-to-mu-tree-imput.py $DIR/$b-tree.txt $DIR/$b-clones.txt $DIR/input-scite.txt $DIR/mutation-info-empty.txt $DIR/cell-names.txt $DIR/output-$b-mut-dcpos --mark-mutations direct-children:positive:ignore-leaf
	done
done
