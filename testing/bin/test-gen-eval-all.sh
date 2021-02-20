#!/bin/bash

set -u

RM=false
SERVER=false

STEP_MIN=100
STEP_MAX=100
STEP_INC=100
CLONE_MIN=20
CLONE_MAX=20
CLONE_INC=20
SAMPLE_MIN=50
SAMPLE_MAX=50
SAMPLE_INC=50
LOCUS_MIN=100
LOCUS_MAX=100
LOCUS_INC=100
MISSING_VALUE_RATE_RANGE=0.5
ZORATE_RANGE=0.185
OZRATE_RANGE=0.08
REPEAT_COUNT=10
STEPMUTATIONRATE_MIN=20
STEPMUTATIONRATE_MAX=20
STEPMUTATIONRATE_INC=20
STEPMUTATIONRATE_RANGE=20
UNOBSERVED=0
GEN=
#E:Beast O:Onconem S:Scite I:Impute B:BitPhylogeny U: SCBeast without locaion L: SCBeast with location
RUN=
RUN2=
SEED=0
KEEPTESTS=false
REMOVE_DUP_OPT=" --remove-dup false "
AIC=1
ADC=2
AKC=10
DEATH_RATE=0.5
DIV_RATE=1.5
OBS_PROB=0.001
ROOT_HEIGHT=4
MOTION_STEP_LENGTH=1
BSCFOLDER=./
GEN_SC_SAMPLE_BOUND_LB=0
GEN_SC_SAMPLE_BOUND_UB=1000000
DATA_FILES=NONE
TIMEOUT="timeout 15m "
TIMEOUT_SECONTS=880
#TIMEOUT=
MATRIX_MATRIX_DISTANCE_NORMALIZATION=
MATRIX_DISTANCE_NORMALIZATION="-m-distance-norm-no"
EXIT_ON_ERROR=1
OPT_GEN=
RUN_OPT_SICLONEFIT=
RUN_OPT_SIFIT=

SDIR=$1
shift

while :
do
	case "$1" in
	-gen)
		GEN=$2
		if [ "$GEN" = "onco" ]; then
			:
		elif [ "$GEN" = "synth" ]; then
			:
		elif [ "$GEN" = "spsimul-load" ]; then
			DIRFILES=$3
			shift 1
			:
		elif [ "$GEN" = "spsimul" ]; then
			:
		elif [ "$GEN" = "data" ]; then
			:
		elif [ "$GEN" = "none" ]; then
			:
		else
			echo "Invalid generation method $GEN"
			exit 1;
		fi
		shift 2
		;;
	-rm)
		RM=true;
		shift
		;;
	-rm-no)
		RM=false;
		shift
		;;
	-run-only)
		GEN="none"
		RM=false
		shift
		;;
	-server)
		SERVER=true
		shift
		;;
	-unobs)
		UNOBSERVED=$2
		shift 2
		;;
	-rep)
		REPEAT_COUNT=$2
		shift 2
		;;
	-mut)
		STEPMUTATIONRATE_MIN=$2
		STEPMUTATIONRATE_MAX=$3
		STEPMUTATIONRATE_INC=$4
		STEPMUTATIONRATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($STEPMUTATIONRATE_MIN, $STEPMUTATIONRATE_MAX, $STEPMUTATIONRATE_INC)]));"`
		shift 4
		;;
	-ozr)
		OZRATE_MIN=$2
		OZRATE_MAX=$3
		OZRATE_INC=$4
		OZRATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($OZRATE_MIN, $OZRATE_MAX, $OZRATE_INC)]));"`
		shift 4
		;;
	-zor)
		ZORATE_MIN=$2
		ZORATE_MAX=$3
		ZORATE_INC=$4
		ZORATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($ZORATE_MIN, $ZORATE_MAX, $ZORATE_INC)]));"`
		shift 4
		;;
	-mvr)
		MISSING_VALUE_RATE_MIN=$2
		MISSING_VALUE_RATE_MAX=$3
		MISSING_VALUE_RATE_INC=$4
		MISSING_VALUE_RATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($MISSING_VALUE_RATE_MIN, $MISSING_VALUE_RATE_MAX, $MISSING_VALUE_RATE_INC)]));"`
		shift 4
		;;
	-locus)
		LOCUS_MIN=$2
		LOCUS_MAX=$3
		LOCUS_INC=$4
		shift 4
		;;
	-sample)
		SAMPLE_MIN=$2
		SAMPLE_MAX=$3
		SAMPLE_INC=$4
		shift 4
		;;
	-clone)
		CLONE_MIN=$2
		CLONE_MAX=$3
		CLONE_INC=$4
		shift 4
		;;
	-step)
		STEP_MIN=$2
		STEP_MAX=$3
		STEP_INC=$4
		shift 4
		;;
	-clone)
		CLONE_MIN=$2
		CLONE_MAX=$3
		CLONE_INC=$4
		shift 4
		;;
	-run)
		RUN=$2
		shift 2
		;;
	-run2)
		RUN2=,$2,
		shift 2
		;;
	-seed)
		SEED=$2
		shift 2
		;;
	-aic)
		AIC=$2
		shift 2
		;;
	-adc)
		ADC=$2
		shift 2
		;;
	-akc)
		AKC=$2
		shift 2
		;;
	-keep)
		KEEPTESTS=true
		shift 1
		;;
	-remove-dup)
		REMOVE_DUP_OPT=" --remove-dup true "
		shift 1
		;;
	-bscfolder)
		BSCFOLDER=$2
		shift 2
		;;
	-deathrate)
		DEATH_RATE=$2
		shift 2
		;;
	-divrate)
		DIV_RATE=$2
		shift 2
		;;
	-obsprob)
		OBS_PROB=$2
		shift 2
		;;
	-samplebound)
		GEN_SC_SAMPLE_BOUND_LB=$2
		GEN_SC_SAMPLE_BOUND_UB=$3
		shift 3
		;;
	-rootheight)
		ROOT_HEIGHT=$2
		shift 2
		;;
	-motionstep)
		MOTION_STEP_LENGTH=$2
		shift 2
		;;
	-data-file)
		# a comma separated list of files
		DATA_FILES=$2
		shift 2
		;;
	-m-distance-norm-pair)
		MATRIX_DISTANCE_NORMALIZATION="-m-distance-norm-pair"
		shift
		;;
	-m-distance-norm-no)
		MATRIX_DISTANCE_NORMALIZATION="-m-distance-norm-no"
		shift
		;;
	-mm-distance-norm)
		MATRIX_MATRIX_DISTANCE_NORMALIZATION=
		shift
		;;
	-mm-distance-no-norm)
		MATRIX_MATRIX_DISTANCE_NORMALIZATION="-no-normal"
		shift
		;;
	-exit-on-error-no)
		EXIT_ON_ERROR=0
		shift
		;;
	-timeout)
		TIMEOUT="timeout $2 "
		TIMEOUT_SECONTS=$3
		if ! [[ $TIMEOUT_SECONTS =~ '^[0-9]+$' ]] ; then
			echo "error: TIMEOUT in SECONDS Not a number" >&2; exit 1
		fi
		shift 3
		;;
	-opt-gen)
		OPT_GEN=$2
		shift 2
		;;
	-cpu-limit)
		RUN_OPT_SICLONEFIT='-opt -XX:ActiveProcessorCount=1'
		RUN_OPT_SIFIT='-opt -XX:ActiveProcessorCount=1'
		echo "--CPUT-LIMIT $RUN_OPT_SIFIT"
		shift
		;;
	--) # End of all options
		shift
		break;
		;;
	-*)
		echo "Error: Unknown option: $1" >&2
		exit 1
		;;
	*)  # No more options
		break
		;;
	esac
done


if [ "$EXIT_ON_ERROR" = "1" ]; then
	# exit when any command fails
	set -e

	# keep track of the last executed command
	trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
	# echo an error message before exiting
	trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
fi

timestamp() { echo -n "timestamp:"; date +"%d/%m/%Y %H:%M:%S"; }

echo "parameters: RM:$RM SERVER:$SERVER STEP_MIN:$STEP_MIN STEP_MAX:$STEP_MAX STEP_INC:$STEP_INC CLONE_MIN:$CLONE_MIN CLONE_MAX:$CLONE_MAX CLONE_INC:$CLONE_INC SAMPLE_MIN:$SAMPLE_MIN SAMPLE_MAX:$SAMPLE_MAX SAMPLE_INC:$SAMPLE_INC LOCUS_MIN:$LOCUS_MIN LOCUS_MAX:$LOCUS_MAX LOCUS_INC:$LOCUS_INC MISSING_VALUE_RATE_MIN:$MISSING_VALUE_RATE_MIN MISSING_VALUE_RATE_MAX:$MISSING_VALUE_RATE_MAX MISSING_VALUE_RATE_INC:$MISSING_VALUE_RATE_INC ZORATE_MIN:$ZORATE_MIN ZORATE_MAX:$ZORATE_MAX ZORATE_INC:$ZORATE_INC OZRATE_MIN:$OZRATE_MIN OZRATE_MAX:$OZRATE_MAX OZRATE_INC:$OZRATE_INC REPEAT_COUNT:$REPEAT_COUNT STEPMUTATIONRATE_MIN:$STEPMUTATIONRATE_MIN STEPMUTATIONRATE_MAX:$STEPMUTATIONRATE_MAX STEPMUTATIONRATE_INC:$STEPMUTATIONRATE_INC SDIR:$SDIR ZORATE_RANGE: $ZORATE_RANGE OZRATE_RANGE: $OZRATE_RANGE MISSING_VALUE_RATE_RANGE: $MISSING_VALUE_RATE_RANGE MISSING_VALUE_RATE_RANGE: $MISSING_VALUE_RATE_RANGE RUN: $RUN RUN2: $RUN2 EXIT_ON_ERROR: $EXIT_ON_ERROR" 


if [ "$SERVER" = true ]; then
	echo "SERVER parameter not supported"
	exit 1
fi

if $RM ; then
	rm -rf $SDIR/*
fi
mkdir -p $SDIR

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
	echo "DIR=$SDIR SAMPLE=$SAMPLE LOCUS=$LOCUS STEP=$STEP MVR=$MISSING_VALUE_RATE 0->1=$ZORATE 1->0:$OZRATE STEPMUT:$STEPMUTATIONRATE CLONE:$CLONE ALLCNT:$ALLCNT DIV_RATE=$DIV_RATE DEATH_RATE=$DEATH_RATE "

	for ((CNT=0; CNT<REPEAT_COUNT; CNT++)); do

		# REVERT IT TO THE NEXT LINE FAST!
		DIR=$SDIR/$ALLCNT/$CNT/
		#DIR=$SDIR/$CNT
		if $KEEPTESTS; then
			DIR=$SDIR/$ALLCNT-$CNT
		fi
		mkdir -p $DIR

		if [ "$GEN" = "onco" ]; then
			$SCRIPT/gen-onco.sh $DIR $SAMPLE $CLONE $LOCUS $((7+77777*CNT+777*ALLCNT+SEED)) $UNOBSERVED $ZORATE $OZRATE $MISSING_VALUE_RATE 2>&1
		elif [ "$GEN" = "synth" ]; then
			$SCRIPT/gen-synthesis.sh $DIR $SAMPLE $LOCUS $STEP $MISSING_VALUE_RATE $ZORATE $OZRATE $STEPMUTATIONRATE $((7+77777*CNT+777*ALLCNT+SEED)) $REMOVE_DUP_OPT $AIC $ADC $AKC "$OPT_GEN" 2>&1
		elif [ "$GEN" = "spsimul" ]; then
			$SCRIPT/gen-spsimul.sh $DIR $LOCUS $DEATH_RATE $DIV_RATE $OBS_PROB $ROOT_HEIGHT $MOTION_STEP_LENGTH $MISSING_VALUE_RATE $ZORATE $OZRATE $STEPMUTATIONRATE $((7+77777*CNT+777*ALLCNT+SEED)) $GEN_SC_SAMPLE_BOUND_LB $GEN_SC_SAMPLE_BOUND_UB $MATRIX_DISTANCE_NORMALIZATION 2>&1
			if [ $? -ne 0 ]; then
				echo "Gen failed"
				continue
			fi
		elif [ "$GEN" = "spsimul-load" ]; then
			$SCRIPT/gen-spsimul-load.sh $DIR $DIRFILES/$CNT 2>&1
		elif [ "$GEN" = "data" ]; then
			cp $DATA_FILE $DIR/input-scite.txt
			python $SCRIPT/convert-input.py $DIR/input-scite.txt $DIR/input-impute.txt $DIR/input-bp.txt 
		elif [ "$GEN" = "none" ]; then
			:
		else
			echo "Gen not a correct method $GEN"
		fi

		CMP_STR="$DIR/true-distance-matrix True$CNT "
		CMP_PART_STR="$DIR/true-part pTrue$CNT "
		OUT_STR="True$CNT $DIR/true-tree.txt $DIR/true-clones.txt"

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
			$TIMEOUT $SCRIPT/run-sifit.sh $DIR 10000 $MATRIX_DISTANCE_NORMALIZATION $RUN_OPT_SIFIT 2>&1
			CMP_STR="$CMP_STR $DIR/sifit-distance-matrix sifit$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/sifit-part psifit$CNT "
			OUT_STR="$OUT_STR SiFit$CNT $DIR/sifit-tree.txt $DIR/sifit-clones.txt"
			echo "RUN_END SIFIT"
		fi

#Running siclonefit
		if [[ $RUN2 =~ .*,siclonefit,.* ]]; then
			timestamp
			$TIMEOUT $SCRIPT/run-siclonefit.sh $DIR 1 $MATRIX_DISTANCE_NORMALIZATION $RUN_OPT_SICLONEFIT 2>&1
			CMP_STR="$CMP_STR $DIR/siclonefit-distance-matrix siclonefit$CNT "
			CMP_PART_STR="$CMP_PART_STR $DIR/siclonefit-part psiclonefit$CNT "
			OUT_STR="$OUT_STR SiCloneFit$CNT $DIR/siclonefit-tree.txt $DIR/siclonefit-clones.txt"
			echo "RUN_END SICLONEFIT"
		fi


# COMPARISON
		echo "STR: $CMP_STR $CMP_PART_STR"
		python $SCRIPT/calc-matrix-distance.py $MATRIX_MATRIX_DISTANCE_NORMALIZATION $CMP_STR
		python $SCRIPT/calc-part-distance.py --match --normalize $CMP_PART_STR
	done


done 
done 
done 
done 
done 
done 
#done 
done 
done 
done 
done

echo "Done TEST-GEN-ALL"
