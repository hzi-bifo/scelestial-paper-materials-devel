DIR=$1
MAXK=$2
shift 2

MATRIX_DISTANCE_NORMALIZATION="-no-normal"
while :
do
	case "$1" in
	-m-distance-norm-pair)
		MATRIX_DISTANCE_NORMALIZATION="-normal-mean"
		shift
		;;
	-m-distance-norm-no)
		MATRIX_DISTANCE_NORMALIZATION="-no-normal"
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

time $SCRIPT/scelestial -min-k 3 -max-k $MAXK < $DIR/input-impute.txt > $DIR/impute-tree-clone.txt 
python $SCRIPT/steiner-to-clone-tree.py $DIR/input-scite.txt $DIR/impute-tree-clone.txt $DIR/impute-tree.txt $DIR/impute-clones.txt 2>/dev/null
python $SCRIPT/steiner-to-seq.py $DIR/input-scite.txt $DIR/impute-tree-clone.txt > $DIR/impute-seq.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/impute-tree.txt $DIR/impute-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/impute-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/impute-tree.txt $DIR/impute-clones.txt > $DIR/impute-part
