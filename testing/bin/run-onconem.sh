DIR=$1
shift 1

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

time Rscript $SCRIPT/run-onconem.R $DIR/input-scite.txt `pwd` $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt 
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/onconeme-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt > $DIR/onconeme-part
#python visualize-tree.py Onco $DIR/onconeme-tree.txt $DIR/onconeme-clones.txt $DIR/output-onconeme

