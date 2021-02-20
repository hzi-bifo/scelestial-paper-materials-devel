DIR=$1
TIMEOUT=$2
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

python $SCRIPT/convert-scite-input-2-phiscs.py $DIR/input-scite.txt $DIR/input-phiscs.txt
time python $SCRIPT/PhISCS/PhISCS-B/csp_z3.py -o $DIR/phiscs -f $DIR/input-phiscs.txt -p 0.1 -n 0.1 -t 1  -m 0 --timeout $TIMEOUT -w 0

python $SCRIPT/phiscs-2-tree-clone.py $DIR/phiscs/input-phiscs.CFMatrix $DIR/input-scite.txt $DIR/phiscs-tree.txt $DIR/phiscs-clones.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/phiscs-tree.txt $DIR/phiscs-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/phiscs-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/phiscs-tree.txt $DIR/phiscs-clones.txt > $DIR/phiscs-part
#python visualize-tree.py SCITE $DIR/phiscs-tree.txt $DIR/phiscs-clones.txt $DIR/output-phiscs

