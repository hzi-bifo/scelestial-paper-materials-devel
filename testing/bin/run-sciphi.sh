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


python $SCRIPT/convert-scite-input-2-sciphi.py $DIR/input-scite.txt $DIR/input-sciphi.mpileup $DIR/sciphi-cellnames.txt 1
time $SCRIPT/SCIPhI/build/sciphi -o $DIR/sciphi-out --in $DIR/sciphi-cellnames.txt --seed 4 $DIR/input-sciphi.mpileup

python $SCRIPT/sciphi-2-tree-clone.py $DIR/sciphi-out/best_index/tree.gv $DIR/input-scite.txt $DIR/sciphi-tree.txt $DIR/sciphi-clones.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/sciphi-tree.txt $DIR/sciphi-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/sciphi-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/sciphi-tree.txt $DIR/sciphi-clones.txt > $DIR/sciphi-part
#python visualize-tree.py SCITE $DIR/sciphi-tree.txt $DIR/sciphi-clones.txt $DIR/output-sciphi

