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

mkdir -p $DIR/temp
time $SCRIPT/SCITE/scite -i $DIR/input-scite.txt -n `cat $DIR/input-scite.txt | wc -l` -m `head -1 $DIR/input-scite.txt | sed 's/[^,]//g' | wc -c` -r 1 -l 200000 -fd 0.00005 -ad 0.21545 0.1 -o $DIR/temp/scite- -max_treelist_size 1 -a -seed 7 
python $SCRIPT/scite-2-tree-clone.py $DIR/temp/scite-_ml0.gv $DIR/input-scite.txt $DIR/scite-tree.txt $DIR/scite-clones.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/scite-tree.txt $DIR/scite-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/scite-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/scite-tree.txt $DIR/scite-clones.txt > $DIR/scite-part
#python visualize-tree.py SCITE $DIR/scite-tree.txt $DIR/scite-clones.txt $DIR/output-scite

