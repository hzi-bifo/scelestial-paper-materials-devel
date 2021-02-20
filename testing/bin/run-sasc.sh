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


sed 's/2/3/g' < $DIR/input-scite.txt > $DIR/input-sasc.txt

time $SCRIPT/sasc/sasc -l -i $DIR/input-sasc.txt -n `cat $DIR/input-scite.txt | wc -l` -m `head -1 $DIR/input-scite.txt | grep -oE [[:digit:]] | wc -l` -a 0.1 -b 0.1 -k 1 > $DIR/sasc.txt

python $SCRIPT/sasc-2-tree-clone.py $DIR/sasc.txt $DIR/input-scite.txt $DIR/sasc-tree.txt $DIR/sasc-clones.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/sasc-tree.txt $DIR/sasc-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/sasc-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/sasc-tree.txt $DIR/sasc-clones.txt > $DIR/sasc-part
#python visualize-tree.py SCITE $DIR/sasc-tree.txt $DIR/sasc-clones.txt $DIR/output-sifit

